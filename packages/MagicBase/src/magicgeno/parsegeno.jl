
function readgenodf(genofile::AbstractString,pedinfo::Union{Integer,AbstractString};
    isfounderinbred::Bool=true,
    formatpriority::AbstractVector=["AD","GT"],
    isdelmultiallelic::Bool=true,
    commentstring::AbstractString="##",
    missingstring=["NA","missing"],
    logfile::Union{Nothing, AbstractString,IO} = nothing,
    workdir::AbstractString=pwd(),
    verbose::Bool=false)
    ext = last(split_allext(genofile))
    ext in [".csv",".csv.gz", ".vcf",".vcf.gz"] || @error string("genofile ext must be .csv, .csv.gz, .vcf or .vcf.gz")
    genofile2=getabsfile(workdir,genofile)
    isfile(genofile2) || @error string(genofile2, " does not exit")
    isvcf = ext in [".vcf",".vcf.gz"]    
    lastcomment = findlastcomment(genofile2; commentstring)
    if isvcf
        lastcomment == 0 && @warn string("no comment lines for vcfile=",genofile,", commentstring=",commentstring)
        lastcomment -= 1 # vcf last comment line: #CHROM, ..., is the headrow of dataframe
    end
    commentlines = read_headlines(genofile2, lastcomment)    
    outstem = tempname(workdir,cleanup=true)
    outext = ".csv.gz"
    tempgenofile = outstem*outext
    try 
        if isvcf            
            parsevcf(genofile2, pedinfo;isfounderinbred, formatpriority,isdelmultiallelic,
                commentstring, outstem,outext, workdir, logfile, verbose
            ) # "NA" as misingstring in output of parsevcf            
        end        
        if isvcf || ext == ".csv.gz"
            genofile3 = isvcf ? tempgenofile : genofile2            
            genodf = readdataframe(genofile3;delim=',',missingstring,commentstring)
        else
            # ext == ".csv"
            genodf = CSV.read(genofile2,DataFrame;delim=',',comment = commentstring,missingstring)
        end
        genodf,commentlines
    finally
        isvcf && rm(tempgenofile; force=true)        
    end
end

function readgenodf(genofile::AbstractString;
    formatpriority::AbstractVector=["AD","GT"],
    isdelmultiallelic::Bool=true,
    commentstring::AbstractString="##",
    missingstring=["NA","missing"],
    logfile::Union{Nothing, AbstractString,IO} = nothing,
    workdir::AbstractString=pwd(),
    verbose::Bool=false)
    nfounder = 0
    readgenodf(genofile,nfounder;
        isfounderinbred=false,formatpriority,isdelmultiallelic,
        commentstring, missingstring, workdir,
        logfile,verbose
    )
end


function parsegenodf!(mapdf::DataFrame, fgenodf::Union{Nothing,DataFrame},
    offgenodf::Union{Nothing,DataFrame};
    isfounderinbred::Bool=true)
    snpsls =splitindex(mapdf[!,:linkagegroup])
    markermap = [mapdf[snps,:] for snps in snpsls]    
    if isnothing(fgenodf)
        foundergeno = nothing
    else
        formatls = mapdf[:,:founderformat]
        foundergeno = [Matrix{Any}(fgenodf[snps,:]) for snps in snpsls]   
        empty!(fgenodf) 
        for chr in eachindex(snpsls)
            chrfgeno = foundergeno[chr]
            chrformatls = formatls[snpsls[chr]]            
            isnonmiss = .!ismissing.(chrfgeno);
            chrfgeno[isnonmiss] .= string.(chrfgeno[isnonmiss])                                    
            parsegenomtx!(chrfgeno,chrformatls)
            isfounderinbred && toinbredgeno!(chrfgeno, chrformatls)
            markermap[chr][!,:founderformat] .= chrformatls
        end                 
    end    
    if isnothing(offgenodf)
        offspringgeno = nothing
    else
        formatls = mapdf[:,:offspringformat]
        offspringgeno = [Matrix{Any}(offgenodf[snps,:]) for snps in snpsls]     
        empty!(offgenodf)   
        GC.gc()
        for chr in eachindex(snpsls)
            chroffgeno = offspringgeno[chr]
            chrformatls = formatls[snpsls[chr]]
            isnonmiss = .!ismissing.(chroffgeno);
            chroffgeno[isnonmiss] .= string.(chroffgeno[isnonmiss])                        
            parsegenomtx!(chroffgeno,chrformatls) # takes lots of time for GP                        
        end        
    end        
    empty!(mapdf)
    markermap, foundergeno, offspringgeno
end

function parsecol_pos_error!(genodf::AbstractDataFrame)
    colls = propertynames(genodf)
    intersect!(colls, [:poscm,:physposbp, :foundererror, :offspringerror, :peroffspringerror, :seqerror, :allelebalancemean, :allelebalancedisperse,:alleledropout])
    for col in colls
        missingstringls = ["","missing","NA"]
        postype = col == :physposbp ? Int : Float64
        if eltype(genodf[!,col]) <: AbstractString    
            genodf[!,col] .= [in(i,missingstringls) ? missing : parse(postype,i) for i in genodf[!,col]]
        end
    end
end

function splitgenodf(genodf::DataFrame,founderid::AbstractVector,
    offspringid::AbstractVector;
    isphysmap::Bool=false,
    recomrate::Real=1.0)
    #  [:marker, :linkagegroup,:poscm, :physposbp, :info,
    #     :founderformat,:offspringformat,:foundererror,:offspringerror,
    #     :seqerror,:allelebalancemean,:allelebalancedisperse,:alleledropout]    
    sampleid = names(genodf)[_col_1stsample:end]
    # println("names(genodf)=",names(genodf), ";founderid=",founderid,",offspringid=",offspringid)    
    diff = setdiff(sampleid,founderid,offspringid)
    if !isempty(diff)
        msg = string("genotyped individuals in genofile but not in pedfile: ",diff)
        @warn msg
    end
    # mapdf
    mapdf = genodf[!,1:_col_1stsample-1]
    parsemarkermap!(mapdf)
    iscmmiss= ismissing.(mapdf[!,:poscm])
    if isphysmap
        b = ismissing.(mapdf[!,:physposbp])
        if all(b)
            msg = string("physical map is missing")
            @warn msg
        else
            if any(b)
                msg = string(sum(b), " out of ", length(iscmmiss), " markers' physical positions are missing")
                @warn msg
            end
            maxloc = maximum(skipmissing(mapdf[!,:physposbp]))
            if maxloc<10^6.0
                @warn string("marker positions in physical map must be in base pair")
            end
            all(iscmmiss) || @warn string("genetic positions are overwritten by transforming physical position using recomrate =",recomrate," cM/Mbp")
            mapdf[!,:poscm] .= [ismissing(i) ? i : i*recomrate*10^(-6.0) for i in mapdf[!,:physposbp]]
            mapdf[!,:linkagegroup] .= mapdf[!,:physchrom]            
        end
    else
        if all(iscmmiss)            
            msg = string("genetic map is missing")
            @info msg
        else
            if any(iscmmiss)
                msg = string(sum(iscmmiss), " out of ", length(iscmmiss), " markers' genetic positions are missing")
                @warn msg
            end
            maxloc = maximum(skipmissing(mapdf[!,:poscm]))
            if maxloc>10^6.0
                @warn string("marker positions in genetic map must be in centiMorgan")
            end
        end        
    end
    coldict = Dict(names(genodf)[_col_1stsample:end] .=> _col_1stsample:size(genodf,2))
    # fgenodf    
    newfounderid = copy(founderid)
    if isempty(founderid)
        fgenodf = nothing
    else
        cols = [get(coldict,i,nothing) for i=founderid]
        missind = founderid[isnothing.(cols)]
        if !isempty(missind)
            msg = string(length(missind), " founders in pedfile but not in genofile: ",missind)
            @warn msg
            b = .!isnothing.(cols)
            cols = cols[b]
            newfounderid = founderid[b]
        end
        fgenodf = genodf[!, cols]
    end
    # offgenodf
    newoffspringid = copy(offspringid)
    if isempty(offspringid)
        offgenodf = nothing
    else
        cols = [get(coldict,i,nothing) for i=offspringid]
        b = isnothing.(cols)        
        if any(b)
            msg = string(sum(b), " offspring in pedfile but not in genofile: ",offspringid[b])
            @warn msg            
            cols = cols[.!b]
            newoffspringid = offspringid[.!b]
        end
        offgenodf = genodf[!, cols]
    end    
    mapdf, fgenodf, newfounderid, offgenodf,newoffspringid
end

function toinbredgeno!(geno::AbstractMatrix,format::AbstractVector)
    b = format .== "GT_phased"
    if sum(b)>0
        nhetero =  sum([i in [["1","2"],["2","1"]] for i=geno[b,:]])
        if nhetero > 0
            @warn string(nhetero, " out of ", length(geno), " heterozygous genotypes are set to missing")
        end
        dict = Dict([["1","1"],["1","2"],["2","1"],["2","2"],
            ["1","N"],["N","1"],["2","N"],["N","2"],["N","N"]] .=>
            ["1","N","N","2","1","1","2","2","N"])
        geno[b,:] = [get(dict,i,i) for i=geno[b,:]]
        format[b] .= "GT_haplo"
    end
    b = format .== "GT_unphased"
    if sum(b)>0
        nhetero =  sum([i in ["12","21"] for i=geno[b,:]])
        if nhetero > 0
            @warn string(nhetero, " out of ", length(geno), " heterozygous genotypes are set to missing")
        end
        dict = Dict(["11","12","21","22","1N","N1","2N","N2","NN"] .=>
            ["1","N","N","2","1","1","2","2","N"])
        geno[b,:] = [get(dict,i,i) for i=geno[b,:]]
        format[b] .= "GT_haplo"
    end
    geno, format
end


function tostringfgl(fgl::Vector{Vector{T}} where T <: AbstractFloat;
    digits::Integer=6)
    y = reduce(hcat,fgl)
    breaks = join(round.(y[1,:],digits=digits),"|")
    labels = string(join(Int.(y[2,1:end-1]),"|"),"|",y[2,end])
    string(breaks,"=>",labels)
end

function parsestringfgl(strfgl::AbstractString)
    # tostringfgl(getfgl(strfgl)) == strfgl
    fgl = [parse.(Float64,i) for i=split.(split(strfgl,"=>"),"|")]
    [[fgl[1][i],fgl[2][i]] for i=1:length(fgl[1])]
end

function parsegenomtx!(genomtx::AbstractMatrix, formatvec::AbstractVector)
    res = genomtx    
    formatset = unique(formatvec)
    if any(ismissing.(formatvec))
        println("there exist unknown genoformat: ",formatset)
    end
    for format = formatset
        ii = findall(formatvec .== format)
        subres =view(res, ii,:)
        if format == "AD"
            ismiss = ismissing.(subres)
            subres[ismiss] .= [Int16[0,0] for i in 1:sum(ismiss)]
            subres[.!ismiss] .= [begin
                ad = tryparse.(Int16,i)
                length(ad)>2 && (ad = [ad[1],sum(ad[2:end])])
                ad
            end for i in split.(subres[.!ismiss],"&")]
        elseif format == "GT_haplo"
            ismiss = ismissing.(subres)
            subres[ismiss] .= "N"
        elseif format == "GT_phased"
            ismiss = ismissing.(subres)
            subres[ismiss] .= [["N","N"] for i in 1:sum(ismiss)]
            subres[.!ismiss] .= [string.(strip.(split(i,"|"))) for i in subres[.!ismiss]]
        elseif format == "GT_unphased"
            ismiss = ismissing.(subres)
            subres[ismiss] .= "NN"
        elseif format == "GT"
            ismiss = ismissing.(subres)
            subres[ismiss] .= "NN"
            subres[.!ismiss] .= [occursin("|",i) ? string.(strip.(split(i,"|"))) : i for i in subres[.!ismiss]]
        elseif format == "GL"
            ismiss = ismissing.(subres)
            subres[ismiss] .= [zeros(Float32,3) for i in 1:sum(ismiss)]
            subres[.!ismiss] .= [tryparse.(Float32,i) for i in split.(subres[.!ismiss],"&")]
        elseif format == "GP"
            # genotype posterior probability (GP)
            # TODO: multi-allles, or phased genotypes, 4-element prob vector
            ismiss = ismissing.(subres)
            isnonmiss = .!ismiss
            subres[ismiss] .= [Float32[0.25,0.5,0.25] for _ in 1:sum(ismiss)]
            subres[isnonmiss] .= [tryparse.(Float32,split(i,"&")) for i in subres[isnonmiss]];            
        elseif format == "discretefgl"
            # no missing for the truevalue resulting from simulation
            # discretefgl must be stored in csv file
            subres .= [tryparse.(Int16,i) for i=split.(subres,"|")]
        elseif format == "contfgl"
            # no missing for the truevalue resulting from simulation            
            # discretefgl must be stored in csv file
            subres .= parsestringfgl.(subres)
        else
            error(string("unknown geno format: ",format))
        end
        # res[ii,:] .= subres
    end
    res
end

function inverse_parsegenomtx!(genomtx::AbstractMatrix, formatvec::AbstractVector,
    isvcf::Bool,missingstring::AbstractString)
    formatset = unique(formatvec)
    res = genomtx
    if isvcf
        for format = formatset
            ii = formatvec .== format
            subres = view(res, ii,:)            
            if format == "GT_haplo"
                rule = ("N"=>".|.","1"=>"0|0","2"=>"1|1")
                uni_subres = unique(subres)       
                d = setdiff(uni_subres, first.(rule))
                if !isempty(d) 
                    error(string("unknown genotypes = ",d, " for format=", format))
                end
                d = setdiff(uni_subres, first.(rule))
                replace!(subres,rule...)
            elseif format == "GT_phased"
                rule = ("N|N"=>".|.","1|1"=>"0|0","1|2"=>"0|1","2|1"=>"1|0","2|2"=>"1|1",
                        "1|N"=>"0|.","N|1"=>".|0","2|N"=>"1|.","N|2"=>".|1")                                
                subres .= join.(subres,"|")
                uni_subres = unique(subres)                
                d = setdiff(uni_subres, first.(rule))
                if !isempty(d) 
                    error(string("unknown genotypes = ",d, " for format=", format))
                end
                replace!(subres,rule...)
            elseif format == "GT_unphased"
                rule = ("NN"=>"./.","11"=>"0/0","12"=>"0/1","21"=>"0/1","22"=>"1/1",
                        "1N"=>"0/.","N1"=>"0/.","2N"=>"1/.","N2"=>"1/.")         
                uni_subres = unique(subres)       
                d = setdiff(uni_subres, first.(rule))
                if !isempty(d) 
                    error(string("unknown genotypes = ",d, " for format=", format))
                end
                replace!(subres,rule...)
            elseif format in ["AD", "GL", "GP"]
                b = ismissing.(subres)
                if any(b)                    
                    subres[b] .= missingstring
                    subres[.!b] .= join.(subres[.!b],",")
                else
                    subres .= join.(subres,",")
                end
            elseif format in ["discretefgl", "contfgl"]
                error(string(format," must be stored in csv file"))
            else
                error(string("unknown geno format: ",format))
            end
        end
    else
        for format = formatset
            if format in ["GT", "GT_phased","AD", "GL", "GP","discretefgl","contfgl"]
                ii = formatvec .== format
                subres = view(res, ii,:)
                if format in ["AD", "GL", "GP"]                    
                    b = ismissing.(subres)
                    if any(b)                    
                        subres[b] .= missingstring
                        subres[.!b] .= join.(subres[.!b],"&")
                    else
                        subres .= join.(subres,"&")
                    end
                elseif format in ["GT_phased","discretefgl"]                    
                    subres .= join.(subres,"|")
                elseif format in ["GT"]                    
                    b = ismissing.(subres)
                    if any(b)                    
                        subres[b] .= missingstring
                        subres[.!b] .= [isa(i,AbstractVector) ? join(i, "|") : i for i in subres[.!b]]
                    else
                        subres .= [isa(i,AbstractVector) ? join(i, "|") : i for i in subres]
                    end
                else
                    #  format == "contfgl"
                    subres .= tostringfgl.(subres;digits=6)
                end
            else
                in(format, ["GT_haplo",  "GT_unphased"]) || error(string("unknown geno format: ",format))
            end
        end
    end
    res
end

function parsemarkermap!(markermap::DataFrame)
    colnames = [:marker, :linkagegroup,:poscm, :physchrom, :physposbp, :info, :founderformat,:offspringformat,
        :foundererror,:offspringerror,:seqerror,:allelebalancemean,:allelebalancedisperse,:alleledropout]
    ncol = size(markermap,2)
    if ncol >=3 
        rename!(markermap,colnames[1:ncol])
    else
        @error string("number of columns = ", ncol, ",< 3")
    end    
    snpls = Vector{Any}(markermap[!,:marker])
    b = ismissing.(snpls)
    if any(b)
        snpls[b].= string.("snp",findall(b))
        markermap[!,:marker] .= string.(snpls)
    end            
    allunique(markermap[!,:marker]) || @error ("marker IDs are not unique")
    for j = [:marker, :linkagegroup]        
        markermap[!,j] .= [ismissing(i) ? i : string(strip(string(i))) for i in markermap[!,j]]
    end
    b = .!ismissing.(markermap[!,:linkagegroup])
    if any(b)
        snpindexlist=splitindex(markermap[b,:linkagegroup])
        chridls = [markermap[first(i),:linkagegroup] for i=snpindexlist]
        if !allunique(chridls)
            @warn ("markers in the same chromosome are not consecutive")
        end
    else
        snpindexlist = [1:size(markermap,1)]
        chridls = markermap[1:1,:linkagegroup]
    end
    if !any(ismissing.(markermap[!,:poscm]))
        poscmls = [isa(i,AbstractString) ? tryparse(Float32,i) : i for i in markermap[!,:poscm]]
        b = isnothing.(poscmls)
        if any(b)
            msg = string("delete ", sum(b), " markers with poscm: ", unique(markermap[b,:poscm]))
            @warn msg                        
            deleteat!(markermap, b)
            markermap[!,:poscm] = Vector{Float32}(poscmls[.!b])
            snpindexlist=splitindex(markermap[!,:linkagegroup])
        end
        for chr in eachindex(snpindexlist)
            ii = snpindexlist[chr]
            A=markermap[ii,:poscm]
            b = A[1:end-1] .<= A[2:end]
            if !all(b)
                posls = findall(.!b)
                union!(posls,posls .+ 1)
                sort!(posls)
                ncol = min(5, size(markermap,2))
                msg = string("Non-descreasing genetic positions in linkagegroup=",chridls[chr], "\n", markermap[ii[posls],1:ncol])
                @error msg
            end
        end        
    end
    markermap
end
