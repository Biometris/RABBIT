"""
    MagicGeno

mutable struct that stores genotypic data

    MagicGeno(magicped,markermap,foundergeno,offspringgeno,deletion,correction)

inner constructor. See also [`formmagicgeno`](@ref).

# Fields

`magicped::MagicPed`: breeding pedigree information. See also [`MagicPed`](@ref).

`markermap::Vector{DataFrame}`: marker map for each chromosome. `markermap[c]`
gives the markermap of chromosome c.

`foundergeno::Vector{Matrix}`: genotypic data in founders. `foundergeno[ch][m,i]`
gives the genotype of founder i at marker m of chromosome ch.

`offspringgeno::Vector{Matrix}`: genotypic data in offspring. `offspringgeno[ch][m,i]`
 gives the genotype of offspring i at marker m of chomosome ch.

`misc::Dict{String, DataFrame}`: contains information such as (1) "deletion" dataframe
for markers that were removed from markermap, (2) "correction" dataframe for
parental error correction.

"""
mutable struct MagicGeno
    magicped::MagicPed
    markermap::Union{Nothing,Vector{DataFrame}}
    # foundergeno: [chr][marker,founder]
    foundergeno::Union{Nothing,Vector{Matrix}}
    # offspringgeno: [chr][marker,offspring]
    offspringgeno::Union{Nothing,Vector{Matrix}}
    misc::Dict{String, DataFrame}
    function MagicGeno(magicped::MagicPed,markermap,foundergeno,offspringgeno,misc)
        nchrls = Vector{Int}()
        isnothing(markermap) || push!(nchrls,length(markermap))
        isnothing(foundergeno) || push!(nchrls,size(foundergeno,1))
        isnothing(offspringgeno) || push!(nchrls,size(offspringgeno,1))
        if length(union(nchrls))>1
            error("inconsistent number of linkage groups ",nchrls)
        end
        nsnpls = Vector{Vector{Int}}()
        isnothing(markermap) || push!(nsnpls,size.(markermap,1))
        isnothing(foundergeno) || push!(nsnpls,size.(foundergeno,1))
        isnothing(offspringgeno) || push!(nsnpls,size.(offspringgeno,1))
        if length(union(nsnpls))>1
            error("inconsistent number of markers within linkage groups ",nsnpls)
        end
        if !isnothing(offspringgeno)
            noffls =  unique(size.(offspringgeno,2))
            if length(noffls)!=1
                error("inconsistent number of offspring ",noffls)
            end
        end
        new(magicped,markermap,foundergeno,offspringgeno,misc)
    end
end


function insert_rabbitdupe!(genodf::AbstractDataFrame, magicped::MagicPed)
    dupels = get_rabbitdupels(magicped)
    namels = names(genodf)
    setdiff!(dupels,namels)
    for i in eachindex(dupels)
        dupe = replace(dupels[i],"_rabbitdupe"=>"")
        col = findfirst(==(dupe),namels)
        if isnothing(col) 
            @error string("For ",dupels[i], ", col of ", dupe, " is missing in genodf")
        else
            insertcols!(genodf,size(genodf,2)+1,dupels[i] => copy(genodf[!,col]))
        end
    end
end

function get_rabbitdupels(magicped::MagicPed)
    offls = magicped.offspringinfo[!,1]
    b = occursin.("_rabbitdupe",offls) 
    offls[b]
end

const _col_1stsample = 15

function formmagicgeno!(genodf::DataFrame,magicped::MagicPed;
    isfounderinbred::Bool=true,
    isphysmap::Bool=false,
    recomrate::Real=1.0)
    insert_rabbitdupe!(genodf,magicped)
    founderid = isnothing(magicped.founderinfo) ? nothing : magicped.founderinfo[!,:individual]
    offspringid = isnothing(magicped.offspringinfo) ? nothing : magicped.offspringinfo[!,:individual]    
    if isnothing(founderid) && isnothing(offspringid)
        msg = string("no founders and offspring")
        @error msg
    elseif !isnothing(founderid) && !isnothing(offspringid)
        indls = names(genodf)[_col_1stsample:end]
        d= setdiff(indls,founderid,offspringid)
        if !isempty(d)
            msg = string(length(d), " individuals in genofile by not in pedfile: ",d)
            @warn msg
        end
    else
        indls = names(genodf)[_col_1stsample:end]
        isnothing(founderid) && (founderid = setdiff(indls, offspringid))
        isnothing(offspringid) && (offspringid = setdiff(indls, founderid))
    end    
    # fgenodf consistent with newfounderid; offgenodif consistent with newoffspringid    
    mapdf, fgenodf, newfounderid, offgenodf,newoffspringid = splitgenodf(genodf,
        founderid,offspringid; isphysmap,recomrate)            
    if newoffspringid != offspringid
        issubset(newoffspringid,offspringid) || @error "unexpected newoffspringid"
        dict = Dict(offspringid .=> 1:length(offspringid))
        indices = [dict[i] for i in newoffspringid]
        magicped.offspringinfo =  magicped.offspringinfo[indices,:]
    end
    if newfounderid != founderid
        issubset(newfounderid,founderid) || @error "unexpected newfounderid"
        dict = Dict(founderid .=> 1:length(founderid))
        indices = [dict[i] for i in newfounderid]
        magicped.founderinfo =  magicped.founderinfo[indices,:]
    end
    ped = magicped.designinfo
    if isa(ped,Pedigree)
        memls = Vector{String}(unique(magicped.offspringinfo[!,:member]))
        subped = Pedigrees.getsubped(ped,memls)
        magicped.designinfo = subped
        # to keep consistent among subped's founders, magicped.founderinfo, and fgenodf
        nf = subped.nfounder
        pedfounders = subped.member[1:nf]
        if pedfounders != newfounderid
            d = setdiff(pedfounders, newfounderid)
            if isempty(d)
                dict = Dict(newfounderid .=> 1:length(newfounderid))
                indices = [dict[i] for i in pedfounders]
                magicped.founderinfo =  magicped.founderinfo[indices,:]
                fgenodf = fgenodf[:, indices]
            else
                @error string("founders in pedfile but not in genofile: ",d)
            end
        end
    end            
    markermap, foundergeno, offspringgeno = parsegenodf!(mapdf, fgenodf,
        offgenodf;isfounderinbred)    
    # mapdf, fgenodf, offgenodf pointing genodf; empty!(genodf) will empty! fgeno and offgenodf    
    # empty!(genodf)    
    misc = Dict{String, DataFrame}()    
    magicgeno = MagicGeno(magicped,markermap,foundergeno,offspringgeno,misc)
    isphysmap && sortmarker!(magicgeno; byphysmap=true)
    magicgeno
end

"""
    formmagicgeno(genofile, nfounder;
        isfounderinbred, formatpriority, isphysmap, recomrate, commentstring)

form magicgeno::MagicGeno from the `genofile` and the number `nfounder` of founders.
Assume that founders' columns are on the left of offspring columns.

"""
function formmagicgeno(genofile::AbstractString,nfounder::Integer;
    isfounderinbred::Bool=true,
    formatpriority::AbstractVector=["AD","GT"],
    isphysmap::Bool=false,
    recomrate::Real=1.0,    
    delmultiallelic::Bool=true,
    commentstring::AbstractString="##",
    missingstring=["NA","missing"],
    workdir::AbstractString=pwd())
    designinfo = parsedesign(string("nfounder=",nfounder))
    formmagicgeno(genofile,designinfo;
        isfounderinbred, formatpriority,isphysmap,recomrate,
        delmultiallelic,
        commentstring,missingstring,workdir)
end

function formmagicgeno(genofile::AbstractString,designinfo::DesignInfo;
    isfounderinbred::Bool=true,
    formatpriority::AbstractVector=["AD","GT"],    
    isphysmap::Bool=false,
    recomrate::Real=1.0,
    delmultiallelic::Bool=true,
    commentstring::AbstractString="##",
    missingstring=["NA","missing"],
    workdir::AbstractString=pwd())        
    # read genofile
    genodf,commentlines = readgenodf(genofile,designinfo.designcode; isfounderinbred,formatpriority,
        delmultiallelic,commentstring, missingstring, workdir)
    # set magicped    
    indidls = names(genodf)[_col_1stsample:end]
    magicped = formmagicped!(designinfo,indidls)    
    # formmagicgeno
    magicgeno = formmagicgeno!(genodf, magicped;isfounderinbred, isphysmap,recomrate)
    if !isempty(commentlines)
        filecomment = DataFrame(sourcefile=getabsfile(workdir,genofile),
            commentlines=commentlines)
        pushmisc!(magicgeno,"filecomment"=>filecomment)
    end
    magicgeno
end

"""
    formmagicgeno(genofile, pedinfo;
        formatpriority, isphysmap, recomrate, commentstring)

form magicgeno from the `genofile` and the pedigree information `pedinfo`.

# Positional arguments

`genofile::AbstractString`: genotypic data file with extension ".vcf" or ".vcf.gz".

`pedinfo::AbstractString`: `designcode` or `pedfile`

* `pedfile`: See [`readmagicped`](@ref) for the format of a pedigree file.

* `designcode`: a string designcode for a breeding population. 
  See [`parsedesign`](@ref) for details.

# Keyword arguments

`isfounderinbred::Bool=true`: if true, founders are inbred, and otherwise outbred.

`formatpriority::AbstractVector=["AD","GT"]`: the priority of genotype 
formats when parasing input vcf genofile.  

`isphysmap::Bool=false`: if ture, transform physical map into genetic map using recomrate and overwrite the exist genetic map. 
If false, keep input physical and/or genetic map."

`recomrate::Real=1.0`, average recombation rate in cM per Mbp. Valid only if
isphysmap = true.

`commentstring::AbstractString="##"`: the lines beginning with commentstring are ignored
in genofile or pedfile.

`workdir::AbstractString=pwd()`: directory for reading genofile and pedfile.

"""
function formmagicgeno(genofile::AbstractString,pedinfo::AbstractString;
    isfounderinbred::Bool=true,
    formatpriority::AbstractVector=["AD","GT"],
    isphysmap::Bool=false,
    recomrate::Real=1.0,
    delmultiallelic::Bool=true,
    commentstring::AbstractString="##",
    missingstring=["NA","missing"],
    workdir::AbstractString=pwd())
    if last(splitext(pedinfo))==".csv"
        # pedinfo: pedfile
        magicped = readmagicped(pedinfo; commentstring, workdir)
        # read genofile
        genodf,commentlines = readgenodf(genofile,pedinfo; isfounderinbred,formatpriority,
            delmultiallelic,commentstring, missingstring, workdir)
        #form magicgeno
        magicgeno = formmagicgeno!(genodf, magicped;isfounderinbred, isphysmap,recomrate)
        if !isempty(commentlines)
            filecomment = DataFrame(sourcefile=getabsfile(workdir,genofile),
                commentlines=commentlines)
            pushmisc!(magicgeno,"filecomment"=>filecomment)
        end
        magicgeno
    else        
        # pedinfo: designcode            
        designinfo = parsedesign(pedinfo)
        magicgeno = formmagicgeno(genofile,designinfo; 
            isfounderinbred, formatpriority, isphysmap,recomrate,
            delmultiallelic,commentstring, missingstring, workdir)        
    end
end

function readmagicgeno(magicgenofile::AbstractString;
    delim::AbstractChar=',',
    commentstring::AbstractString="##",
    workdir::AbstractString=pwd())
    magicgenofile2=getabsfile(workdir,magicgenofile)    
    # issparse = true, CSV.read is used to parse string, when cells with too long string will be truncated
    res = readmultitable(magicgenofile2; delim,commentstring,isparse=false)    
    haskey(res, "genodata") || error(string("key genodata is not in ", keys(res)))
    genodf = res["genodata"]    
    parsecol_pos_error!(genodf)        
    magicped = parsemagicped(res)    
    formatls = unique(genodf[!,:founderformat])
    isfounderinbred = isempty(intersect(["GT_unphased","GT_phased"], formatls))
    if !isfounderinbred && in("GT_haplo",formatls)
        @warn string("markers with genotypes and haplotypes. founderformat=",formatls)
    end
    magicgeno = formmagicgeno!(genodf,magicped; isfounderinbred,isphysmap=false)
    misc = Dict{String, DataFrame}()
    for (key,val) in res
        key[1:5]=="misc|" && push!(misc, key[6:end]=>val)
    end
    magicgeno.misc = misc
    magicgeno
end

function pushmisc!(magicgeno::MagicGeno,miscpair::Pair{String,DataFrame})
    strkey, df = miscpair
    if haskey(magicgeno.misc,strkey)
        magicgeno.misc[strkey] = vcat(magicgeno.misc[strkey],df)
    else
        push!(magicgeno.misc, miscpair)
    end
end

function formfhaplo(fhaplofile::AbstractString;
    nfounder::Union{Nothing,Integer}=nothing,
    isfounderinbred::Bool=true,
    formatpriority::AbstractVector=["AD","GT"],
    isphysmap::Bool=false,
    recomrate::Real=1.0, #
    commentstring::AbstractString="##",
    missingstring=["NA","missing"],
    workdir::AbstractString=pwd())
    ext = last(split_allext(fhaplofile))
    ext in [".csv",".csv.gz",".vcf",".vcf.gz"] || @error string("fhaplofile ext must be in [.csv,.csv.gz,.vcf,.vcf.gz]")
    myopen = ext in [".vcf.gz",".csv.gz"] ? GZip.open : open
    isvcf = ext in [".vcf",".vcf.gz"]
    fhaplofile2=getabsfile(workdir,fhaplofile)
    lastcomment = findlastcomment(fhaplofile2; commentstring)
    nheader = lastcomment + 1
    if isnothing(nfounder)
        nfounder = 0
        myopen(fhaplofile2,"r") do io
            for _ in 1:nheader-1
                readline(io,keep=true)
            end
            if isvcf
                titlerow = split(readline(io,keep=false),"\t")
                nfounder = length(titlerow) - 9
            else
                readline(io,keep=true)
                titlerow = split(readline(io,keep=false),",")
                nfounder = length(titlerow) - _col_1stsample + 1
            end
        end
    end
    fhaplo = formmagicgeno(fhaplofile, nfounder;
        isfounderinbred,formatpriority, isphysmap, recomrate,
        missingstring, commentstring,workdir)
    fformat = unique(reduce(vcat,fhaplo.markermap)[!,:founderformat])
    if isfounderinbred
        issetequal(fformat,["GT_haplo"]) || @error string("unexpected founderformat: ",
            fformat, ", fouderformat must be GT_haplo for inbred founders!")
    else
        in("GT_haplo",fformat) && @error string("unexpected founderformat: ",
            fformat, ", fouderformat GT_haplo is not allowed for outbred founders!")
    end
    fhaplo
end

function setmagicped!(magicgeno::MagicGeno,magicped::MagicPed)
    # founders
    sample_fid = magicgeno.magicped.founderinfo[!,:individual]    
    founderid = magicped.founderinfo[!,:individual]    
    diff = setdiff(founderid, sample_fid)
    if !isempty(diff)
        msg = string("founders in magicped but not in magicgeno: ",dff)
        @error msg
    end
    diff = setdiff(sample_fid, founderid)
    if !isempty(diff)
        msg = string("founders in magicgeno but not in magicped: ",dff)
        @warn msg
    end    
    sampledict = Dict(sample_fid .=> 1:length(sample_fid))
    indices = [get(sampledict, i, nothing) for i=founderid]
    if indices != 1:length(founderid)
        for i in magicgeno.foundergeno
            i = i[:,indices] 
        end
    end    
    # offspring
    sample_offid = magicgeno.magicped.offspringinfo[!,:individual]
    offid = magicped.offspringinfo[!,:individual]
    diff = setdiff(offid, sample_offid)
    if !isempty(diff)
        msg = string("offspring in magicped but not in magicgeno: ",dff)
        @error msg
    end
    diff = setdiff(sample_offid,offid)
    if !isempty(diff)
        msg = string("offspring in magicgeno but not in magicped: ",dff)
        @warn msg
    end    
    sampledict = Dict(sample_offid .=> 1:length(sample_offid))
    indices = [get(sampledict, i, nothing) for i=offid]
    if indices != 1:length(offid)
        for i in magicgeno.offspringgeno
            i = i[:,indices] 
        end
    end
    magicgeno.magicped = magicped
    magicgeno
end


"""
    setmagicped!(magicgeno::MagicGeno,juncdist::JuncDist)

set magicgeno.designinfo from `juncdist`.

"""
function setmagicped!(magicgeno::MagicGeno,juncdist::JuncDist)
    # checkjuncdist(juncdist) # TODO
    # TODO: reset foundergeno
    juncdist.nfounder == size(magicgeno.magicped.founderinfo,1) || @error "inconsitent nfounder"
    noff=size(magicgeno.magicped.offspringinfo,1)
    magicgeno.magicped.offspringinfo[:member] = repeat(["juncdist"],noff)
    magicgeno.magicped.offspringinfo[:isfglexch] = trues(noff)
    magicgeno.magicped.offspringinfo[:gender] = repeat(["notapplicable"],noff)
    magicgeno.magicped.designinfo = juncdist
    magicgeno
end

"""
    savegenodata(outfile,magicgeno; missingstring="NA",workdir=pwd(),delim=',')

save genotypic data of magicgeno::MagicGeno into outfile

# Positional arguments

`outfile::AbstractString`: output genofile for saving genotypic data.

`magicgeno::MagicGeno`: a struct returned by [`formmagicgeno`](@ref).

# Keyward arguments

`workdir::AbstractString = pwd()`: directory for writing outfile.

`delim::AbstractChar=','`: delimitor character

"""
function savegenodata(sink::Union{IO,AbstractString},magicgeno::MagicGeno;
    target::AbstractString="all",
    keepcomment::Bool=true,
    workdir::AbstractString = pwd(),
    delim::AbstractChar = ',',
    commentstring::AbstractString="##",
    missingstring::AbstractString="NA")
    
    in(target, ["all","founder","offspring"]) || @error string("unknown target=",
        target, ", target must be in [all, founder, offspring]")
    delim2 = delim
    if typeof(sink) <: AbstractString
        _, ext = split_allext(sink)
        myopen = ext in [".vcf.gz",".csv.gz"] ? GZip.open : open
        if ext in [".vcf",".vcf.gz"]
            delim2 = '\t'
        elseif ext in [".csv",".csv.gz"]
            delim2 = ','
        else
            @error string("outfile extension = ", ext, ", not .vcf, .vcf.gz, .csv, or .csv.gz")
        end
        outputfile2 = getabsfile(workdir,sink)
        io = myopen(outputfile2, "w")
    elseif typeof(sink) <: IO
        io = sink
    else
        error(string("unknow sink: ",sink))
    end
    isvcf = delim2 == '\t'
    if isvcf
        if typeof(sink) <: AbstractString
            if haskey(magicgeno.misc,"filecomment")
                commentlines = join(magicgeno.misc["filecomment"][!,:commentlines],"\n")
            else
                commentlines = ""
            end
            cc = commentstring # comment string = "##" for vcf genofile
            
            occursin("fileformat=",commentlines) || write(io, cc*"fileformat=VCFv4.3\n")
            filedate = string(cc*"filedate=",string(now()))
            write(io, filedate,"\n")
            write(io, cc*"source=RABBIT\n")
            msg = cc*"FORMAT=<ID=GP,Number=.,Type=Integer,Description=\"Probabilities of bi-allelic genotypes [0/0,1/1] (or [0/0,0/1,1/1] or [0/0,0/1,1/0,1/1]) for the vector length of 2 (or 3 or 4)\">"
            occursin("<ID=GP",commentlines) || write(io, msg,"\n")
            msg = cc*"FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">"
            occursin("<ID=AD",commentlines) || write(io, msg,"\n")
            msg = cc*"FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
            occursin("<ID=GT",commentlines) || write(io, msg,"\n")
            msg = cc*"INFO=<ID=LINKAGEGROUP,Number=1,Type=String,Description=\"Linkage group in genetic map\">"
            occursin("<ID=LINKAGEGROUP",commentlines) || write(io, msg,"\n")
            msg = cc*"INFO=<ID=POSCM,Number=1,Type=Float,Description=\"Genetic marker position in centiMorgan\">"
            occursin("<ID=POSCM",commentlines) || write(io, msg,"\n")
            msg = cc*"INFO=<ID=FOUNDERERROR,Number=1,Type=Float,Description=\"Founder allelic error rate\">"
            occursin("<ID=FOUNDERERROR",commentlines) || write(io, msg,"\n")
            msg = cc*"INFO=<ID=FOUNDERERROR,Number=1,Type=Float,Description=\"Offspring allelic error rate\">"
            occursin("<ID=FOUNDERERROR",commentlines) || write(io, msg,"\n")
            msg = cc*"INFO=<ID=SEQERROR,Number=1,Type=Float,Description=\"sequencing base error rate\">"
            occursin("<ID=SEQERROR",commentlines) || write(io, msg,"\n")
            msg = cc*"INFO=<ID=SEQERROR,Number=1,Type=Float,Description=\"sequencing allele balance at each marker ~ Beta(alpha,beta) where allelebalancemean = alpha/(alpha+beta) and allelebalancedisperse = 1/(alpha+beta)\">"
            occursin("<ID=SEQERROR",commentlines) || write(io, msg,"\n")
            msg = cc*"INFO=<ID=ALLELEBALANCEDISPERSE,Number=1,Type=Float,Description=\"sequencing allele balance at each marker ~ Beta(alpha,beta) where allelebalancemean = alpha/(alpha+beta) and allelebalancedisperse = 1/(alpha+beta)\">"
            occursin("<ID=ALLELEBALANCEDISPERSE",commentlines) || write(io, msg,"\n")
            msg = cc*"INFO=<ID=ALLELEBALANCEDISPERSE,Number=1,Type=Float,Description=\"probability of one of alleles being dropout for a heterzygous genotype at each marker\">"
            occursin("<ID=ALLELEBALANCEDISPERSE",commentlines) || write(io, msg,"\n")
        end
    else
        if typeof(sink) <: AbstractString
            write(io, commentstring*"fileformat=csv\n")
            filedate = string(commentstring*"fileDate=",replace(string(Date(now())),"-"=>""))
            write(io, filedate,"\n")
            write(io, commentstring*"source=RABBIT\n")
        end
    end
    if keepcomment
        if haskey(magicgeno.misc,"filecomment")
            commentlines = join(magicgeno.misc["filecomment"][!,:commentlines],"\n")
            write(io,commentlines,"\n")
        end
    end
    colnames = geno_colnames(magicgeno,isvcf,target)
    write(io, join(colnames,delim2),"\n")
    nchr = length(magicgeno.markermap)
    for chr in 1:nchr
        mtx = togenomtx(magicgeno, chr; isvcf, missingstring, target)
        writedlm(io,mtx,delim2)
        flush(io)
    end
    typeof(sink) <: AbstractString && close(io)
    sink
end

function geno_colnames(magicgeno::MagicGeno,isvcf::Bool,target::AbstractString)
    if isvcf
        colnames = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL","FILTER", "INFO", "FORMAT"]
    else
        colnames = ["marker","linkagegroup","poscm","physchrom", "physposbp","info", "founderformat","offspringformat",
            "foundererror","offspringerror","seqerror","allelebalancemean","allelebalancedisperse","alleledropout"]
    end
    if !isnothing(magicgeno.foundergeno) && in(target,["all","founder"])
        founderid = magicgeno.magicped.founderinfo[!,:individual]
        append!(colnames, founderid)
    end
    if !isnothing(magicgeno.offspringgeno) && in(target,["all","offspring"])
        offid = magicgeno.magicped.offspringinfo[!,:individual]
        append!(colnames, offid)
    end
    colnames
end

function togenomtx(magicgeno::MagicGeno,chr::Integer;
    isvcf::Bool,
    missingstring::AbstractString="NA",
    target::AbstractString="all")
    in(target, ["all","founder","offspring"]) || @error string("unknown target=",
        target, ", target must be in [all, founder, offspring]")
    if in(target, ["all","founder"])
        if isnothing(magicgeno.foundergeno)
            foundermtx = nothing
        else
            foundermtx = Matrix{Any}(magicgeno.foundergeno[chr])
            formatvec = magicgeno.markermap[chr][!,:founderformat]
            inverse_parsegenomtx!(foundermtx,formatvec,isvcf,missingstring)
        end
    else
        foundermtx = nothing
    end
    if in(target, ["all","offspring"])
        if isnothing(magicgeno.offspringgeno)
            offmtx = nothing
        else
            offmtx = Matrix{Any}( magicgeno.offspringgeno[chr])
            formatvec = magicgeno.markermap[chr][!,:offspringformat]
            inverse_parsegenomtx!(offmtx,formatvec,isvcf,missingstring)
        end
    else
        offmtx = nothing
    end
    #  [:marker, :linkagegroup,:poscm, :physchrom, :physposbp, :info,:founderformat,:offspringformat,
    #   :foundererror,:offspringerror, :seqerror,:allelebalancemean,:allelebalancedisperse,:alleledropout]
    mapmtx = Matrix{Any}(magicgeno.markermap[chr]);   
    submapmtx = view(mapmtx,:, vcat([3],9:_col_1stsample-1))    
    submapmtx .= [ismissing(i) ? i : round(i,digits=8) for i in submapmtx]
    if isvcf
        mapmtx =  tovcfgenomtx!(mapmtx, foundermtx, offmtx)
    else
        submapmtx = view(mapmtx,:, vcat(1:6,9:_col_1stsample-1))
        submapmtx[ismissing.(submapmtx)] .= missingstring
        if target == "founder"
            # col8=:offspringformat, col10=:offspringerror            
            mapmtx[:,8] .= missingstring
            mapmtx[:,10] .= missingstring
        elseif target == "offspring"
            # col7=:founderformat, col9=:foundererror
            mapmtx[:,7] .= missingstring
            mapmtx[:,9] .= missingstring
        end
    end
    b = [isnothing(foundermtx),isnothing(offmtx)]
    if b == [true,true]
        mapmtx
    elseif b == [true,false]
        hcat(mapmtx,offmtx)
    elseif b == [false,true]
        hcat(mapmtx,foundermtx)
    else
        hcat(mapmtx, foundermtx, offmtx)
    end
end

# modify foundermtx and offspingmtx, return new vcfmapmtx
function tovcfgenomtx!(mapmtx::AbstractMatrix,
    foundermtx::Union{Nothing,AbstractMatrix},
    offspringmtx::Union{Nothing,AbstractMatrix})    
    # input mapmtx columns: 
    #  ["marker","linkagegroup","poscm","physchrom", "physposbp","info", "founderformat","offspringformat", 
    #   "foundererror","offspringerror", "seqerror","allelebalancemean","allelebalancedisperse","alleledropout"]
    # vcfmap colnames = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL","FILTER", "INFO", "FORMAT"]
    # mid: col4-col9
    nsnp = size(mapmtx,1)
    vcfmapmtx = Matrix(undef, nsnp, 9)
    vcfmapmtx[:,1:3] .= mapmtx[:,[4,5,1]]
    submapmtx = view(vcfmapmtx,:, 1:3)
    submapmtx[ismissing.(submapmtx)] .= "."
    isnofounder = isnothing(foundermtx)
    isnooff = isnothing(offspringmtx)
    for snp in 1: nsnp
        # set vcf 4-7 columns: "REF","ALT","QUAL","FILTER","INFO"
        info = mapmtx[snp,6]
        vcfmapmtx[snp,4:7] .= "."
        if ismissing(info)            
            vcfinfo = ""
        else
            infols = split.(split(info,";"),"=")
            otherinfo = infols[length.(infols) .== 1]
            vcfinfo = isempty(otherinfo) ? "" : join(reduce(vcat,otherinfo),";")
            infols = infols[length.(infols) .== 2]
            infols = [strip.(i) for i in infols]  
            infodict = OrderedDict([Pair(i...) for i in infols])
            keyls = ["REF","ALT","QUAL","FILTER"]
            for i in eachindex(keyls)
                if haskey(infodict,keyls[i])
                    vcfmapmtx[snp,3+i]= replace(infodict[keyls[i]],"|"=>",")
                    delete!(infodict,keyls[i])
                end
            end
        end        
        linkagegroup = mapmtx[snp,2]
        if !ismissing(linkagegroup)
            vcfinfo *= string(";LINKAGEGROUP=",linkagegroup)
        end
        poscm = mapmtx[snp,3]
        if !ismissing(poscm)
            vcfinfo *= string(";POSCM=",poscm)
        end
        if !ismissing(mapmtx[snp,9])
            vcfinfo *= string(";FOUNDERERROR=",mapmtx[snp, 9])
        end
        if !ismissing(mapmtx[snp,10])
            vcfinfo *= string(";OFFSPRINGERROR=",mapmtx[snp, 10])
        end
        if !ismissing(mapmtx[snp,11])
            vcfinfo *= string(";SEQERROR=",mapmtx[snp, 11])
        end
        if !ismissing(mapmtx[snp,12])
            vcfinfo *= string(";ALLELEBALANCEMEAN=",mapmtx[snp, 12])
        end
        if !ismissing(mapmtx[snp,13])
            vcfinfo *= string(";ALLELEBALANCEDISPERSE=",mapmtx[snp, 13])
        end
        if !ismissing(mapmtx[snp,14])
            vcfinfo *= string(";ALLELEDROPOUT=",mapmtx[snp, 14])
        end
        if !ismissing(info) && !isempty(infodict)
            vcfinfo *=";"*join(map(i->string(i[1],"=",replace(i[2],"|"=>",")),collect(infodict)),";")
        end
        vcfinfo = strip(vcfinfo,';')
        isempty(vcfinfo) && (vcfinfo = ".")
        vcfmapmtx[snp,8] = vcfinfo
        # set vcf 9th column: "FORMAT"
        genoformat = [ismissing(i) ? i : i[1:2] for i in mapmtx[snp,7:8]]
        if isnofounder && isnooff
            vcfmapmtx[snp,9] = "."
        elseif isnofounder
            vcfmapmtx[snp,9] = genoformat[2]
        elseif isnooff
            vcfmapmtx[snp,9] = genoformat[1]
        else
            any(ismissing.(genoformat)) && @error string("missing genoformat: ",genoformat)
            if genoformat[1] == genoformat[2]
                vcfmapmtx[snp,9] =  genoformat[1]
            else
                if genoformat[2] == "GT"
                    foundermtx[snp,:] .= [".:"*i for i=foundermtx[snp,:]]
                    offspringmtx[snp,:] .= [i*":." for i=offspringmtx[snp,:]]
                    vcfmapmtx[snp,9] = join(reverse(genoformat),":")
                else
                    foundermtx[snp,:] .= [i*":." for i=foundermtx[snp,:]]
                    offspringmtx[snp,:] .= [".:"*i for i=offspringmtx[snp,:]]
                    vcfmapmtx[snp,9] = join(genoformat,":")
                end
            end
        end
    end
    # VCFv4.3 requires that REF must be one of A, T, G, C, N
    vcfmapmtx[vcfmapmtx[:,4] .== ".",4] .= "N" 
    vcfmapmtx
end


"""
    savemagicgeno(outfile,magicgeno,workdir=pwd(),delim=',')

save genotypic data of magicgeno::MagicGeno into outfile

# Positional arguments

`outfile::AbstractString`: filename for saving results.

`magicgeno::MagicGeno`: a struct returned by [`formmagicgeno`](@ref).

# Keyward arguments

`missingstring::AbstractString = "NA"`: string for missing values
`workdir::AbstractString = pwd()`: directory for writing outfile.

`delim::AbstractChar=','`: delimitor character

"""
function savemagicgeno(sink::Union{IO,AbstractString},magicgeno::MagicGeno;
    missingstring::AbstractString = "NA",
    workdir::AbstractString = pwd(),
    delim::AbstractChar = ',')
    if typeof(sink) <: AbstractString
        ext = last(splitext(sink))
        myopen = ext == ".gz" ? GZip.open : open
        outputfile2 = getabsfile(workdir,sink)
        io = myopen(outputfile2, "w")
    elseif typeof(sink) <: IO
        io = sink
    else
        error(string("unknow sink: ",sink))
    end
    savemagicped(io,magicgeno.magicped; workdir,delim)
    write(io,"RABBIT",delim,"genodata","\n")
    savegenodata(io,magicgeno; keepcomment=false,missingstring, workdir,delim)
    misc = magicgeno.misc
    if !isempty(misc)
        for (key,val) in misc
            key == "filecomment" && continue
            appenddf(io, val; delim, initial="RABBIT",dfname=string("misc|",key))
        end
    end
    flush(io)
    typeof(sink) <: AbstractString && close(io)
    sink
end

function submagicgeno!(magicgeno::MagicGeno;
    chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    progeny_subset::Union{Nothing,AbstractRange,AbstractVector}=nothing)
    # default nothing === full set
    magicped = magicgeno.magicped
    noff = size(magicped.offspringinfo,1)
    if isnothing(progeny_subset)
        offii = 1:noff
    else
        offii = progeny_subset[progeny_subset .<= noff]
        magicped.offspringinfo = magicped.offspringinfo[offii, :]
    end
    markermap = magicgeno.markermap
    nchr=length(markermap)
    nsnp = maximum(size.(markermap,1))
    chrsubset2 = isnothing(chrsubset) ? (1:nchr) : chrsubset[chrsubset.<= nchr]    
    snpsubset2 = isnothing(snpsubset) ? (1:nsnp) : snpsubset
    if isempty(chrsubset2)
        msg = string("chrsubset=",chrsubset, " becomes empty after excluding index > #chromsomes(=",nchr,")")
        error(msg)
    end
    snpsubsetlist = [snpsubset2[snpsubset2 .<= size(markermap[ch],1)] for ch=chrsubset2]    
    magicgeno.markermap= map((x,y)->markermap[x][y,:],chrsubset2,snpsubsetlist)
    magicgeno.foundergeno=map((x,y)->magicgeno.foundergeno[x][y,:],chrsubset2,snpsubsetlist)
    magicgeno.offspringgeno=map((x,y)->magicgeno.offspringgeno[x][y,offii],chrsubset2,snpsubsetlist)
    magicgeno
end

function submagicgeno(magicgeno::MagicGeno;
    chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    progeny_subset::Union{Nothing,AbstractRange,AbstractVector}=nothing)
    # default nothing === full set
    magicped = deepcopy(magicgeno.magicped)
    noff = size(magicped.offspringinfo,1)
    if isnothing(progeny_subset)
        offii = 1:noff
    else
        offii = progeny_subset[progeny_subset .<= noff]
        magicped.offspringinfo = magicped.offspringinfo[offii, :]
    end
    markermap = magicgeno.markermap
    nchr=length(markermap)
    nsnp = maximum(size.(magicgeno.markermap,1))
    chrsubset2 = isnothing(chrsubset) ? (1:nchr) : chrsubset[chrsubset.<= nchr]
    snpsubset2 = isnothing(snpsubset) ? (1:nsnp) : snpsubset
    if isempty(chrsubset2)
        msg = string("chrsubset=",chrsubset, " becomes empty after excluding index > #chromsomes(=",nchr,")")
        error(msg)
    end
    snpsubsetlist = [snpsubset2[snpsubset2 .<= size(markermap[ch],1)] for ch=chrsubset2]
    markermap= map((x,y)->markermap[x][y,:],chrsubset2,snpsubsetlist)
    foundergeno=map((x,y)->magicgeno.foundergeno[x][y,:],chrsubset2,snpsubsetlist)
    offspringgeno=map((x,y)->magicgeno.offspringgeno[x][y,offii],chrsubset2,snpsubsetlist)
    MagicGeno(magicped,markermap,foundergeno,
        offspringgeno,deepcopy(magicgeno.misc))
end

# function getoffsex(chroffgeno::AbstractMatrix)
#     gset = kron(["1|","2|","N|"],["1","2","N","Y"])
#     push!(gset,kron(["1","2","N"],["1","2","N","Y"])...)
#     if issubset(unique(chroffgeno),gset)
#         malegset = ["1Y","2Y","NY","1|Y","2|Y","N|Y"]
#         offsex = [issubset(unique(i),malegset) ? "male" : "female" for i=eachcol(chroffgeno)]
#     else
#         offsex = repeat(["notapplicable"],size(chroffgeno,2))
#     end
#     offsex
# end


function setunphasedgeno!(magicgeno::MagicGeno)
    founderformat = unique(reduce(vcat,[unique(i[!,:founderformat]) for i=magicgeno.markermap]))
    if !issubset(founderformat,["GT_haplo", "GT_unphased","GT_phased","GP","AD"])
        error(string("founderformat=", founderformat, ", not in [GT_haplo, GT_unphased, GT_phased, GP, AD]"))
    end
    if in("GT_unphased",founderformat)
        for chr in 1:length(magicgeno.foundergeno)
            fgeno = magicgeno.foundergeno[chr]
            fgeno[fgeno .== "21"] .= "12"
            fgeno[fgeno .== "N1"] .= "1N"
            fgeno[fgeno .== "N2"] .= "2N"
        end
    end
    offformat = unique(reduce(vcat, [unique(i[!,:offspringformat]) for i=magicgeno.markermap]))
    if !issubset(offformat,["GT_haplo", "GT_unphased","GT_phased","GP","AD"])
        error(string("offspringformat=", offformat, ", not in [GT_haplo, GT_unphased, GT_phased, GP, AD]"))
    end
    if in("GT_unphased",offformat)
        for chr in 1:length(magicgeno.offspringgeno)
            offgeno = magicgeno.offspringgeno[chr]
            offgeno[offgeno .== "21"] .= "12"
            offgeno[offgeno .== "N1"] .= "1N"
            offgeno[offgeno .== "N2"] .= "2N"
        end
    end
    founderformat, offformat
end

function setmarkermap!(magicgeno::MagicGeno,markermap::DataFrame)
    merge_chromosome!(magicgeno)
    oldsnps = only(magicgeno.markermap)[!,:marker]
    oldsnpdict = Dict(oldsnps .=> 1:length(oldsnps))
    common_snps = intersect(oldsnps,markermap[!,:marker])
    b = [in(i,common_snps) for i in markermap[!,:marker]]
    groupmarkermap = groupby(markermap[b,:],:linkagegroup)
    snplsls = [[get(oldsnpdict,i,nothing) for i in chrmarkermap[!,:marker]]
        for chrmarkermap in groupmarkermap]
    b = length.(snplsls) .> 0
    if sum(b) < length(snplsls)
        snplsls = snplsls[b]
        groupmarkermap = groupmarkermap[b]
    end
    magicgeno.markermap = [begin
        chrmarkermap = groupmarkermap[i]
        oldmap = only(magicgeno.markermap)[snplsls[i],:]
        chrmarkermap[!,:marker] == oldmap[!,:marker] || @error "inconsistent markers"
        hcat(chrmarkermap[:,1:3],oldmap[:,4:end])
    end for i in 1:length(groupmarkermap)]
    magicgeno.foundergeno = [only(magicgeno.foundergeno)[i,:]  for i in snplsls]
    magicgeno.offspringgeno = [only(magicgeno.offspringgeno)[i,:]  for i in snplsls]
    magicgeno
end

function saveby_chromosome(magicgeno::MagicGeno;     
    nworker::Integer, 
    workdir::AbstractString=pwd(),    
	outstem::AbstractString="outstem")
    genofilels = genofileby_chromosome(magicgeno; nworker, outstem, workdir)
    for chr in eachindex(genofilels)
        subgeno = submagicgeno(magicgeno,chrsubset=[chr])
        savemagicgeno(genofilels[chr],subgeno)
    end
    genofilels
end

function saveby_chromosome(magicgeno::MagicGeno, outfiles::AbstractVector)    
    nchr = length(magicgeno.markermap)
    nchr == length(outfiles) || error(string("number of outfiles != nchr"))
    for chr in 1:nchr
        subgeno = submagicgeno(magicgeno,chrsubset=[chr])
        savemagicgeno(outfiles[chr],subgeno)
    end
    outfiles
end

function get_chroo(nsnpls::AbstractVector)
    chroo = sortperm(nsnpls;rev=true)
	reverse!(view(chroo,2:2:length(chroo)))
    chroo
end


function genofileby_chromosome(magicgeno::MagicGeno; 
    nworker::Integer, 
	outstem::AbstractString,
    workdir::AbstractString=pwd())    
    chridls = [i[1,:linkagegroup] for i in magicgeno.markermap]
	nsnpls = [size(i,1) for i in magicgeno.markermap]
    chroo = MagicBase.get_chroo(nsnpls)			
	noff = size(magicgeno.magicped.offspringinfo,1)	
	genofilels = [getabsfile(workdir,string(outstem*"_",chrid)) for chrid in chridls]
	nwork = min(nworker,length(genofilels))	
	tsleep = mean(nsnpls[chroo[1:nwork]]) * noff * (3e-6)
	for i in eachindex(chroo)
		if i<= nwork
			genofilels[chroo[i]] = string(genofilels[chroo[i]], "_tsleep", round(Int,(i-1)*tsleep),".csv.gz")
		else
			genofilels[chroo[i]] = string(genofilels[chroo[i]],"_tsleep0.csv.gz")
		end
	end
	genofilels
end	


function readby_chromosome!(magicgeno::MagicGeno, genofiles::AbstractVector;
    dropoffspringgeno::Bool = false,
    workdir::AbstractString=pwd())
    if isempty(genofiles) 
        @warn "genofiles are empty"
        return magicgeno
    end
    # assume that magicped and misc are the same among genofiles
    for i in eachindex(genofiles)        
        magicgeno2 = readmagicgeno(genofiles[i]; workdir)
        for (key, val) in magicgeno2.misc
            pushmisc!(magicgeno, key => val)
        end
        if i == 1
            magicgeno.magicped = magicgeno2.magicped 
            magicgeno.markermap = magicgeno2.markermap
            magicgeno.foundergeno = magicgeno2.foundergeno            
            magicgeno.offspringgeno = dropoffspringgeno ? nothing : magicgeno2.offspringgeno            
        else
            push!(magicgeno.markermap,magicgeno2.markermap[1])
            push!(magicgeno.foundergeno,magicgeno2.foundergeno[1])            
            dropoffspringgeno || push!(magicgeno.offspringgeno,magicgeno2.offspringgeno[1])
            offcols = names(magicgeno.magicped.offspringinfo)
            offcols2 = names(magicgeno2.magicped.offspringinfo)
            d = setdiff(offcols2,offcols)
            if !isempty(d)
                col1 = offcols[1]
                col1 == offcols2[1] || @error "inconsisent col1"
                magicgeno.magicped.offspringinfo[!,col1] == magicgeno2.magicped.offspringinfo[!,col1] || @error "inconsisent col1"
                magicgeno.magicped.offspringinfo = hcat(magicgeno.magicped.offspringinfo,magicgeno2.magicped.offspringinfo[:,d])
            end
        end
    end
    magicgeno
end

function readby_chromosome(genofiles::AbstractVector;
    dropoffspringgeno::Bool = false,
    workdir::AbstractString=pwd())
    if isempty(genofiles) 
        @warn "genofiles are empty"
        return nothing
    end
    # assume that magicped and misc are the same among genofiles
    magicgeno  = readmagicgeno(genofiles[1]; workdir)
    dropoffspringgeno && (magicgeno.offspringgeno = nothing)
    for i in 2:length(genofiles)        
        magicgeno2 = readmagicgeno(genofiles[i]; workdir)
        push!(magicgeno.markermap,magicgeno2.markermap[1])
        push!(magicgeno.foundergeno,magicgeno2.foundergeno[1])
        dropoffspringgeno || push!(magicgeno.offspringgeno,magicgeno2.offspringgeno[1])
    end
    magicgeno
end

function splitby_chromosome!(magicgeno::MagicGeno)
    b = length(magicgeno.markermap) == length(magicgeno.foundergeno) == length(magicgeno.offspringgeno)
    b || @error "#chromosomes in magicgeno  are not consistent"
    chrls = only(magicgeno.markermap)[!,:linkagegroup]
    snplsls = [chrls .== i for i in unique(chrls)]
    magicgeno.markermap = [only(magicgeno.markermap)[i,:] for i in snplsls]
    magicgeno.foundergeno = [only(magicgeno.foundergeno)[i,:]  for i in snplsls]
    magicgeno.offspringgeno = [only(magicgeno.offspringgeno)[i,:]  for i in snplsls]
    magicgeno
end

function merge_chromosome!(magicgeno::MagicGeno)
    magicgeno.markermap = [reduce(vcat,magicgeno.markermap)]
    magicgeno.foundergeno=[reduce(vcat,magicgeno.foundergeno)]
    magicgeno.offspringgeno=[reduce(vcat,magicgeno.offspringgeno)]
    magicgeno
end

function get_offspringformats(magicgeno::MagicGeno)
    if isnothing(magicgeno.markermap)
        nothing
    else
        unique(reduce(vcat,[unique(i[!,:offspringformat]) for i in magicgeno.markermap]))
    end
end
function info_magicgeno(magicgeno::MagicGeno;
    io::Union{Nothing,IO}=nothing,
    verbose::Bool=true)
    if !isnothing(magicgeno.markermap)
        founderformat = unique(reduce(vcat,[unique(i[!,:founderformat])
            for i in magicgeno.markermap]))
        offformat = unique(reduce(vcat,[unique(i[!,:offspringformat])
            for i in magicgeno.markermap]))
        printconsole(io,verbose,string("founderformat = ", join(founderformat,",")))
        printconsole(io,verbose,string("offspringformat = ", join(offformat,",")))
    end
    nchr=length(magicgeno.markermap)
    noff=size(magicgeno.magicped.offspringinfo,1)
    nfounder=size(magicgeno.magicped.founderinfo,1)
    nsnp=sum(size.(magicgeno.markermap,1))
    nsubpop = length(unique(magicgeno.magicped.offspringinfo[!,:member]))
    msg = string("#founder=",nfounder,", #offspring=",noff,
        ", #subpop=", nsubpop, ", #chr=",nchr, ", #snp=",nsnp)
    printconsole(io,verbose,msg)
end

function check_required_geneticmap(magicgeno::MagicGeno; io::Union{IO,Nothing}=nothing)
    nmarker = sum(size(i,1) for i in magicgeno.markermap)
    nmiss_poscm = sum([sum(.!(isa.(i[!,:poscm],Real))) for i in magicgeno.markermap])
    if nmiss_poscm == nmarker
        msg = string("genetic map is missing! suggest to set isphysmap and recomrate if physmap exist, or infer genetic map using magicmap!")
        printconsole(io, false, "ERROR: "*msg)
        error(msg)
    else
        if nmiss_poscm > 0 
            msg = string(nmiss_poscm, " out of ",nmarker, " markers have missing genetic position")
            printconsole(io, false, "ERROR: "*msg)
            error(msg)
        end
    end
end

function get_isorderedls(chrmapls; 
    isphysmap::Bool = false)
    poscol = isphysmap ? :physposbp : :poscm    
    isorderedls = [begin 
        ismiss = all(ismissing.(chrmap[!,poscol]))
        if ismiss 
            missing
        else
            issorted(skipmissing(chrmap[!,poscol]))
        end
    end for chrmap in chrmapls]    
    isorderedls
end

function check_markerorder(magicgeno::MagicGeno;
    io::Union{IO,Nothing}=nothing, verbose::Bool=true)
    res = true    
    is_miss_bp = any([any(ismissing.(i[!,:physposbp])) for i in magicgeno.markermap])
    is_miss_cm = any([any(ismissing.(i[!,:poscm])) for i in magicgeno.markermap])
    if is_miss_bp && is_miss_cm 
        @warn string("there exist missing posistions in both physposbp and poscm!")
        return false
    end
    for col in [:poscm, :physposbp]
        # colchr = col == :poscm ? :linkagegroup : :physchrom
        # chridls = [i[1,colchr] for i in magicgeno.markermap]
        boolmiss = col==:poscm ? is_miss_cm : is_miss_bp
        if !boolmiss
            bool = issorted.(magicgeno.markermap,col)
            poskind = col == :poscm ? "Genetic positions" : "Physical positions"
            if all(bool)
                msg = string(poskind, " of markers in each chromosome are in non-descreasing order")
                MagicBase.printconsole(io,verbose,msg)
            else                
                msg = string(poskind, " of markers are not in non-descreasing order!")
                @warn msg
                MagicBase.printconsole(io,false,"warning: "*msg)
                res = false
            end
        end
    end    
    res
end


function get_errorls(magicgeno::MagicGeno,chr::Integer;
    errorname::Symbol,
    io::Union{Nothing,IO}=nothing,
    verbose::Bool=true)
    if errorname == :peroffspringerror
        chrid = magicgeno.markermap[chr][1,:linkagegroup]
        errcol = string(errorname, "_",chrid)
        offinfo = magicgeno.magicped.offspringinfo
        if in(errcol,names(offinfo))
            res = Vector(offinfo[!, errcol])
            msg = string("chr=", chrid, ", ", errorname, "=", round(mean(skipmissing(res)),digits=4),"(",
                join(round.(extrema(skipmissing(res)),digits=4),"~"),")")
            printconsole(io, verbose, msg)
        else
            return nothing
        end
    else
        if eltype(magicgeno.markermap[chr][!,errorname]) === Missing
            res = nothing
        else
            res = Vector(magicgeno.markermap[chr][!,errorname])
            chrid = magicgeno.markermap[chr][1,:linkagegroup]
            msg = string("chr=", chrid, ", ", errorname, "=", round(mean(skipmissing(res)),digits=4),"(",
                join(round.(extrema(skipmissing(res)),digits=4),"~"),")")
            printconsole(io, verbose, msg)
        end
    end
    res
end


function sortmarker!(magicgeno::MagicGeno; byphysmap::Bool=false)
    if byphysmap 
        chrcol,poscol = :physchrom, :physposbp
    else
        chrcol,poscol = :linkagegroup, :poscm
    end
    b = [any(ismissing.(i[!,chrcol]) .|| ismissing.(i[!,poscol])) for i in magicgeno.markermap]
    if any(b)
        @warn string("physical map is missing for some markers; sortmarker! is not perform")
    else
        merge_chromosome!(magicgeno)
        chridls = sort(unique(magicgeno.markermap[1][!,chrcol]))
        ls = [begin 
            posls = findall(magicgeno.markermap[1][!,chrcol] .== chrid)
            posls[sortperm(magicgeno.markermap[1][posls,poscol])]
            # all(markermap[1][posls,:physchrom] .== chrid)
        end for chrid in chridls]
        magicgeno.markermap = [magicgeno.markermap[1][i,:] for i in ls]
        magicgeno.foundergeno =  [magicgeno.foundergeno[1][i,:] for i in ls]
        magicgeno.offspringgeno = [magicgeno.offspringgeno[1][i,:] for i in ls]        
    end
    magicgeno
end


function splitby_connectedped(magicgeno::MagicGeno)
    ccpops = get_connected_pops(magicgeno.magicped)
    isempty(ccpops) && error("unexpected connected pops=",ccpops)
    if length(ccpops) ==1 
        [magicgeno]
    else
        [submagicgeno(magicgeno,pops) for pops in ccpops]
    end
end

function submagicgeno(magicgeno::MagicGeno, subpops::AbstractVector)
    magicped = magicgeno.magicped
    memls = unique(magicped.offspringinfo[!,:member])
    issubset(subpops,memls) || error(string("subpops=",subpops, " is not a subset of subpopulations = ",memls))
    if isa(magicped.designinfo, Dict{String, DesignInfo})        
        newdesigninfo = Dict([convert(String,subpop)=>magicped.designinfo[subpop] for subpop in subpops])         
        founders = unique(reduce(vcat, [newdesigninfo[subpop].founders for subpop in subpops]))            
    elseif isa(magicped.designinfo, Pedigree)
        newdesigninfo = Pedigrees.getsubped(magicped.designinfo, subpops)    
        founders = Pedigrees.getfounderid(newdesigninfo)    
    else
        @error string("Could not extract ped from designinfo=",magicped.designinfo)
        return nothing        
    end
    # new founder
    iskeep_founder = [in(i,founders) for i in magicped.founderinfo[!,:individual]]
    newfounderinfo = magicped.founderinfo[iskeep_founder,:]
    newfoundergeno = [i[:, iskeep_founder] for i in magicgeno.foundergeno]
    # new offspring
    iskeep_off = [in(i,subpops) for i in magicped.offspringinfo[!,:member]]
    newoffinfo = magicped.offspringinfo[iskeep_off,:]
    newoffgeno = [i[:, iskeep_off] for i in magicgeno.offspringgeno]
    newmagicped = MagicPed(newdesigninfo,newfounderinfo,newoffinfo)
    MagicGeno(newmagicped,deepcopy(magicgeno.markermap),newfoundergeno,newoffgeno,deepcopy(magicgeno.misc))
end
