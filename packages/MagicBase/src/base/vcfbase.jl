
function vcf_count_markers(vcffile::AbstractString;  commentstring = "##",workdir::AbstractString=pwd())
    inext = last(MagicBase.split_allext(vcffile))
    inext in [".vcf",".vcf.gz"] || @error string("vcffile ext must be .vcf or .vcf.gz")
    in_open = inext in [".vcf.gz"] ? GZip.open : open
    vcffile2 = getabsfile(workdir,vcffile)
    lastcomment = MagicBase.findlastcomment(vcffile2; commentstring)
    nlines = in_open(vcffile2,"r") do io; countlines(io) end
    nmarkers =  nlines - lastcomment - 1 # 1 denotes title line
    nmarkers
end

function vcf_count_samples(vcffile::AbstractString;  commentstring = "##",workdir::AbstractString=pwd())
    length(vcf_get_samples(getabsfile(workdir,vcffile);commentstring))
end

function vcf_get_markers(vcffile::AbstractString;  commentstring = "##",workdir::AbstractString=pwd())
    inext = last(MagicBase.split_allext(vcffile))
    inext in [".vcf",".vcf.gz"] || @error string("vcffile ext must be .vcf or .vcf.gz")
    in_open = inext in [".vcf.gz"] ? GZip.open : open
    vcffile2 = getabsfile(workdir,vcffile)
    lastcomment = MagicBase.findlastcomment(vcffile2; commentstring)
    in_open(vcffile2,"r") do io
        for _ in 1:lastcomment
            readline(io)
        end
        # read title row
        readline(io)
        # read marker id    
        markerls = Vector{String}()
        while !eof(io)
            # vcfcols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL","FILTER", "INFO", "FORMAT"]
            line = split(readline(io;keep=false),"\t")
            push!(markerls,strip(line[3]))
        end    
        markerls
    end
end

function csv_get_samples(csvgenofile::AbstractString; commentstring = "##",workdir::AbstractString=pwd())
    inext = last(MagicBase.split_allext(csvgenofile))
    inext in [".csv",".csv.gz"] || @error string("csvgenofile ext must be .csv or .csv.gz")
    in_open = inext in [".csv.gz"] ? GZip.open : open
    csvgenofile2 = getabsfile(workdir,csvgenofile)
    lastcomment = MagicBase.findlastcomment(csvgenofile2; commentstring)
    in_open(csvgenofile2,"r") do io
        for _ in 1:lastcomment
            readline(io)
        end
        # parse title row    
        string.(split(readline(io,keep=false),",")[_col_1stsample:end])
    end
end


function vcf_get_samples(vcffile::AbstractString; commentstring = "##",workdir::AbstractString=pwd())
    inext = last(MagicBase.split_allext(vcffile))
    inext in [".vcf",".vcf.gz"] || @error string("vcffile ext must be .vcf or .vcf.gz")
    in_open = inext in [".vcf.gz"] ? GZip.open : open
    vcffile2 = getabsfile(workdir,vcffile)
    isfile(vcffile2) || error(string(vcffile2, " does not exist"))
    lastcomment = MagicBase.findlastcomment(vcffile2; commentstring)    
    in_open(vcffile2,"r") do io
        for _ in 1:lastcomment
            readline(io)
        end
        # parse title row
        string.(split(readline(io,keep=false),"\t")[10:end])
    end
end

"""
    vcf_pad_samples(vcffile; keyargs...)

pad samples into vcffiles with all genotypes being missing. 

# Positional arguments

`vcffile::AbstractString`: vcf genofile. 

# Keyword arguments

`padsamples::AbstractVector`: a list of samples to be padded 

`commentstring::AbstractString="##"`: the lines beginning with  are ignored

`outstem::AbstractString=popid`: stem of output filename.

`workdir::AbstractString=pwd()`: directory for reading and saving files.

"""
function vcf_pad_samples(vcffile::AbstractString; 
    padsamples::AbstractVector,
    commentstring::AbstractString="##",
    outstem::AbstractString="oustem",
    workdir::AbstractString=pwd())
    inext = last(MagicBase.split_allext(vcffile))
    inext in [".vcf",".vcf.gz"] || @error string("vcffile ext must be .vcf or .vcf.gz")
    in_open = inext in [".vcf.gz"] ? GZip.open : open    
    lastcomment = MagicBase.findlastcomment(vcffile; commentstring)    
    in_open(getabsfile(workdir,vcffile),"r") do io
        GZip.open(getabsfile(workdir,outstem*".vcf.gz"), "w") do outio
            for _ in 1:lastcomment
                line = readline(io;keep=true)
                write(outio, line)
            end
            # parse title row    
            titlerow = split(readline(io,keep=false),"\t")
            commonsamples = intersect(titlerow[10:end],padsamples)
            if !isempty(commonsamples) 
                setdiff!(padsamples, commonsamples)
                msg = string("samples already in vcffile: ",commonsamples, ", reset padsamples=",padsamples)
                @warn msg
            end
            append!(titlerow,padsamples)
            write(outio, join(titlerow,"\t"),"\n")
            # parse data
            while !eof(io)
                line = readline(io;keep=false)
                write(outio, line*"\t.\n")
            end
        end
    end
end

function vcf_del_samples(vcffile::AbstractString;
    delsamples::AbstractVector,
    commentstring = "##",
    outstem::AbstractString="outstem",
    logfile::Union{Nothing,AbstractString} = nothing,
    workdir::AbstractString=pwd(),
    verbose::Bool=true)
    if isempty(delsamples)
        @warn "empty delsamples"
        return nothing
    end
    starttime = time()
    logio = MagicBase.set_logfile_begin(logfile, workdir, "vcf_del_samples"; verbose,delim="-")    
    msg = string("list of args: \n",
        "vcffile = ", vcffile, "\n",        
        "delsamples = ", delsamples, "\n",        
        "commentstring = ", commentstring, "\n",
        "outstem = ", outstem, "\n",        
        "logfile = ", logfile, "\n",
        "workdir = ", workdir, "\n",
        "verbose = ", verbose)
    printconsole(logio,verbose,msg)
    inext = last(MagicBase.split_allext(vcffile))
    inext in [".vcf",".vcf.gz"] || @error string("vcffile ext must be .vcf or .vcf.gz")
    in_open = inext in [".vcf.gz"] ? GZip.open : open    
    lastcomment = MagicBase.findlastcomment(vcffile; commentstring)    
    outfile = getabsfile(workdir,outstem*"_geno.vcf.gz")
    in_open(getabsfile(workdir,vcffile),"r") do io
        GZip.open(outfile, "w") do outio
            for _ in 1:lastcomment
                line = readline(io;keep=true)
                write(outio, line)
            end
            # parse title row    
            titlerow = split(readline(io,keep=false),"\t")
            samples = string.(strip.(titlerow[10:end]))
            if isempty(delsamples)
                delsamples2 = []
            elseif eltype(delsamples) <: Integer 
                delsamples2 = samples[delsamples]
            elseif eltype(delsamples) <: AbstractString
                delsamples2 = string.(strip.(delsamples))
                d = setdiff(delsamples2, samples)
                if !isempty(d)
                    @warn string("delsamples that are not in vcffile: ",d)
                end
                intersect!(delsamples2,samples)
            else
                @error string("unknown type=",eltype(delsamples), " for input delsamples=",delsamples)
            end               
            dict = Dict(samples .=> 1:length(samples))
            cols = [get(dict,old,nothing) for old in delsamples2]
            any(isnothing.(cols)) && error(string("unxecpected cols=",cols))
            cols .+= 9
            deleteat!(titlerow, cols)
            write(outio, join(titlerow,"\t"),"\n")
            # parse data
            while !eof(io)                
                line = split(readline(io,keep=false),"\t")
                deleteat!(line, cols)
                write(outio, join(line,"\t"),"\n")
            end
        end
    end
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, logio, tused,"vcf_del_samples"; verbose,delim="-")    
    outfile    
end

function vcf_rename_samples(vcffile::AbstractString; 
    renamerule::Union{AbstractDict, Vector{T} where T <:Pair},
    commentstring = "##",
    outstem::AbstractString="outstem",
    logfile::Union{Nothing,AbstractString} = nothing,
    workdir::AbstractString=pwd(),
    verbose::Bool=true)
    if isempty(renamerule)
        @warn "empty renamerule"
        return nothing
    end
    starttime = time()
    logio = MagicBase.set_logfile_begin(logfile, workdir, "vcf_rename_samples"; verbose,delim="-")    
    msg = string("list of args: \n",
        "vcffile = ", vcffile, "\n",        
        "renamerule = ", renamerule, "\n",        
        "commentstring = ", commentstring, "\n",
        "outstem = ", outstem, "\n",        
        "logfile = ", logfile, "\n",
        "workdir = ", workdir, "\n",
        "verbose = ", verbose)
    printconsole(logio,verbose,msg)
    inext = last(MagicBase.split_allext(vcffile))
    inext in [".vcf",".vcf.gz"] || @error string("vcffile ext must be .vcf or .vcf.gz")
    in_open = inext in [".vcf.gz"] ? GZip.open : open    
    vcffile2 = getabsfile(workdir,vcffile)
    lastcomment = MagicBase.findlastcomment(vcffile2; commentstring)    
    outfile = getabsfile(workdir,outstem*".vcf.gz")
    in_open(vcffile2,"r") do io
        GZip.open(outfile, "w") do outio
            for _ in 1:lastcomment
                line = readline(io;keep=true)
                write(outio, line)
            end
            # parse title row    
            titlerow = split(readline(io,keep=false),"\t")
            samples = string.(strip.(titlerow[10:end]))
            dict = Dict(samples .=> 1:length(samples))
            col_new = [[get(dict,old,nothing),old, new] for (old, new) in renamerule]
            cols = first.(col_new)
            news = [i[3] for i in col_new]
            b = isnothing.(cols)
            if any(b)
                msg = string("samples in renamedict but not in vcffile: ",[i[2] for i in col_new[b]])                
                @warn msg
                printconsole(logio, false, "Warning: "*msg)
                cols = cols[.!b]
                news = news[.!b]
            end
            if isempty(cols)            
                msg = string("no samples were renamed!")                
                @warn msg
                printconsole(logio, false, "Warning: "*msg)
            else
                msg = string("rename samples: ", samples[cols] .=> news)
                printconsole(logio, verbose, msg)
                cols .+= 9
                titlerow[cols] .= news            
            end
            write(outio, join(titlerow,"\t"),"\n")
            # parse data
            while !eof(io)
                line = readline(io;keep=true)
                write(outio, line)
            end
        end
    end
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, logio, tused,"vcf_rename_samples"; verbose,delim="-")    
    outfile    
end

function vcf_get_subgeno(vcffile::AbstractString; 
    subsamples::AbstractVector=[],
    commentstring = "##",    
    workdir::AbstractString=pwd())
    inext = last(MagicBase.split_allext(vcffile))
    inext in [".vcf",".vcf.gz"] || @error string("vcffile ext must be .vcf or .vcf.gz")
    in_open = inext in [".vcf.gz"] ? GZip.open : open    
    vcffile2 = getabsfile(workdir,vcffile)
    lastcomment = MagicBase.findlastcomment(vcffile2; commentstring)    
    in_open(vcffile2,"r") do io        
        for _ in 1:lastcomment
            readline(io;keep=true)            
        end
        # parse title row    
        titlerow = split(readline(io,keep=false),"\t")
        samples = strip.(titlerow[10:end])
        if isempty(subsamples)
            subsamples2 = []
        elseif eltype(subsamples) <: Integer 
            subsamples2 = samples[subsamples]
        elseif eltype(subsamples) <: AbstractString
            subsamples2 = string.(strip.(subsamples))
            d = setdiff(subsamples2, samples)
            if !isempty(d)
                @warn string("subsamples that are not in vcffile: ",d)
            end
            intersect!(subsamples2,samples)
        else
            @error string("unknown type=",eltype(subsamples), " for input subsamples=",subsamples)
        end           
        vcfcols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL","FILTER", "INFO", "FORMAT"]
        leftcols= 1:9
        leftpairs = [Symbol(i)=>String[] for i in vcfcols[leftcols]]
        if isempty(subsamples2)
            resdf = DataFrame(leftpairs...)
            cols = leftcols
        else
            resdf = DataFrame(leftpairs..., [Symbol(i) => String[] for i in subsamples2]...)            
            cols = [findfirst(==(i),samples) + 9 for i in subsamples2]
            pushfirst!(cols,leftcols...) 
        end        
        while !eof(io)            
            line = split(readline(io,keep=false),"\t")
            push!(resdf, line[cols])            
        end
        resdf
    end
end

"""
    vcffilter(vcffile; keyargs...)

filter markers line by line. 

# Positional arguments

`vcffile::AbstractString`: vcf genofile. 

# Keyword arguments

`setmarkerid::Union{Nothing,Bool}=nothing`: if true, set markerid. If it is nothing, setmarkerid = true only if markerid is missing

`delsamples::Union{Nothing,AbstractVector}=nothing`: list of sample IDs to be deleted. If it is nothing, no deletion of samples 

`deldupe::Bool=false`: if true, delete sucessive markers that have exactly duplicated genotypes in format of GT

`isdelmultiallelic::Bool=true`: if true, delete markers with >2 alleles. 

`isdelmonomorphic::Bool=true`: if true, delete markers with single allele. 

`seqstretch::Integer=0`: delete non-initial markers in a sequence stretch of length <= seqstretch (in bp), assuming marker are ordered by physical positions. If it is not positive, no filtering for short streches.

`maxmiss::Real = 0.99`: delete markers with missing fraction > maxmiss

`minmaf::Real = 0.01`: delete markers with minor allele frequency < minmaf

`commentstring::AbstractString="##"`: the lines beginning with  are ignored

`outstem::AbstractString=popid`: stem of output filename.

`workdir::AbstractString=pwd()`: directory for reading and saving files.

`logfile::AbstractString = outstem*"_vcffilter.log"`: log filename. 

`verbose::Bool=true`: if true, print details on the stdout.

"""
function vcffilter(vcffile::AbstractString;
    setmarkerid::Union{Nothing,Bool}=nothing,  # if it is nothign, set it in format of CHROM_POS only if markerid is missing
    delsamples::Union{Nothing,AbstractVector}=nothing,     
    deldupe::Bool=false, 
    isdelmultiallelic::Bool=true,    
    isdelmonomorphic::Bool=true,    
    keeponlypass::Bool=true, 
    seqstretch::Integer=0, # 64 bp for length of sequence tags
    maxmiss::Real = 0.99,
    minmaf::Real = 0.01,
    commentstring::AbstractString="##",
    outstem::AbstractString = first(split_allext(basename(vcffile))),      
    workdir::AbstractString=pwd(),
    logfile::AbstractString = outstem*"_vcffilter.log",
    verbose::Bool=true)    
    starttime = time()
    logio = MagicBase.set_logfile_begin(logfile, workdir, "vcffilter"; verbose,delim="-")
    inext = last(split_allext(vcffile))
    inext in [".vcf",".vcf.gz"] || @error string("vcffile ext must be .vcf or .vcf.gz")
    msg = string("list of args: \n",
        "vcffile = ", vcffile, "\n",        
        "setmarkerid = ", setmarkerid, "\n",        
        "delsamples = ", delsamples, "\n",        
        "deldupe = ", deldupe, " \n",        
        "isdelmonomorphic = ", isdelmonomorphic, "\n",        
        "isdelmultiallelic = ", isdelmultiallelic, "\n",                
        "keeponlypass = ", keeponlypass, "\n",                
        "seqstretch = ", seqstretch, " \n",        
        "maxmiss = ", maxmiss, "\n",        
        "minmaf = ", minmaf, "\n",        
        "commentstring = ", commentstring, "\n",
        "outstem = ", outstem, "\n",        
        "logfile = ", logfile, "\n",
        "workdir = ", workdir, "\n",
        "verbose = ", verbose)
    printconsole(logio,verbose,msg)    
    msg = deldupe ? "delete sucessive duplicate markers\n" : ""
    seqstretch > 0 && (msg *= string("delete noninitial markers in a sequence stretch of length < ", seqstretch, "\n"))
    msg *= string("filter sequentially for biallelic, non-monomorphic, missing <= ", maxmiss, ", and maf >= ",minmaf)
    printconsole(logio,verbose,msg)        
    in_open = inext in [".vcf.gz"] ? GZip.open : open
    vcffile2=getabsfile(workdir,vcffile)
    lastcomment = findlastcomment(vcffile2; commentstring)    
    outext = inext
    in_open(vcffile2,"r") do inio        
        out_open = outext in [".vcf.gz"] ? GZip.open : open
        outfile = outstem*"_vcffilter"*outext
        outfile2 = getabsfile(workdir,outfile)
        out_open(outfile2, "w") do outio
            # commentstring = "##" for vcf file
            filedate = string(commentstring*"fileDate=",replace(string(Date(now())),"-"=>""))
            write(outio, filedate,"\n")
            write(outio, commentstring*"source=RABBIT\n")
            for _ in 1:lastcomment
                line = readline(inio,keep=true)                
                write(outio, line)
            end            
            # parse title row    
            titlerow = split(readline(inio,keep=false),"\t")
            samples = string.(strip.(titlerow[10:end]))
            if isnothing(delsamples)
                delsamples2 = []
            elseif eltype(delsamples) <: Integer 
                delsamples2 = samples[delsamples]
            elseif eltype(delsamples) <: AbstractString
                delsamples2 = string.(strip.(delsamples))
                d = setdiff(delsamples2, samples)
                if !isempty(d)
                    @warn string("delsamples that are not in vcffile: ",d)
                end
                intersect!(delsamples2,samples)
            else
                @error string("unknown type=",eltype(delsamples), " for input delsamples=",delsamples)
            end               
            dict = Dict(samples .=> 1:length(samples))
            cols = [get(dict,old,nothing) for old in delsamples2]
            any(isnothing.(cols)) && error(string("unexpected cols=",cols))
            cols .+= 9
            deleteat!(titlerow, cols)
            write(outio, join(titlerow,"\t"),"\n")
            nind = length(titlerow)-9
            msg = string("#samples=", nind, ",sample (up to 10): ", join(titlerow[10:min(19,length(titlerow))],","))
            printconsole(logio,verbose,msg)                      

            # read marker line
            startt = time()            
            nmarker_incl = 0
            nmarker = 0
            ndupe = 0
            nnotpass = 0
            nstretch = 0            
            nmultia = 0            
            nmono = 0
            nmiss = 0
            nmaf = 0
            binchr = ""
            binpos = -10*seqstretch
            bingeno = ["" for _ in 1:nind]
            while !eof(inio)
                # rem(nmarker, 1000) == 0 && (startt = time())
                nmarker += 1
                rowgeno = split(readline(inio,keep=false),"\t")
                deleteat!(rowgeno, cols)           
                # vcf column ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL","FILTER", "INFO", "FORMAT"]
                formatls = split(rowgeno[9],":")  # col9 = FORAMT
                gt = findfirst(x->x=="GT",formatls)
                if isnothing(gt) 
                    @error string("GT format does not exist! marker_row=",rowgeno) maxlog = 10
                end
                genols = [i[gt] for i in split.(rowgeno[10:end],":")]                
                nowchr = rowgeno[1]
                nowpos = tryparse(Int, rowgeno[2])                                
                nowgeno = genols
                ismultia =length(split(rowgeno[5],",")) > 1  # col5=alternative      
                if keeponlypass && (rowgeno[7] != "." && !occursin("PASS", uppercase(rowgeno[7]))) #col7 = FILTER
                    isdel = true
                    nnotpass += 1                    
                elseif nowchr == binchr && all(nowgeno .== bingeno) && deldupe
                    isdel = true 
                    ndupe += 1
                else                    
                    if nowchr == binchr && !isnothing(binpos) && !isnothing(nowpos) && abs(nowpos - binpos)+1 <= seqstretch
                        isdel = true 
                        nstretch += 1
                    else
                        binchr = nowchr
                        binpos = nowpos
                        bingeno = nowgeno
                        if ismultia 
                            nmultia += 1
                            isdel = isdelmultiallelic                    
                        else
                            isdel = false
                        end            
                    end                
                end                
                if !isdel
                    # filter for non-monomorphic        
                    # vcf column ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL","FILTER", "INFO", "FORMAT"]                                               
                    alleles = replace(string(unique(genols)...),"/"=>"","|"=>"","."=>"")
                    alleleset = unique(split(alleles,""))
                    ismono = isempty(alleles) || length(alleleset) == 1    # isempty(alleles) == 100% missing
                    if ismono
                        nmono += 1
                        isdel = isdelmonomorphic
                    end                
                    if length(alleleset) >=3  && !ismultia && !deldupe
                        msg = string("inconsistent multi-allelic: marker index = ", nmarker, ", allele_set=",alleleset, ", ALT=",alt)
                        @warn msg maxlog = 20
                        printconsole(logio, false, msg)                            
                    end                        
                    if !isdel                        
                        # filter missing                        
                        if ismono
                            if minmaf > 0 && maxmiss < 1.0
                                nmaf += 1
                                isdel = true
                            end
                        else
                            freqmiss = mean([in(i,[".","./.",".|."]) for i in genols])
                            if freqmiss > maxmiss 
                                nmiss += 1
                                isdel = true
                            else
                                # filter for maf
                                if minmaf > 0 
                                    alleles = replace(string(genols...),"/"=>"","|"=>"","."=>"")
                                    freqls = [count(i,alleles)/length(alleles) for i in alleleset]
                                    maf = 1.0- max(freqls...)
                                    if maf < minmaf || maf > 1-minmaf
                                        nmaf += 1
                                        isdel = true
                                    end
                                end
                            end
                        end
                    end
                end
                if !isdel 
                    # vcf column ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL","FILTER", "INFO", "FORMAT"]                                               
                    isimputeid = setmarkerid
                    if isnothing(setmarkerid) && strip(rowgeno[3]) == "."
                        isimputeid = true
                    end
                    if  !isnothing(isimputeid) && isimputeid
                        physchr, physpos = strip.(rowgeno[1:2])
                        if physchr == "." || physpos == "."
                            rowgeno[3] = string("snp",nmarker)
                        else
                            rowgeno[3] = string(physchr,"_", physpos)
                        end
                    end
                    write(outio, join(rowgeno,"\t"),"\n")                    
                    nmarker_incl += 1
                end
                if rem(nmarker, 1000) == 0
                    msg = string("#marker=", nmarker, 
                        ", #marker(keep)=", nmarker_incl, 
                        keeponlypass ? string(", #del_notpass=", nnotpass) : "",                         
                        deldupe ? string(", #del_dupe=", ndupe) : "",                         
                        seqstretch > 0 ? string(", #del_stretch=", nstretch) : "",                         
                        string(", #multi", isdelmultiallelic ? "=" : "(keep)=",nmultia), 
                        string(", #mono", isdelmonomorphic ? "=" : "(keep)=",nmono), 
                        ", #>maxmiss=", nmiss, ", #<minmaf=", nmaf,                         
                        ", tused=", round(time()-startt,digits=1),"s")
                    printconsole(logio,verbose,msg)
                    flush(outio)
                end
            end            
            msg = string("#individuals=", nind, ", #markers=",nmarker, 
                ", #marker(keep)=", nmarker_incl, 
                keeponlypass ? string(", #del_notpass=", nnotpass) : "",                         
                deldupe ? string(", #del_dupe=", ndupe) : "",                         
                seqstretch > 0 ? string(", #del_stretch=", nstretch) : "",                         
                string(", #multi", isdelmultiallelic ? "=" : "(keep)=",nmultia), 
                string(", #mono", isdelmonomorphic ? "=" : "(keep)=",nmono), 
                ", #>maxmiss=", nmiss, ", #<minmaf=", nmaf,                                         
            )
            printconsole(logio,verbose,msg)
        end
        msg = string("output vcffile: ", outfile)
        printconsole(logio,verbose, msg)
    end
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, logio, tused,"vcffilter"; verbose,delim="-")
    return 
end


"""
    merge_vcffiles(vcffiles; outstem, workdir)

merge vcffiles into a single vcf genofile. 

# Positional arguments

`vcffiles::AbstractVector`: a list of vcf genofile. 

# Keyword arguments

`outstem::AbstractString="outstem"`: stem of output filename.

`workdir::AbstractString=pwd()`: directory for reading and saving files.

"""
function merge_vcffiles(vcffiles::AbstractVector;
    outstem::AbstractString = "outstem",    
    workdir::AbstractString = pwd())
    vcffilels = [getabsfile(workdir,i) for i in vcffiles]
    genodf,comments = MagicBase.readgenofile(vcffilels[1])
    if all(genodf[!,:ID] .== ".")
        genodf[!,:ID] .= string.(genodf[!,"#CHROM"], "_",genodf[!,"POS"])
    end
    allunique(genodf[!,:ID]) || @warn "marker IDs are not unique"
    commentls = split(comments,"\n")
    filter!(x->!occursin("filedate",x),commentls)
    popfirst!(vcffilels)
    for vcffile in vcffilels        
        @info string("merging vcffile=",vcffile)
        genodf2,comments2 = MagicBase.readgenofile(vcffile)
        if all(genodf2[!,:ID] .== ".")
            genodf2[!,:ID] .= string.(genodf2[!,"#CHROM"], "_",genodf2[!,"POS"])
        end
        allunique(genodf2[!,:ID]) || @warn "marker IDs are not unique"
        # update genodf2 to have consistent REF and ALT with genodf
        snpdict = Dict(genodf[!,:ID] .=> 1:size(genodf,1))
        nswitch = 0
        nmismatch = 0
        for row in eachrow(genodf2)
            snp = row[:ID]
            if haskey(snpdict,snp)
                snpindex = snpdict[snp]
                # check ref and alt
                ref,alt = genodf[snpindex,[:REF,:ALT]]
                ref2,alt2 = row[[:REF,:ALT]]
                if ref == ref2 == "N" || alt == alt2 == "."
                    @warn string("consistency of alleles for marker=",snp,
                        " is not checked between current mergedgeno and vcffile=",vcffile) maxlog=5
                else                    
                    if (ref,alt) == (ref2,alt2)
                        0
                    elseif (ref,alt) == (alt2,ref2)
                        nswitch += 1
                        row[[:REF,:ALT]] .= [alt2,ref2]
                        for j in 10:length(row)
                            row[j] = replace(row[j],"1"=>"0","0"=>"1")
                        end                        
                    else
                        nmismatch += 1                        
                        for j in 10:length(row)
                            row[j] =  "."
                        end      
                    end
                end
                # check format
                format = genodf[snpindex, :FORMAT]
                format2 = row[:FORMAT]
                if format != format2
                    formatls = split(format, ":")
                    formatls2 = split(format2, ":")
                    maxformatmatch = 0
                    for i in 1:min(length(formatls), length(formatls))
                        if formatls[i] == formatls2[i] 
                            maxformatmatch += 1
                        else
                            break
                        end
                    end
                    if maxformatmatch == 0 
                        @error string("no format match; format=",format, "; format2=",format2) maxlog=5
                    else
                        newformat = join(formatls[1:maxformatmatch],":")
                        row[:FORMAT] = newformat
                        genodf[snpindex, :FORMAT] = newformat
                    end
                end
            end
        end
        if nswitch > 0
            @info string("switch REF and ALT for ", nswitch, " out of ", size(genodf2,1), " markers from vcffile=",vcffile)            
        end
        if nmismatch > 0
            msg = string("set genotypes missing for ", nmismatch, " markers from vcffile = ",vcffile)
            msg = string("because of inconsistent (ref,alt) between upto mergedgeno and the vcffile")
            @warn msg
        end    
        # TODO
        # outerjoin vcf
        genodf = outerjoin(genodf,genodf2; on = :ID,makeunique=true)
        commentls2 = split(comments2,"\n")
        filter!(x->!occursin("filedate",x),commentls2)
        union!(commentls,commentls2)
        # vcf column ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL","FILTER", "INFO", "FORMAT"]
        # combine columns 1:9; col=:ID
        cols = names(genodf, vcat(1:2,4:9))
        for c in cols
            cls = filter(x->occursin(c,x), names(genodf))
            setdiff!(cls,[c])
            geno_cls =  Matrix(genodf[!,vcat(cls,[c])])
            geno_c = [begin 
                temp = collect(skipmissing(unique(i)))
                deleteat!(temp, temp .== ".")
                temp
            end for i in eachrow(geno_cls)]            
            genodf[!,c] .= [length(i)==1 ? only(i) : "." for i in geno_c]
            if c in ["#CHROM","POS","FORMAT"]
                b = length.(geno_c) .> 1
                if any(b)
                    msg = string(sum(b), " markers have conflict values at col=",c)
                    @info msg
                end
            end
            select!(genodf,Not(cls))
        end
        # combine  columns for duplicated samples
        samplels = intersect(names(genodf)[10:end],names(genodf2)[10:end])
        for c in samplels            
            cls = filter(x->occursin(c,x), samplels)            
            intersect!(cls,[string(c,"_",i) for i in 1:length(cls)]) # resulting from makeunique=true            
            setdiff!(cls,[c])
            if !isempty(cls)                
                geno_cls =  Matrix(genodf[!,vcat(cls,[c])])
                geno_c = [collect(skipmissing(unique(i))) for i in eachrow(geno_cls)]
                genodf[!,c] .= [length(i)==1 ? i[1] : missing for i in geno_c]
                select!(genodf,Not(cls))
                ninconsist = sum(length.(geno_c) .> 1)
                msg = string("merging samples: ",vcat(cls,[c]), ", #inconsistent=", ninconsist, " (set to missing)")
                printconsole(logio,verbose,msg)        
            end
        end        
    end
    pushfirst!(commentls,string("##filedate=",string(now())))
    outfile = outstem*"_geno.vcf.gz"
    outfile2 = getabsfile(workdir,outfile)
    GZip.open(outfile2,"w") do file
        write(file,join(commentls,"\n"),"\n")
        CSV.write(file,genodf; 
            writeheader=true,
            append=true, 
            missingstring=".",
            delim = "\t"
        )
    end
    outfile2
end


"""
    resetmap(vcffile,mapfile;
        missingstring, commentstring, outstem, workdir)

exports a new vcf file with marker map replaced with mapfile. 

# Positional arguments

`vcffile::AbstractString`: genotypic data file with extension ".vcf" or ".vcf.gz".

`mapfile::AbstractString`: file for marker map, it can either be in VCF format or in CSV format. 
 For CSV-format, it must contain at least five columns: marker, linkagegroup, poscm, physchrom, physposbp. 
 The values are represented by missingstring. 

# Keyword arguments

`missingstring::AbstractString="NA"`: string representing missing value. 

`commentstring::AbstractString="##"`: the lines beginning with commentstring are ignored
in vcffile or mapfile.

`outstem::Union{Nothing,AbstractString}="outstem"`: stem of output filenames.

`workdir::AbstractString=pwd()`: directory for reading genfile and pedfile.    

"""
function resetmap(vcffile::AbstractString, mapfile::AbstractString;
    missingstring::AbstractString="NA",    
    commentstring = "##",
    outstem::AbstractString = "outstem",
    workdir::AbstractString=pwd())
    # it works also for genofile(vcffile) in CSV fromat
    # mapfile must have columns: marker, chromosome poscm, physposbp
    # genodf
    ext = last(MagicBase.split_allext(vcffile))
    isvcf = ext in [".vcf",".vcf.gz"]
    genodf,commentlines = readgenofile(vcffile;
        commentstring,missingstring, workdir)
    snpcol = isvcf ? :ID : :marker
    genosnps = string.(strip.(string.(genodf[!,snpcol])))    
    # mapdf: markers with missing chromosomes have been dropped
    mapdf = readmarkermap(mapfile; del_ungrouped=true, commentstring,missingstring=unique(vcat(missingstring,["missing",".","NA"])), workdir)
    mapsnps = string.(strip.(mapdf[!,:marker]))
    # snpdiff = setdiff(mapsnps, genosnps)
    # isempty(snpdiff) || @warn string("all markers in mapfile are not in vcffile")
    dict= Dict(genosnps .=> 1:length(genosnps))
    ii = [get(dict, i, nothing) for i in mapsnps]
    b = .!isnothing.(ii)
    sum(b) == 0 && @error string("all markers in mapfile are not in vcffile")
    resdf = genodf[ii[b],:]
    mapdf = mapdf[b,:]
    inputcolls = propertynames(mapdf)
    if isvcf
        #  for vcf, resdf columns: ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL","FILTER", "INFO", "FORMAT"]
        resdf[!,:INFO] .= map((x,y,z)->get_info_poscm(x,y,z),resdf[!,:INFO], mapdf[!,:linkagegroup], mapdf[!,:poscm])        
        inputcolls = propertynames(mapdf)
        if in(:physchrom,inputcolls)
            resdf[!,1] .= mapdf[!,:physchrom]
        end
        if in(:physposbp,inputcolls)
            resdf[!,2] .= mapdf[!,:physposbp]
        end
        resdf[!,3] == mapdf[!,:marker] || @error "inconsistent markers"        
    else
        # for csv, columns: [marker, linkagegroup,poscm,physchrom,physposbp]
        if issubset([:linkagegroup,:poscm,:physchrom, :physposbp],inputcolls)
            resdf[!,2:5] .= mapdf[!,[:linkagegroup,:poscm,:physchrom, :physposbp]]
        elseif  issubset([:linkagegroup,:poscm],inputcolls)
            resdf[!,2:3] .= mapdf[!,[:linkagegroup,:poscm]]
        else
            @error string("unknown cols=",inputcolls, " in mapfile=",mapfile)
        end
        resdf[!,:marker] == mapdf[!,:marker] || @error "inconsistent markers"
    end
    # save resdf
    outfile =string(outstem,"_resetmap", ext)
    outfile2 = MagicBase.getabsfile(workdir,outfile)
    myopen = ext in [".vcf.gz",".csv.gz"] ? GZip.open : open
    myopen(outfile2,"w") do io        
        cc = commentstring
        if isvcf            
            write(io, cc*"fileformat=VCFv4.3\n")
        else
            write(io, cc*"fileformat=CSV\n")
        end
        filedate = string(cc*"##filedate=",string(now()))
        write(io, filedate,"\n")
        write(io, cc*"source=RABBIT\n")
        isempty(commentlines) || write(io, commentlines,"\n")        
        CSV.write(io, resdf; delim= isvcf ? '\t' : ',',
            header=true, missingstring = isvcf ? "." : missingstring,append=true)
    end
    outfile2
end

function get_info_poscm(vcfinfo::AbstractString,linkagegroup, poscm::Union{Missing,Real})
    if vcfinfo != "."
        infols = split.(split(vcfinfo,";"),"=")
        oldinfo = ""
        b = length.(infols) .!= 2
        if any(b)
            oldinfo = join([length(i) > 2 ? join(i,"=") : i for i in infols[b]],";")
            @warn string("info without pairing: ", oldinfo) maxlog=10
            deleteat!(infols, b)
        end
        # unique(length.(infols)) == [2] || @error string("info with unexpected format, info = ",vcfinfo)
        infodict = OrderedDict([Pair(i...) for i in infols])
        if !ismissing(linkagegroup) && !ismissing(poscm)
            infodict["LINKAGEGROUP"] = string(linkagegroup) # insert or modify
            infodict["POSCM"] = string(poscm)        
        end
        newinfo = join(map(i->string(i[1],"=",i[2]),collect(infodict)),";")
        isempty(oldinfo) ? newinfo : string(newinfo, ";", oldinfo)
    else
        newinfo = string("LINKAGEGROUP=", linkagegroup, ";POSCM=",poscm)
    end
end


function readmarkermap(mapfile::AbstractString;
    missingstring=["NA","missing"],
    del_ungrouped::Bool = true,
    commentstring::AbstractString="##",
    workdir::AbstractString = pwd())
    ext = last(split_allext(mapfile))
    ext in [".csv",".csv.gz",".vcf",".vcf.gz"] || @error string("genofile ext must be .csv, .csv.gz, .vcf or .vcf.gz")
    if ext in [".vcf",".vcf.gz"]
        # readmarkermap_vcf is much faster than readgenodf 
        # since readmarkermap_vcf read only non-geno on the left part of dataframe
        mapdf = readmarkermap_vcf(mapfile; commentstring,missingstring, workdir)
    else
        mapdf, _ = readgenodf(mapfile;isdelmultiallelic=false,commentstring,missingstring, workdir)
    end
    if size(mapdf,2) < 3
        error(string(mapfile, " must contain #column >= 3"))
    end
    # mapdf cols:"marker", "linkagegroup","poscm","physchrom", "physposbp"
    mapdf = mapdf[!,1:min(_col_1stsample-1,size(mapdf,2))]
    if del_ungrouped        
        isallmiss = false
        b = .!ismissing.(mapdf[!,2]) 
        if any(b)
            mapdf = mapdf[b,:]
        else    
            if size(mapdf,2) <5
                isallmiss = true
            else
                b = .!ismissing.(mapdf[!,4]) 
                if any(b)
                    mapdf = mapdf[b,:]
                else
                    isallmiss = true
                end
            end
        end
        isallmiss && @warn string("marker grouping is missing")
    end
    parsemarkermap!(mapdf)
    mapdf
end


function readmarkermap_vcf(genofile::AbstractString;
    missingstring=["NA","missing"],
    commentstring::AbstractString="##",
    workdir::AbstractString = pwd())
    # vcf column ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL","FILTER", "INFO", "FORMAT"]
    markermap = MagicBase.vcf_get_subgeno(genofile; subsamples = [], commentstring, workdir);
    colnames= ["marker", "linkagegroup","poscm","physchrom", "physposbp", "info", "founderformat","offspringformat",
        "foundererror","offspringerror","baseerror","allelicbias","allelicoverdispersion","allelicdropout"]
    df = DataFrame([String[] for _ in 1:length(colnames)],colnames)
    for rowgeno in eachrow(markermap)
        physchrom,physposbp,snpid = [i== "." ? "NA" : i for i in rowgeno[1:3]]       
        info = rowgeno[8]
        # res: newinfo, linkagegroup, poscm, foundererror, offspringerror, baseerror, 
        #       allelicbias, allelicoverdispersion, allelicdropout
        resinfo = MagicBase.parse_vcf_info(info)
        res = vcat([snpid,resinfo[2],resinfo[3],physchrom, physposbp,resinfo[1]], ["NA", "NA"], resinfo[4:end])
        push!(df, res)
    end    
    misscodes = isa(missingstring,AbstractString) ? [missingstring] : missingstring
    in(misscodes,"") || push!(misscodes,"")
    df[!,:poscm] .= [in(i,misscodes) ? missing : parse(Float64,i) for i in df[!,:poscm]]
    df[!,:physposbp] .= [in(i,misscodes)  ? missing : parse(Int,i) for i in df[!,:physposbp]]    
    for col in vcat([2,4],6:8)
        df[!,col] .= [in(i,misscodes) ? missing : i for i in df[!,col]]
    end
    for col in 9:length(colnames)
        df[!,col] .= [in(i,misscodes) ? missing : parse(Float64,i) for i in df[!,col]]
    end
    df
end



function readgenofile(genofile::AbstractString;
    commentstring::AbstractString="##",
    missingstring=["NA","missing"],
    workdir::AbstractString=pwd())
    ext = last(split_allext(genofile))
    in(ext, [".csv",".csv.gz",".vcf",".vcf.gz"]) || @error "genofile ext must be .csv or .csv.gz"    
    isvcf = ext in [".vcf",".vcf.gz"]
    delim = isvcf ? '\t' : ','
    genofile2=getabsfile(workdir,genofile)
    isfile(genofile2) || @error string(genofile2, " does not exit")
    lastcomment = findlastcomment(genofile2; commentstring)
    commentlines = read_headlines(genofile2, lastcomment)
    if in(ext, [".csv.gz",".vcf.gz"])
        genodf = readdataframe(genofile2; delim, missingstring, commentstring)
    else
        genodf = CSV.read(genofile2,DataFrame; delim, missingstring, comment=commentstring)
    end
    genodf, commentlines
end


"""
    vcf_extract_pedfile(vcffile; keyargs...)

extract pedfile from vcffile. Work only for a non-subdivided population.  

# Positional arguments

`vcffile::AbstractString`: vcf genofile. 

# Keyword arguments

`designcode::AbstractString`: designcode

`ishomozygous::Bool=false`: specify if offspring are completely homozygous.

`isfglexch::Bool=true`: specify if founders are exchangeable.

`popid::AbstractString="pop"`: population id. 

`outstem::AbstractString="outstem"`: stem of output filename.

`workdir::AbstractString=pwd()`: directory for reading and saving files.

"""
function vcf_extract_pedfile(vcffile::AbstractString;         
    designcode::AbstractString,    
    ishomozygous::Bool=false,
    isfglexch::Bool=true,
    popid::AbstractString = "pop",
    outstem::AbstractString= "outstem",
    workdir::AbstractString=pwd())
    samples = MagicBase.vcf_get_samples(vcffile; workdir)
    designinfo = parsedesign(designcode; popid)
    if in(designinfo.designtype,[:commoncross,:juncdist])
        nfounder = getnfounder(designinfo)    
        length(samples) >= nfounder || @error string("#indviiduals in arrayfiles < nfounder by designcode=",designcode)
        founders = samples[1:nfounder]
        setfounders!(designinfo,founders)
    elseif designinfo.designtype == :breedcross
        founders = designinfo.founders
        d = setdiff(founders,samples)
        isempty(d) || @error string("designcode founders=",d, " are not in arrayfile")        
    else
        @error string("unknown designtype=",designtype)
    end    
    @info string("extracted founders = ",founders)        
    d = setdiff(founders,samples)
    isempty(d) || @error string("founders that are not genotyped in vcffile: ",d)
    pedfile = getabsfile(workdir,outstem*"_ped.csv")
    open(pedfile,"w") do io
        write(io,"RABBIT,designinfo\n")
        write(io, "member,founders,designcode\n")
        member = popid
        founders2 = designinfo.designtype == :breedcross ? "NA" : join(founders,"||")
        write(io, join([member, founders2, designcode],","),"\n")
        write(io,"RABBIT,offspringinfo\n")        
        offls = setdiff(samples, founders)
        offinfo = DataFrame(individual=offls,member=member,ishomozygous=ishomozygous,isfglexch=isfglexch)
        @info offinfo[1:min(10,size(offinfo,1)),:]
        CSV.write(io, offinfo; append=true,writeheader=true)
    end
    pedfile
end


function vcf_plink_map(vcffile::AbstractString;    
    isphysmap::Bool=true,
    recomrate::Real=1, # 1 cm?Mbp
    commentstring = "##",
    outstem::AbstractString="outstem",
    logfile::Union{Nothing,AbstractString} = nothing,
    workdir::AbstractString=pwd(),
    verbose::Bool=true)    
    starttime = time()
    logio = MagicBase.set_logfile_begin(logfile, workdir, "vcf_plink_map"; verbose,delim="-")    
    msg = string("list of args: \n",
        "vcffile = ", vcffile, "\n",                
        "isphysmap = ", isphysmap, "\n",                
        "recomrate = ", recomrate, " cM/Mbp \n",                

        "commentstring = ", commentstring, "\n",
        "outstem = ", outstem, "\n",        
        "logfile = ", logfile, "\n",
        "workdir = ", workdir, "\n",
        "verbose = ", verbose)
    printconsole(logio,verbose,msg)
    inext = last(MagicBase.split_allext(vcffile))
    inext in [".vcf",".vcf.gz"] || @error string("vcffile ext must be .vcf or .vcf.gz")
    in_open = inext in [".vcf.gz"] ? GZip.open : open    
    lastcomment = MagicBase.findlastcomment(vcffile; commentstring)
    nheader = lastcomment + 1
    # The PLINK MAP file must has exactly four columns with the following information (the columns should be separated by a whitespace or a tab):
    # Chromosome ID (e.g. Chr1 for Chromosome 1)
    # Unique SNP identifier
    # Genomic distance (if unknown use 0)
    # SNP Position
    # if there exists genetic map in vcf, chrid = linkageroup, and otherwise chrid = CHROM
    outfile = getabsfile(workdir,outstem*"_plink_map.tsv")
    in_open(getabsfile(workdir,vcffile),"r") do io
        open(outfile, "w") do outio
            for _ in 1:nheader-1
                readline(io)
            end
            # skip title row                
            readline(io)
            # parse data
            while !eof(io)                
                line = split(readline(io,keep=false),"\t")                
                # vcfcols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL","FILTER", "INFO", "FORMAT"]
                # resinfo: [newinfo, linkagegroup, poscm, foundererror, offspringerror, baseerror, allelicbias, allelicoverdispersion, allelicdropout]
                resinfo = MagicBase.parse_vcf_info(line[8])
                if in("NA", resinfo[2:3])                 
                    chrid = line[3]   
                    if isphysmap && !in(".",line[1:2])
                        posbp = parse(Int, line[2])
                        poscm = posbp * 1e-6 * recomrate
                    else
                        poscm = "."
                    end
                else
                    chrid, poscm = resinfo[2:3]                    
                end
                line2 = [line[1],chrid, poscm, line[2]]
                write(outio, join(line2,"\t"),"\n")                
            end
        end
    end
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, logio, tused,"vcf_plink_map"; verbose,delim="-")    
    outfile    
end


function plink_map(mapfile::AbstractString;
    missingstring=["NA","missing"],
    outstem::AbstractString = "outstem",     
    commentstring::AbstractString="##",
    workdir::AbstractString = pwd())
    mapdf = readmarkermap(mapfile; commentstring,missingstring, workdir, del_ungrouped = true)
    # The PLINK MAP file must has exactly four columns with the following information (the columns should be separated by a whitespace or a tab):
    # Chromosome ID (e.g. Chr1 for Chromosome 1)
    # Unique SNP identifier
    # Genomic distance (if unknown use 0)
    # SNP Position
    # if there exists genetic map in vcf, chrid = linkageroup, and otherwise chrid = CHROM
    if all(ismissing.(mapdf[!,:linkagegroup]))
        res = mapdf[!,[4,1,3,5]]
    else
        res = mapdf[!,[2,1,3,5]]
    end
    outfile = getabsfile(workdir,outstem*"_plink_map.tsv")
    CSV.write(outfile,res; header=false,delim='\t',missingstring=".")
end