
# output arraygenofile: called genotypic data matrix: 
#   first row is sample id, first column is marker id, 
#   entries are called genotypes, one of the possible pairwise combinatons of A, T, G, and C (- denotes missing allele)
function toarrayfile(rawfile::AbstractString; 
    delim::AbstractChar = '\t',
    genoformat::AbstractString = "Top Alleles",
    outstem::AbstractString  = "outstem",
    outext::AbstractString = ".csv.gz", 
    workdir::AbstractString = pwd(),
    logfile::Union{Nothing,AbstractString,IO}= nothing,
    verbose::Bool=true)
    starttime = time()
    logio = MagicBase.set_logfile_begin(logfile, workdir, "toarrayfile"; verbose)   
    in(outext,[".csv",".csv.gz"]) || @error string("outext=",outext, " is not in [.csv,.csv.gz]")
    msg = string("list of args: \n",
        "rawfile = ", rawfile, "\n",        
        "delim = ", string(delim), "\n",
        "genoformat = ", genoformat, "\n",
        "workdir = ",workdir,"\n",        
        "outstem = ", isnothing(outstem) ? "no output files" : outstem,"\n",
        "outext = ",outext,"\n",
        "logfile = ",logfile,"\n",
        "verbose = ",verbose)
    printconsole(logio,verbose,msg)    
    rawfile2 = getabsfile(workdir,rawfile)
    in_open = last(splitext(rawfile)) == ".gz" ? GZip.open : open
    nmarker = 0
    nind = [0]
    in_open(rawfile2,"r") do inio    
        outfile = getabsfile(workdir,outstem*outext)
        out_open = outext == ".csv.gz" ? GZip.open : open
        out_open(outfile, "w") do outio            
            title = split(readline(inio;keep=false),delim)
            msg = string("original ", length(title), " columns (print upto 10): ", join(title[1:min(10,length(title))],","))
            printconsole(logio, verbose, msg)
            cols = findall([occursin(genoformat,i) for i in title])
            pushfirst!(cols,1)
            title = replace.(title[cols],"."*genoformat =>"","_"*genoformat =>"", genoformat =>"")
            nind[1] = length(title) -1
            msg = string("new ", length(title), " columns (print upto 10): ", join(title[1:min(10,length(title))],","),
                " corresponding to oriingal cols=",cols[1:min(10,length(cols))])
            printconsole(logio, verbose, msg)
            write(outio, join(title,","),"\n")
            while !eof(inio)
                nmarker += 1
                line = split(readline(inio;keep=false),delim)
                write(outio, join(line[cols],","),"\n")
            end
        end        
        msg = string("save arraygenofile: ", outfile, "")
        printconsole(logio,verbose,msg)
        msg = string("#markers=",nmarker, ", #individuals=",nind[1])
        printconsole(logio,verbose,msg)
        tused = round(time()-starttime,digits=1)
        MagicBase.set_logfile_end(logfile, logio, tused,"toarrayfile"; verbose)
        outfile
    end
end

"""
    arrayfile2vcf(arrayfile; keyargs...)

extract pedfile from arrayfile. Work only for a non-subdivided population.  

# Positional arguments

`arrayfile::AbstractString`: SNP array genofile. 

# Keyword arguments

`delmultiallelic::Bool = true,`: if true, delete markers with #alleles >= 3. 

`delim= ","`: text delimiter. 

`outstem::AbstractString="outstem"`: stem of output filename.

`workdir::AbstractString=pwd()`: directory for reading and saving files.

`logfile::Union{Nothing,AbstractString,IO}= outstem*"_arrayfile2vcf.log"`: log file or IO for writing log. If it is nothing, no log file.

`verbose::Bool=true`: if true, print details on the stdout.

"""
function arrayfile2vcf(arrayfile::AbstractString; 
    missingallele = "-", 
    delmultiallelic::Bool = true,
    delim = ',',        
    outstem::AbstractString = first(MagicBase.split_allext(basename(arrayfile))),
    workdir::AbstractString = pwd(),
    logfile::Union{Nothing,AbstractString,IO}= outstem*"_arrayfile2vcf.log",
    verbose::Bool=true)
    starttime = time()
    logio = MagicBase.set_logfile_begin(logfile, workdir, "arrayfile2vcf"; verbose)       
    msg = string("list of args: \n",
        "arrayfile = ", arrayfile, "\n",        
        "missingallele = ", missingallele, "\n",
        "delim = ", string(delim), "\n",        
        "workdir = ",workdir,"\n",        
        "outstem = ", isnothing(outstem) ? "no output files" : outstem,"\n",        
        "logfile = ",logfile,"\n",
        "verbose = ",verbose)
    printconsole(logio,verbose,msg)        
    arrayfile2 = getabsfile(workdir,arrayfile)
    in_open = last(splitext(arrayfile2)) == ".gz" ? GZip.open : open
    startt = time()
    nmarker = 0
    nmulti = 0
    nsample = [0]
    outfile = getabsfile(workdir,outstem*"_geno.vcf.gz")    
    in_open(arrayfile2,"r") do inio    
        GZip.open(outfile, "w") do outio    
            title = split(readline(inio;keep=false),delim)
            popfirst!(title) # 1st column is markerid
            nsample[1] = length(title)
            msg = string(nsample[1], " genotyped individuals (print upto 10): ", join(title[2:min(11,nsample[1])],","))
            printconsole(logio, verbose, msg)
            title =vcat(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL","FILTER", "INFO", "FORMAT"],title)
            write(outio, join(title,'\t'),"\n")     
            while !eof(inio)    
                nmarker += 1
                rowgeno = split(readline(inio,keep=false),delim)
                markerid = popfirst!(rowgeno)
                alleles = replace(join(rowgeno),missingallele=>"")
                aset = unique(alleles)
                if length(aset) >= 3 
                    nmulti += 1
                    delmultiallelic && continue
                end
                if isempty(aset)
                    ref = "N"
                    alt = "."
                    rowgeno .= "./."
                else
                    afreq = [count(i, alleles) for i in aset]
                    ref = aset[argmax(afreq)]
                    setdiff!(aset,[ref])
                    if isempty(aset)
                        alt = "."
                        rowgeno .= replace.(rowgeno, missingallele=>".", ref=>"0")
                    elseif length(aset) == 1
                        alt = only(aset)
                        rowgeno .= replace.(rowgeno, missingallele=>".", ref=>"0",alt=>"1")
                    else
                        alt = join(sort(aset))
                        rowgeno .= replace.(rowgeno, missingallele=>".", ref=>"0", (aset .=> "1")...)
                    end    
                    rowgeno .= [join(split(i,""),"/") for i in rowgeno]
                end
                col19 = [".",".",markerid,ref,alt,".",".",".","GT"]
                write(outio, join(col19,"\t"), "\t", join(rowgeno,"\t"),"\n")
                if rem(nmarker, 1000) == 0
                    msg = string("marker=", nmarker,", tused=", round(time()-startt,digits=1),"s")
                    printconsole(logio,verbose,msg)
                    flush(outio)
                end
            end
        end
    end
    msg = string("#individuals=",nsample[1],", #markers=",nmarker)
    if nmulti > 0
        msg *= string(", #markers_multiallelic=",nmulti, delmultiallelic ? "(deleted)" : "" )
    end
    printconsole(logio,verbose,msg)    
    msg = string("ouput vcffile: ", outfile)
    printconsole(logio,verbose,msg)    
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, logio, tused,"arrayfile2vcf"; verbose)        
    outfile
end

function array_get_samples(arrayfile::AbstractString; 
    delim=",",
    workdir::AbstractString=pwd())
    inext = last(MagicBase.split_allext(arrayfile))
    inext in [".csv",".csv.gz"] || @error string("arrayfile ext must be .csv or .csv.gz")
    in_open = inext in [".csv.gz"] ? GZip.open : open        
    arrayfile2 = getabsfile(workdir,arrayfile)
    in_open(arrayfile2,"r") do io        
        # parse title row
        string.(split(readline(io,keep=false),delim)[2:end])
    end
end


function array_rename_samples(arrayfile::AbstractString; 
    old_new::AbstractVector,
    delim::AbstractString=",", 
    outstem::AbstractString="outstem",
    logfile::Union{Nothing,IO,AbstractString} = nothing,
    workdir::AbstractString=pwd(),
    verbose::Bool=true)    
    renamedict = old_new 
    if isempty(renamedict)
        @warn "empty old_new"
        return nothing
    end
    starttime = time()
    logio = MagicBase.set_logfile_begin(logfile, workdir, "array_rename_samples"; verbose,delim="-")    
    msg = string("list of args: \n",
        "arrayfile = ", arrayfile, "\n",        
        "old_new (upto 10)= ", old_new[1:min(10,length(old_new))], "\n",                
        "outstem = ", outstem, "\n",        
        "logfile = ", logfile, "\n",
        "workdir = ", workdir, "\n",
        "verbose = ", verbose)
    printconsole(logio,verbose,msg)
    inext = last(MagicBase.split_allext(arrayfile))
    inext in [".csv",".csv.gz"] || @error string("arrayfile ext must be .csv or .csv.gz")
    in_open = inext in [".csv.gz"] ? GZip.open : open    
    outfile = getabsfile(workdir,outstem*".csv.gz")
    in_open(getabsfile(workdir,arrayfile),"r") do io
        GZip.open(outfile, "w") do outio            
            # parse title row    
            titlerow = split(readline(io,keep=false),delim)
            samples = string.(strip.(titlerow[2:end]))
            dict = Dict(samples .=> 1:length(samples))
            col_new = [[get(dict,old,nothing),old, new] for (old, new) in renamedict]
            cols = first.(col_new)
            news = [i[3] for i in col_new]
            b = isnothing.(cols)
            if any(b)
                msg = string("samples in old_new but not in arrayfile: ",[i[2] for i in col_new[b]])                
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
                # msg = string("rename samples: ", samples[cols] .=> news)
                # printconsole(logio, verbose, msg)
                cols .+= 1
                titlerow[cols] .= news            
            end
            write(outio, join(titlerow,delim),"\n")
            # parse data
            while !eof(io)
                line = readline(io;keep=true)
                write(outio, line)
            end
        end
    end
    msg = string("save samples-renamed arrayfile: ", outfile)
    printconsole(logio, verbose,msg)
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, logio, tused,"array_rename_samples"; verbose,delim="-")    
    outfile    
end

function array_get_subgeno(arrayfile::AbstractString; 
    subsamples::AbstractVector=[],    
    delim=",",
    workdir::AbstractString=pwd())
    inext = last(MagicBase.split_allext(arrayfile))
    inext in [".csv",".csv.gz"] || @error string("arrayfile ext must be .csv or .csv.gz")
    in_open = inext in [".csv.gz"] ? GZip.open : open 
    arrayfile2 = getabsfile(workdir,arrayfile)           
    in_open(arrayfile2,"r") do io                
        # parse title row    
        titlerow = split(readline(io,keep=false),delim)
        length(titlerow) <=1 && @warn string("delim=",delim, " mismatches the title row=",titlerow)
        samples = strip.(titlerow[2:end])
        if isempty(subsamples)
            subsamples2 = []
        elseif eltype(subsamples) <: Integer 
            subsamples2 = samples[subsamples]
        elseif eltype(subsamples) <: AbstractString
            subsamples2 = string.(strip.(subsamples))
            d = setdiff(subsamples2, samples)
            if !isempty(d)
                @warn string("subsamples that are not in arrayfile: ",d)
            end
            intersect!(subsamples2,samples)
        else
            @error string("unknown type=",eltype(subsamples), " for input subsamples=",subsamples)
        end           
        lefttitles = ["marker"]
        leftcols= [1]
        leftpairs = [Symbol(i)=>String[] for i in lefttitles[leftcols]]
        if isempty(subsamples2)
            resdf = DataFrame(leftpairs...)
            cols = leftcols
        else
            resdf = DataFrame(leftpairs..., [Symbol(i) => String[] for i in subsamples2]...)            
            cols = [findfirst(==(i),samples) + 1 for i in subsamples2]
            pushfirst!(cols,leftcols...) 
        end        
        while !eof(io)            
            line = split(readline(io,keep=false),delim)
            push!(resdf, line[cols])            
        end
        resdf
    end
end


function array_recode_geno(arrayfile::AbstractString; 
    allelefreqfile::AbstractString, 
    delim::AbstractChar = ',',        
    outstem::AbstractString = first(MagicBase.split_allext(basename(arrayfile)))*"_recode",    
    outext::AbstractString = ".csv.gz", 
    logfile::Union{Nothing,AbstractString,IO}= nothing,
    workdir::AbstractString = pwd(),    
    verbose::Bool=true)
    starttime = time()
    logio = MagicBase.set_logfile_begin(logfile, workdir, "array_recode_geno"; verbose,delim="-")    
    msg = string("list of args: \n",
        "arrayfile = ", arrayfile, "\n",        
        "allelefreqfile = ",allelefreqfile, "\n",         
        "delim = ",delim, "\n", 
        "outstem = ", outstem, "\n",        
        "logfile = ", logfile, "\n",
        "workdir = ", workdir, "\n",
        "verbose = ", verbose)
    printconsole(logio,verbose,msg)
    freq = CSV.read(getabsfile(workdir, allelefreqfile),DataFrame)
    newcoding = Dict(freq[!,:marker] .=> tuple.(freq[!,"allele_major"],freq[!,"allele_alter"]))
    startt = time()
    nmarker = 0
    nind = [0]
    missingallele = '-'    
    inext = last(MagicBase.split_allext(arrayfile))
    inext in [".csv",".csv.gz"] || @error string("arrayfile ext must be .csv or .csv.gz")
    in_open = inext in [".csv.gz"] ? GZip.open : open        
    outext in [".csv",".csv.gz"] || @error string("outext ext must be .csv or .csv.gz")
    out_open = outext in [".csv.gz"] ? GZip.open : open    
    outfile = outstem*outext
    outfile2 = getabsfile(workdir,outfile)
    in_open(getabsfile(workdir,arrayfile),"r") do inio
        out_open(outfile2, "w") do outio            
            # parse title row    
            title0 = readline(inio;keep=false)
            title = split(title0,delim)
            msg = string(length(title), " columns (print upto 10): ", join(title[1:min(10,length(title))],","))
            nind[1] = length(title) - 1 
            printconsole(logio, verbose, msg)                        
            write(outio, title0, "\n")
            # parse data            
            while !eof(inio)                
                nmarker += 1
                line0 = readline(inio;keep=false)
                line = split(line0,delim)                
                marker = first(line)               
                if !haskey(newcoding,marker)
                  write(outio, line0,"\n")
                  continue
                end
                newmajor, newalter = newcoding[marker]  
                if ismissing(newmajor)
                    write(outio, line0,"\n")
                    continue
                end
                if !ismissing(newalter) && length(newalter) > 1
                    msg = string("Not recode genotypes at marker = ", marker, " with mulitple allele_alter  = ", newalter, " in the allelefreqfile")
                    printconsole(logio, false, "Warning: "*msg)
                    @warn msg
                    write(outio, line0,"\n")
                    continue
                end           
                alleles = replace(join(line[2:end]),missingallele =>"")    
                aset = unique(alleles)
                nallele = length(aset)
                if nallele == 0 
                    write(outio, line0,"\n")
                    continue
                end
                if nallele > 2
                    msg = string("Not recode genotypes at marker = ", marker, " with mulitple alleles = ", aset, " in the arrayfile")
                    printconsole(logio, false, "Warning: "*msg)
                    @warn msg
                    write(outio, line0,"\n")
                    continue
                end
                if nallele == 2 && count(aset[1],alleles)/length(alleles) < 0.5
                    reverse!(aset)
                end
                if ismissing(newalter)
                    if in(newmajor, aset)
                        write(outio, line0,"\n")
                    else
                        line[2:end] .= replace.(line[2:end],aset[1]=>newmajor)
                        write(outio, join(line,delim),"\n")
                    end
                else         
                    if issubset(aset, [newmajor, newalter])    
                        write(outio, join(line,delim),"\n")
                    else
                        if nallele == 1    
                            line[2:end] .= replace.(line[2:end],aset[1]=>newmajor)
                            write(outio, join(line,delim),"\n")
                        elseif nallele == 2                            
                            line[2:end] .= replace.(line[2:end],aset[1]=>newmajor,aset[2]=>newalter)
                            write(outio, join(line,delim),"\n")
                        else
                            msg = string("unexpected alleles = ", aset, " at marker = ", marker, " in the allelefreqfile")
                            printconsole(logio, false, "Warning: "*msg)
                            @warn msg
                        end                    
                    end
                end                
                if rem(nmarker, 1000) == 0
                    msg = string("marker=", nmarker,", tused=", round(time()-startt,digits=1),"s")
                    printconsole(logio,verbose,msg)
                    flush(outio)
                end
            end            
        end
    end
    msg = string("#individuals=", nind[1], ", #markers=",nmarker)
    printconsole(logio,verbose,msg)
    msg = string("save recode arrayfile: ", outfile)
    printconsole(logio, verbose,msg)
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, logio, tused,"array_recode_geno"; verbose,delim="-")    
    outfile  
end


function array_extract_allelefreq(arrayfile::AbstractString; 
    delim::AbstractChar = ',',            
    outstem::AbstractString = first(MagicBase.split_allext(basename(arrayfile)))*"_allelefreq",    
    workdir::AbstractString = pwd(),
    logfile::Union{Nothing,AbstractString,IO}= nothing,
    verbose::Bool=true)
    # input arrayfile: 
    #   first row is sample id, first column is marker id, 
    #   entries are called genotypes, one of the possible pairwise combinatons of A, T, G, and C (- denotes missing allele
    starttime = time()
    logio = MagicBase.set_logfile_begin(logfile, workdir, "array_extract_allelefreq"; verbose)               
    msg = string("list of args: \n",
        "arrayfile = ", arrayfile, "\n",        
        "delim = ", string(delim), "\n",        
        "workdir = ",workdir,"\n",        
        "outstem = ", outstem, "\n",
        "logfile = ",logfile,"\n",
        "verbose = ",verbose)
    printconsole(logio,verbose,msg)    
    arrayfile2 = getabsfile(workdir,arrayfile)
    in_open = last(splitext(arrayfile)) == ".gz" ? GZip.open : open
    outfile = outstem*".csv"
    nmarker = 0
    nind = [0]
    missingallele = '-'
    in_open(arrayfile2,"r") do inio    
        outfile2 = getabsfile(workdir,outfile)
        open(outfile2, "w") do outio            
            title = split(readline(inio;keep=false),delim)
            msg = string(length(title), " columns (print upto 10): ", join(title[1:min(10,length(title))],","))
            nind[1] = length(title) - 1 
            printconsole(logio, verbose, msg)                        
            write(outio, "marker,allele_major,allele_alter,freq_major,freq_alter,allele_count\n")
            while !eof(inio)
                nmarker += 1
                line = split(readline(inio;keep=false),delim)                
                marker = first(line)
                alleles = replace(join(line[2:end]),missingallele =>"")
                write(outio, marker,",")
                if isempty(alleles)
                    write(outio, ",,,, 0\n")
                else
                    aset = unique(alleles)
                    acount = [count(i,alleles) for i in aset]
                    cc, pos = findmax(acount)
                    tot = sum(acount)
                    allele_major = aset[pos]
                    freq_major = round(cc/tot,digits=6)
                    deleteat!(aset, pos)
                    deleteat!(acount, pos)
                    allele_alter = join(aset,"|")
                    freq_alter = join(round.(acount ./ tot,digits=6), "|")
                    write(outio, join([allele_major, allele_alter, freq_major,freq_alter, tot],","), "\n")
                end
            end
        end        
        msg = string("save allelefreq: ", outfile, "")
        printconsole(logio,verbose,msg)
        msg = string("#markers=",nmarker, ", #individuals=",nind[1])
        printconsole(logio,verbose,msg)
        tused = round(time()-starttime,digits=1)
        MagicBase.set_logfile_end(logfile, logio, tused,"toarrayfile"; verbose)
        outfile
    end
end


"""
    array_extract_pedfile(arrayfile; keyargs...)

extract pedfile from arrayfile. Work only for a non-subdivided population.  

# Positional arguments

`arrayfile::AbstractString`: SNP array genofile. 

# Keyword arguments

`designcode::AbstractString`: designcode

`popid::AbstractString`: population id. 

`delim= ","`: text delimiter. 

`outstem::AbstractString=popid`: stem of output filename.

`workdir::AbstractString=pwd()`: directory for reading and saving files.
"""
function array_extract_pedfile(arrayfile::AbstractString;         
    designcode::AbstractString,    
    popid::AbstractString,
    delim= ",",    
    outstem::AbstractString = first(MagicBase.split_allext(basename(arrayfile))),
    workdir::AbstractString=pwd())
    samples = MagicBase.array_get_samples(arrayfile; delim)
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
    isempty(d) || @error string("founders that are not genotyped in arrayfile: ",d)
    pedfile = getabsfile(workdir,outstem*"_ped.csv")
    open(pedfile,"w") do io
        write(io,"RABBIT,designinfo\n")
        member = popid
        if designinfo.designtype == :breedcross
            write(io, "member,designcode\n")
            write(io, join([member, designcode],","),"\n")
        else
            write(io, "member,designcode,founders\n")
            write(io, join([member, designcode, join(founders,"||")],","),"\n")
        end        
        write(io,"RABBIT,offspringinfo\n")        
        offls = setdiff(samples, founders)
        if occursin(r"=>DH$",designcode) || occursin(r"=>FIXED$",designcode)
            ishomo = true
        elseif designinfo.designtype == :juncdist 
            ishomo = designinfo.juncdist.ibd == 1.0        
        else 
            ishomo = false
        end
        isfglexch = designinfo.designtype == :juncdist         
        offinfo = DataFrame(individual=offls,member=member)
        ishomo && insertcols!(offinfo, size(offinfo,2)+1,:ishomozygous=>ishomo)        
        isfglexch && insertcols!(offinfo, size(offinfo,2)+1,:isfglexch=>isfglexch)        
        @info offinfo[1:min(10,size(offinfo,1)),:]
        CSV.write(io, offinfo; append=true,writeheader=true)
    end
    pedfile
end

"""
    merge_arrayfiles(arrayfiles; outstem, workdir)

merge arrayfiles into a single SNP array genofile. 

# Positional arguments

`arrayfiles::AbstractVector`: a list of SNP array genofile. 

# Keyword arguments

`missingallele::AbstractString = "-"`: string for missing allele. 

`outstem::AbstractString="outstem"`: stem of output filename.

`outext::AbstractString=".csv.gz"`: extension of output file.

`workdir::AbstractString=pwd()`: directory for reading and saving files.

`logfile::Union{Nothing,AbstractString,IO}= outstem*"_merge_arrayfiles.log"`: log file or IO for writing log. If it is nothing, no log file.

`verbose::Bool=true`: if true, print details on the stdout.
"""
function merge_arrayfiles(arrayfiles::AbstractVector;
    missingallele::AbstractString = "-", 
    outstem::AbstractString = "outstem",    
    outext::AbstractString = ".csv.gz", 
    workdir::AbstractString = pwd(),
    logfile::Union{Nothing,AbstractString,IO}= outstem*"_merge_arrayfiles.log",
    verbose::Bool=true)
    starttime = time()
    logio = MagicBase.set_logfile_begin(logfile, workdir, "merge_arrayfiles"; verbose)   
    in(outext,[".csv",".csv.gz"]) || @error string("outext=",outext, " is not in [.csv,.csv.gz]")
    msg = string("list of args: \n",
        "arrayfiles = ", arrayfiles, "\n",                                
        "outstem = ", isnothing(outstem) ? "no output files" : outstem,"\n",
        "outext = ",outext,"\n",
        "workdir = ",workdir,"\n",        
        "logfile = ",logfile,"\n",
        "verbose = ",verbose)
    printconsole(logio,verbose,msg)        
    arrayfilels = [getabsfile(workdir,i) for i in arrayfiles]
    genodf = MagicBase.readdataframe(arrayfilels[1])
    rename!(genodf,1=>:marker)
    allunique(genodf[!,:marker]) || @warn "markers are not unique"    
    popfirst!(arrayfilels)
    missingcode = missingallele^2     
    for arrayfile in arrayfilels       
        genodf2 = MagicBase.readdataframe(arrayfile)
        msg = string("accumulated #markers=",size(genodf,1), ", #individuals=",size(genodf,2)-1) 
        msg *= string(", start merging arrayfile=",arrayfile)
        msg *= string(" with #markers=",size(genodf2,1), ", #individuals=",size(genodf2,2)-1) 
        printconsole(logio,verbose,msg)                        
        rename!(genodf2,1=>:marker)
        allunique(genodf2[!,:marker]) || @warn "markers are not unique"  
        # combine  columns for duplicated samples
        dupls = setdiff(intersect(names(genodf),names(genodf2)),["marker"])
        genodf = outerjoin(genodf,genodf2; on = :marker,makeunique=true)      
        samplels = names(genodf)[2:end]
        for c in dupls            
            cls = filter(x->occursin(c,x), samplels)            
            intersect!(cls,[string(c,"_",i) for i in 1:length(cls)]) # resulting from makeunique=true            
            setdiff!(cls,[c])
            if !isempty(cls)                
                geno_cls =  Matrix(genodf[!,vcat(cls,[c])])
                geno_c = [setdiff(skipmissing(unique(i)),[missingcode]) for i in eachrow(geno_cls)]
                genodf[!,c] .= [length(i)==1 ? i[1] : missing for i in geno_c]
                select!(genodf,Not(cls))
                ninconsist = sum(length.(geno_c) .> 1)
                msg = string("merging samples: ",vcat(cls,[c]), ", #inconsistent=", ninconsist, " (set to missing)")
                printconsole(logio,verbose,msg)        
            end
        end
    end    
    outfile = getabsfile(workdir,outstem*"_array"*outext)
    out_open = last(splitext(outext)) == ".gz" ? GZip.open : open
    out_open(outfile,"w") do file        
        CSV.write(file,genodf; 
            writeheader=true,
            append=true, 
            missingstring=missingcode,
            delim = ","
        )
    end
    msg = string("accumulated #markers=",size(genodf,1), ", #individuals=",size(genodf,2)-1)
    printconsole(logio,verbose,msg)    
    msg = string("save merged arrayfile: ", outfile, "")
    printconsole(logio,verbose,msg)    
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, logio, tused,"merge_arrayfiles"; verbose)        
    outfile
end

