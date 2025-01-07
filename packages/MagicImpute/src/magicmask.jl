
function magicmask(genofile::AbstractString,pedinfo::Union{Integer,MagicBase.JuncDist, AbstractString};
    isfounderinbred::Bool=true,
    formatpriority::AbstractVector=["AD","GT"],
    isphysmap::Bool=false,
    recomrate::Real=1.0,
    foundermask::Real = 0.1,
    offspringmask::Real = 0.1,
    skipmarker::Real = 0.99, #skip masking genotypes if offspring missing fraction >= skipmarker    
    minread::Integer = 10,
    commentstring::AbstractString="##",
    outstem::AbstractString= "outstem",
    outext::AbstractString=in(last(MagicBase.split_allext(genofile)),[".vcf",".vcf.gz"]) ? ".vcf.gz" : ".csv.gz",
    logfile::Union{AbstractString,IO} = outstem*"_magicmask.log",
    workdir::AbstractString=pwd(),
    verbose::Bool=true)
    starttime = time()
    logio = MagicBase.set_logfile_begin(logfile, workdir, "magicmask"; verbose,delim="=")
    MagicBase.check_infiles(genofile,pedinfo; isbreedped=false, io=logio, commentstring,workdir,verbose)    
    info_file_arg(genofile, pedinfo, formatpriority,isphysmap, recomrate,
        commentstring,workdir, logio, verbose)
    MagicBase.printconsole(logio,verbose,string("nthreads=",Threads.nthreads()))
    check_mask_arg(foundermask, offspringmask, minread; outext)
    tused = @elapsed magicgeno = formmagicgeno(genofile, pedinfo;
        isfounderinbred, isphysmap, recomrate,
        formatpriority, commentstring, missingstring=["NA","missing"], workdir)
    msg = string("tused=", round(tused,digits=1), " seconds by formmagicgeno")
    MagicBase.printconsole(logio,verbose,msg)
    magicmask!(magicgeno;
        foundermask, offspringmask, minread,
        skipmarker,
        outstem,outext,logfile=logio, workdir, verbose)
    outfile = string(outstem,"_magicmask_geno")*outext
    savegenodata(outfile,magicgeno; commentstring,workdir)
    msg = string("save magicmask genofile: ", outfile)
    MagicBase.printconsole(logio,verbose,msg)
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, logio, tused,"magicmask"; verbose,delim="=")
    magicgeno
end

function magicmask!(magicgeno::MagicGeno;
    foundermask::Real = 0.0,
    offspringmask::Real = 0.1,    
    minread::Integer = 10,
    skipmarker::Real = 0.99,
    commentstring::AbstractString="##",
    outstem::AbstractString= "outstem",
    outext::AbstractString=".vcf.gz",
    logfile::Union{AbstractString,IO} = outstem*"_magicmask.log",
    workdir::AbstractString=pwd(),
    verbose::Bool=true)
    starttime = time()
    logio = MagicBase.set_logfile_begin(logfile, workdir, "magicmask!"; verbose)
    check_mask_arg(foundermask, offspringmask, minread)
    msg = string("list of options: \n",
        "foundermask = ", foundermask, "\n",
        "offspringmask = ", offspringmask, "\n",        
        "minread = ", minread, "\n",
        "skipmarker = ", skipmarker, "\n",
        "workdir = ",workdir,"\n",
        "outstem = ", outstem,"\n",
        "logfile = ",logfile,"\n",
        "verbose = ",verbose)
    MagicBase.printconsole(logio,verbose,msg)    
    MagicBase.info_magicgeno(magicgeno;io=logio,verbose)    
    # deepcopy for a large size of magicgeno is very slow
    # tused = @elapsed undermask = deepcopy(magicgeno) 
    tused = @elapsed undermask = MagicGeno(deepcopy(magicgeno.magicped),
        copy.(magicgeno.markermap), copy.(magicgeno.foundergeno),
        copy.(magicgeno.offspringgeno),deepcopy(magicgeno.misc))
    msg = string("tused=",round(tused,digits=1),"s in copying magicgeno")
    printconsole(logio,verbose,msg)
    # masking founders
    nmask,nnonmiss = mask_founder_offspring!(magicgeno,undermask;
        minread, skipmarker,
        ismaskfounder=true,foundermask,offspringmask,
        io = logio,verbose)
    # masking offspring
    nmask,nnonmiss = mask_founder_offspring!(magicgeno,undermask;
        minread, skipmarker,
        ismaskfounder=false,foundermask,offspringmask,
        io = logio,verbose)    
    outfile =string(outstem,"_magicmask_reversed",outext)
    savegenodata(outfile,undermask; commentstring,workdir)
    msg = string("save reversed masked-genofile (ground truth): ", outfile)
    MagicBase.printconsole(logio,verbose,msg)
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, logio, tused,"magicmask!"; verbose)
    magicgeno
end

function info_file_arg(genofile::AbstractString,
    pedinfo::Union{Integer, MagicBase.JuncDist, AbstractString},
    formatpriority::AbstractVector,
    isphysmap::Bool,
    recomrate::Real,
    commentstring::AbstractString,
    workdir::AbstractString,
    io::Union{Nothing,IO},
    verbose::Bool)
    if !isdir(workdir)
        @error string(workdir, " is not a directory")
    end
    genofile2 = getabsfile(workdir, genofile)
    if !isfile(genofile2)
        @error string(genofile2, " does not exist")
    end
    if isa(pedinfo, AbstractString) && last(splitext(pedinfo))==".csv"
        pedfile2 = getabsfile(workdir, pedinfo)
        if !isfile(pedfile2)
            @error string(pedfile2, " does not exist")
        end
    end
    if !issubset(formatpriority,["GT","AD", "GP"])
        @error string("formatpriority=",formatpriority, ", not a subset of [GT, AD, GP]")
    end
    if isphysmap
        if !(recomrate>0)
            @error string("recomrate=", recomrate, " is not positive")
        end
    end
    msg = string("list of file options: \n",
            "genofile = ", genofile, "\n",
            "pedinfo = ", pedinfo, "\n",
            "formatpriority = ", formatpriority, "\n",
            "isphysmap = ", isphysmap, "\n",
            "recomrate = ", recomrate, "\n",
            "commentstring = ", commentstring, "\n",
            "workdir = ", workdir)
    MagicBase.printconsole(io,verbose,msg)
end

function check_mask_arg(foundermask::Real, offspringmask::Real, minread::Integer;
    outext::AbstractString=".vcf.gz")
    if !(0<=foundermask<=1)
        @error string("foundermask=", foundermask, " is not in [0,1]")
    end
    if !(0<offspringmask<=1)
        @error string("offspringmask=", offspringmask, " is not in (0,1]")
    end
    if !(minread>0)
        @error string("minread=", minread, " is not positive")
    end
    if !in(outext, [".vcf.gz",".csv.gz", ".vcf",".csv"])
        @error string("outext=",outext, ", not in [.vcf.gz, .csv.gz, .vcf, .csv]")
    end
    nothing
end

function mask_founder_offspring!(magicgeno::MagicGeno,undermask::MagicGeno;
    skipmarker::Real,
    ismaskfounder::Bool,minread::Integer,
    foundermask::Real, offspringmask::Real,
    io::Union{Nothing,IO}=nothing,
    verbose::Bool=true)
    starttime = time()
    if ismaskfounder
        foundermask <= 0.0 && return (0,missing)
    else
        offspringmask <= 0.0 && return (0,missing)
    end
    nmask = 0
    nnonmiss = 0
    nchr = length(magicgeno.markermap)
    for chr in 1:nchr           
        if ismaskfounder            
            ismask = rand(Bernoulli(foundermask),size(magicgeno.foundergeno[chr])...)
        else            
            ismask = rand(Bernoulli(offspringmask),size(magicgeno.offspringgeno[chr])...)
        end                
        ismiss_chr,misscode_chr = get_ismissls(magicgeno,chr;minread, isfounder=ismaskfounder)     
        ismask .= map((x,y)-> x && !y, ismask, ismiss_chr)   
        b = mean(ismiss_chr,dims=2)[:,1] .>= skipmarker
        ismask[b,:] .= false
        missmtx = repeat(misscode_chr,1,size(ismask,2))
        isnonmask = .!ismask
        if ismaskfounder
            undermask.foundergeno[chr][isnonmask] .= missmtx[isnonmask]
            magicgeno.foundergeno[chr][ismask] .= missmtx[ismask]
        else
            undermask.offspringgeno[chr][isnonmask] .= missmtx[isnonmask]
            magicgeno.offspringgeno[chr][ismask] .= missmtx[ismask]
        end
        nmask += sum(ismask)
        nnonmiss += length(ismiss_chr) - sum(ismiss_chr)        
    end
    nsnp = sum(size.(magicgeno.markermap,1))
    if ismaskfounder
        nfounder = size(magicgeno.magicped.founderinfo,1)
        miss_freq = round(1.0 - nnonmiss/(nfounder*nsnp),digits=4)
        msg = string("foundermask = ", foundermask)
    else
        noff = size(magicgeno.magicped.offspringinfo,1)
        miss_freq = round(1.0 - nnonmiss/(noff*nsnp),digits=4)       
        msg = string("offspringmask = ", offspringmask)
    end
    tused = round(time()-starttime, digits=1)
    msg *= string(", #mask = ", nmask,
    ", #nonmiss = ", nnonmiss,
    ", miss_freq = ", miss_freq,
    ", tused=", tused,"s")       
    MagicBase.printconsole(io,verbose,msg)
    nmask,nnonmiss
end

function get_ismissls(magicgeno::MagicGeno, chr::Integer;    
    minread::Integer, isfounder::Bool)    
    if isfounder
        chrgeno = magicgeno.foundergeno[chr]
        chrformat = magicgeno.markermap[chr][!,:founderformat]
    else
        chrgeno = magicgeno.offspringgeno[chr]
        chrformat = magicgeno.markermap[chr][!,:offspringformat]
    end
    ismiss_chr = falses(size(chrgeno)...)
    misscode_chr = Vector(undef,length(chrformat))
    formatset = unique(chrformat)
    for format in formatset
        ii = chrformat .== format
        subgeno = view(chrgeno,ii,:)
        if format == "GT_unphased"
            miss_set = ["NN","1N","2N","N1","N2"]                
            ismiss_chr[ii,:] .= [in(i,miss_set) for i in subgeno]
            misscode_chr[ii] .= "NN"
        elseif format == "AD"                
            ismiss_chr[ii,:] .= [sum(i) < minread for i in subgeno]
            for i in findall(ii)
                misscode_chr[i] = [0,0]
            end
        else
            misscode = MagicBase.get_missingcode(format)
            for i in findall(ii)
                misscode_chr[i] = misscode
            end
            ismiss_chr[ii,:] .= [i == misscode for i in subgeno]
        end
    end        
    ismiss_chr,misscode_chr
end
