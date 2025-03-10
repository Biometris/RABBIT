

"""
    magicreconstruct(genofile, pedinfo;
        formatpriority, isphysmap, recomrate, commentstring,kwargs...)

haplotye reconstruction from genofile and pedinfo.

# Positional arguments

`genofile::AbstractString` genotypic data file.

`pedinfo::Union{MagicBase.JuncDist,AbstractString}` specifies pedigree information via a pedigree fille or a string designcode or via a struct juncdist::JuncDist.

# Keyword arguments

See [`formmagicgeno`](@ref) for the arguments (`formatpriority`, `isphysmap`,
`recomrate`,`commentstring`) that are used for formming magicgeno, 
except that formatpriority=["GP", "AD", "GT"] by default. 

See [`magicreconstruct!`](@ref) for the other arguments.

# Examples
```julia-repl
julia> magicreconstruct("geno.vcf.gz","4ril_self3")
```
"""
function magicreconstruct(genofile::AbstractString,
    pedinfo::Union{JuncDist,AbstractString};
    formatpriority::AbstractVector=["GP", "AD", "GT"],
    isphysmap::Bool=false,
    recomrate::Real=1.0,
    model::AbstractString="jointmodel",    
    likeparameters::LikeParameters=LikeParameters(),   
    israndallele::Bool=true, 
    isfounderinbred::Bool=true,        
    chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    hmmalg::AbstractString="forwardbackward",
    isignorephase::Bool=false, 
    isMMA::Bool=false,
    nplot_subpop::Integer = 10,
    posteriordigits::Integer = 4, 
    thincm::Real = 0, 
    isparallel::Bool=true,
    workdir::AbstractString=pwd(),
    tempdirectory::AbstractString = tempdir(),
    commentstring::AbstractString="##",
    outext::AbstractString=".csv.gz",
    outstem::Union{Nothing,AbstractString}="outstem",
    logfile::Union{Nothing,AbstractString,IO}= (isnothing(outstem) ? nothing : string(outstem,"_magicreconstruct.log")),
    verbose::Bool=true)
    starttime = time()
    io = MagicBase.set_logfile_begin(logfile, workdir, "magicreconstruct"; verbose,delim="=")    
    if isa(logfile, AbstractString)
        msg=string("Julia version: ",VERSION)
        MagicBase.printconsole(io,verbose,msg)
        MagicBase.printpkgst(io,verbose,"MagicReconstruct")
    end
    MagicBase.check_infiles(genofile,pedinfo; isbreedped=false, io, commentstring,workdir,verbose)
    isnothing(outstem) || (outstem = outstem*"_magicreconstruct")
    info_file_arg(genofile, pedinfo, formatpriority,isphysmap, recomrate,
        commentstring,workdir, io, verbose)
    (;foundererror, offspringerror, seqerror) = likeparameters
    check_common_arg(model, foundererror, offspringerror, seqerror, workdir, tempdirectory;
        chrsubset,snpsubset)
    check_reconstruct_arg(hmmalg,outext)
    tused = @elapsed begin 
        magicgeno=formmagicgeno(genofile,pedinfo;
            isfounderinbred, formatpriority, isphysmap, recomrate, commentstring,workdir)
        mem1 = round(Int, memoryuse()/10^6)
        GC.gc()
        mem2 = round(Int, memoryuse()/10^6)
    end
    msg = string("formmagicgeno, tused=", round(tused,digits=1), 
        "seconds, mem=",mem1,"|",mem2,"MB")
    MagicBase.printconsole(io,verbose,msg)                
    magicancestry=magicreconstruct!(magicgeno;
        model, likeparameters, israndallele, 
        isfounderinbred,  
        chrsubset, snpsubset,
        isparallel,hmmalg, posteriordigits, isignorephase, isMMA, tempdirectory,
        nplot_subpop, thincm, workdir,outstem,outext,logfile=io, verbose)
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, io, tused,"magicreconstruct"; verbose,delim="=")
    magicancestry
end


"""
    magicreconstruct!(magicgeno::MagicGeno; kwargs...)

haplotype reconstruction from magicgeno.

# Keyword arguments

`model::AbstractString="jointmodel"`:  prior depedence of ancestral prior process
  along the two homologous chromosomes within an offspring. It must be "depmodel",
  "indepmodel", or "jointmodel". 

`israndallele::Bool=true`: if true, genotyping error model follows the random allelic model, and otherwise the random genotypic model. 

`isfounderinbred::Bool=true`: if true, founders are inbred, and otherwise they are outbred.

`chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing`: subset of chromosome indices.
  `nothing` denotes all chromosomes. Delete chromosome indices that are out of range.

`snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing`: subset of marker indices within each chromosome.
  `nothing` denotes all markers. Marker indices that are larger than the number of markers
  within the chromosome are deleted.

`hmmalg::AbstractString="forwardbackward"`: HMM alogrithm for haplotype
  reconstruction, and it must be either "forwardbackward" or "viterbi".

`isignorephase::Bool=false`: if true, the phases of offspring genotypes are ignored. 

`isMMA::Bool=true`: if true, the Mathematica version of RABBIT.

`nplot_subpop::Integer=10`: plots for up to nplot_subpop offspring in each subpopulation. 

`thincm::Real = 0`: thin ancestry results so that inter-marker distances > thincm. 

`isparallel::Bool=true`: if true, multicore computing over chromosomes.

`workdir::AbstractString=pwd()`: working directory for input and output files.

`tempdirectory::AbstractString = tempdir()`: temporary directory for inter-mediate results.

`outstem::Union{Nothing,AbstractString}="outstem"`: stem of output filenames.

`outext::AbstractString=".csv.gz"`: extension of output file for imputed geno.

`logfile::Union{Nothing,AbstractString,IO}= (isnothing(outstem) ? nothing : string(outstem,"_magicimpute.log"))`:
  log file or IO for writing log. If it is nothing, no log file.

`verbose::Bool=true`: if true, print details on the stdout.

# Examples
```julia-repl
julia> magicgeno = formmagicgeno("geno.vcf.gz","ped.csv")
julia> magicreconstruct(magicgeno,model="jointmodel")
```
"""
function magicreconstruct!(magicgeno::MagicGeno;
    model::AbstractString="jointmodel",
    likeparameters::LikeParameters=LikeParameters(),   
    israndallele::Bool=true, 
    isfounderinbred::Bool=true,    
    chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    hmmalg::AbstractString="forwardbackward",
    isignorephase::Bool=false, 
    isMMA::Bool=false,
    nplot_subpop::Integer=10,
    posteriordigits::Integer = 4, 
    thincm::Real = 0, 
    isparallel::Bool=true,
    workdir::AbstractString=pwd(),
    tempdirectory::AbstractString = tempdir(),
    outstem::Union{Nothing,AbstractString}="outstem",
    outext::AbstractString=".csv.gz",
    logfile::Union{Nothing,AbstractString,IO}= (isnothing(outstem) ? nothing : string(outstem,".log")),
    verbose::Bool=true)
    starttime = time()
    io = MagicBase.set_logfile_begin(logfile, workdir, "magicreconstruct"; verbose)        
    isa(logfile, AbstractString) && printpkgst(io,verbose,"MagicReconstruct")
    # check and print args
    (;foundererror,offspringerror,seqerror) = likeparameters
    check_common_arg(model, foundererror,offspringerror,seqerror,
        workdir, tempdirectory; chrsubset,snpsubset)
    check_reconstruct_arg(hmmalg,outext)        
    model = MagicBase.reset_model(magicgeno.magicped,model;io,verbose)    
    if model == "depmodel" && !isignorephase
        msg = "reset isignorephase = true for depmodel"
        MagicBase.printconsole(io,verbose,msg)
    end
    isparallel = isparallel && nprocs() > 1
    hmmalg = lowercase(hmmalg)
    msg = string("list of options: \n",
        "model = ", model, "\n",        
        "likeparameters = ", likeparameters,"\n",                
        "chrsubset = ", isnothing(chrsubset) ? "all chromosomes" : chrsubset,"\n",
        "snpsubset = ", isnothing(snpsubset) ? "all markers" : snpsubset,"\n",
        "hmmalg = ", hmmalg,"\n",
        "isignorephase = ", isignorephase,"\n",        
        "nplot_subpop = ", nplot_subpop, "\n",
        "thincm = ", thincm, "\n",
        "posteriordigits = ", posteriordigits,"\n",
        "isparallel = ", isparallel, isparallel ? string("(nworker=",nworkers(),")") : "", "\n",
        "workdir = ",workdir,"\n",
        "tempdirectory = ",tempdirectory,"\n",
        "outstem = ", isnothing(outstem) ? "no output files" : outstem,"\n",
        "outext = ",outext,"\n",
        "logfile = ",logfile,"\n",
        "verbose = ",verbose)
    MagicBase.printconsole(io,verbose,msg)
    MagicBase.reset_juncdist!(magicgeno.magicped,model;io,verbose,isfounderinbred)    
    submagicgeno!(magicgeno; chrsubset,snpsubset)        
    isignorephase && reset_ignorephase!(magicgeno)
    MagicBase.check_required_geneticmap(magicgeno; io)
    MagicBase.check_markerorder(magicgeno; io,verbose)    
    isnothing(seqerror) && (seqerror = 0.001)
    MagicBase.rawgenoprob!(magicgeno; targets = ["founders"], seqerror, isfounderinbred)
    MagicBase.rawgenocall!(magicgeno; targets = ["founders"], callthreshold = 0.95, isfounderinbred)
    # MagicBase.rawgenocall!(magicgeno; targets = ["founders","offspring"], callthreshold = 0.95, isfounderinbred, isoffspringphased=true)
    founderformat, offspringformat = MagicBase.setunphasedgeno!(magicgeno)
    if in("GP",founderformat)
        @error string("founderformat=",founderformat, ", not GT_haplo or GT_phased")
    end
    if in("GT_unphased",founderformat)
        @error string("founderformat=",founderformat, ", founders are not phased. Use magicimpute to phase founders")
    end
    if in("GT_phased",founderformat)
        if isfounderinbred
            @error string("founderformat=",founderformat, ", inconsistent with isfounderinbred = ",isfounderinbred)        
        end
    end
    if in("GT_haplo",offspringformat)
        @error string("offspringformat=",offspringformat, ", offspring have haplotypes, not genotypes")
    end
    MagicBase.info_magicgeno(magicgeno;io,verbose)
    outtarfile =  reconstruct!(magicgeno;
        model, likeparameters, israndallele, isfounderinbred,
        hmmalg, isMMA, posteriordigits, isparallel,io, outext = ".csv.gz", 
        outstem, workdir, tempdirectory,verbose)        
    magicancestry = try 
        readmagicancestry(outtarfile)      
    catch err
        @warn string(err, ". Could not read ", outtarfile)
        nothing
    end
    if !isnothing(magicancestry)
        # detect offpsing outliers
        tused = @elapsed begin 
            posterior_prior_recom!(magicancestry; minprob=0.7)
            detect_outlier!(magicancestry; tukeyfence=3, istransform=false)       
            if !isnothing(outstem)
                recomfile = save_posterior_recom(magicancestry.magicped; outstem, workdir)                 
                plot_posterior_recom(recomfile; workdir, outstem)
            end
        end
        msg = string("detect outliers, tused=",round(tused,digits=1), "s")
        MagicBase.printconsole(io,verbose,msg)
        magicgeno.magicped.offspringinfo = deepcopy(magicancestry.magicped.offspringinfo)
        outlier0=magicancestry.magicped.offspringinfo[!,:isoutlier]
        outlier = skipmissing(outlier0)
        isempty(outlier) ? noutlier = 0 : noutlier = sum(outlier)
        noutlier == 0 && MagicBase.printconsole(io,verbose,string("no outlier offspring"))
        if noutlier > 0            
            offoutid = magicancestry.magicped.offspringinfo[outlier0 .=== true,:individual]
            msg = string(noutlier, " outlier offspring: ",join(offoutid,", "))
            MagicBase.printconsole(io,verbose,msg)
        end        
        MagicBase.printconsole(io,verbose,"saving large output might take quite a few minutes!")
    end
    if isnothing(outstem)                        
        rm(outtarfile; force=true)     
    else
        if isnothing(magicancestry)
            @warn string("return nothing and results saved in ", outtarfile)
        else
            tarmb = filesize(outtarfile)*1e-6 
            if tarmb < 100  # file size < 100MB       
                rm(outtarfile; force=true)
            end
            # save results        
            tused = @elapsed if outext == last(MagicBase.split_allext(outtarfile))
                outputfile = outtarfile
            else
                outputfile =string(outstem,"_ancestry",outext)        
                savemagicancestry(outputfile,magicancestry; workdir)                                
            end
            msg = string("save in ",outputfile, ", tused=", round(tused,digits=1),"s")
            MagicBase.printconsole(io,verbose,msg)              
            MagicBase.thinmagicancestry(outputfile; thincm,workdir, 
                verbose, outext,io, isdeloutlier = false, 
                outstem = outstem*"_ancestry"
            )            

            # TODO: too large pedidgree results in error in plotting
            if size(magicancestry.magicped.founderinfo,1) <= 50            
                try 
                    fn = string(outstem, "_ped.png")
                    gdesign=plotmagicped(magicancestry.magicped;
                        isfounderinbred,outfile = joinpath(workdir,fn))
                    msg = isa(gdesign, AbstractString) ? gdesign : string("save in ", fn)
                    MagicBase.printconsole(io,verbose,msg)
                catch
                    @warn string("Could not plot magicped")
                end
            end
            # for isperchrom in [true, false]
            #     save_posterior_recom(magicancestry;isperchrom, isfounderinbred,
            #         isautosome=true,workdir,outstem,io, verbose)                
            # end
            if model != "depmodel"
                try 
                    ginbredfile = string(outstem, "_inbredcoef.png")
                    ginbred =plotinbredcoef(magicancestry)
                    savefig(ginbred,getabsfile(workdir, ginbredfile))                
                    msg = string("save in ", ginbredfile)
                    MagicBase.printconsole(io,verbose,msg)
                catch err
                    @warn string(err, ". Could not plot inbredcoef")
                end
            end
            if hmmalg == "forwardbackward"
                figdir = getabsfile(workdir, string(outstem, "_probplots"))
                isdir(figdir) || mkdir(figdir)                  
                try 
                    probtypels = model == "depmodel" ? ["haploprob"] : ["haploprob", "genoprob"]
                    tused = @elapsed  for probtype in probtypels
                        saveprobplot(magicancestry; nplot_subpop, probtype, workdir=figdir, outstem)     
                    end                    
                catch err
                    @warn string(err, ". Could not plot condprob")
                finally 
                    figfiles = readdir(abspath(figdir); join=true)
                    b = length.(figfiles) .>= 256
                    if any(b)                        
                        msg = string("save condprob in figdir = ",figdir,", tused=",round(tused,digits=1), "s")                        
                        MagicBase.printconsole(io,verbose,msg)  
                        msg = string("NOT remove figdir that contains ", sum(b), " files with too long filename such as ", first(figfiles[b]))
                        MagicBase.printconsole(io,false,"Warning: "*msg)  
                        @warn msg
                    else
                        tarplotfile = string(outstem, "_probplots.tar")
                        MagicBase.create_tar(figdir, getabsfile(workdir,tarplotfile))                            
                        msg = string("save condprob in ",tarplotfile)
                        MagicBase.printconsole(io,verbose,msg)  
                        try 
                            rm(figdir;recursive=true)                        
                        catch
                            @warn string("Could not remove ",figdir)
                        end
                    end
                end                    
            else
                # hmmalg == "viterbi"      
                contfgl = MagicBase.tocontfgl(magicancestry)      
                try                   
                    subpop2off = MagicBase.get_subpop2offspring(magicancestry.magicped; isindex=true)    
                    offls = sort(reduce(vcat, [length(i)<nplot_subpop ? i : i[1:nplot_subpop] for i in values(subpop2off)]))                
                    savemosaic(contfgl; 
                        offspring = offls, 
                        outstem=outstem*"_viterbi",
                        workdir, io,verbose)      
                catch
                    @warn string("Could not plot mosaic")
                end
            end    
        end
    end
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, io, tused,"magicreconstruct"; verbose)
    magicancestry
end

function reset_ignorephase!(magicgeno::MagicGeno; isfounderphase=false)
    formatcol = isfounderphase ? "founderformat" : "offspringformat"
    for chr in 1:length(magicgeno.markermap)
        formatls = magicgeno.markermap[chr][!,formatcol]
        b = formatls .== "GT_phased"
        if any(b)
            if isfounderphase
                subgeno = view(magicgeno.foundergeno[chr], b,:)
            else
                subgeno = view(magicgeno.offspringgeno[chr], b,:)
            end
            for i in eachindex(subgeno)                
                subgeno[i] = join(subgeno[i])                                
            end
            magicgeno.markermap[chr][b,formatcol] .= "GT_unphased"
        end
        b = formatls .== "GP"
        if any(b)
            if isfounderphase
                subgeno = view(magicgeno.foundergeno[chr], b,:)
            else
                subgeno = view(magicgeno.offspringgeno[chr], b,:)
            end            
            for i in eachindex(subgeno)
                ls = subgeno[i]
                if length(ls) == 4                    
                    subgeno[i] = [ls[1],ls[2]+ls[3],ls[4]]
                end
            end
        end
    end
    magicgeno
end


function info_file_arg(genofile::AbstractString,
    pedinfo::Union{JuncDist,AbstractString},
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
    MagicBase.printconsole(io,verbose,string("nthreads=",Threads.nthreads()))    
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

function check_reconstruct_arg(hmmalg::AbstractString, outext::AbstractString)
    hmmalg = lowercase(hmmalg)
    if !in(hmmalg, ["forwardbackward","viterbi"])
        @error string("hmmalg=",hmmalg, ", not in [forwardbackward, viterbi]")
    end
    if !in(outext, [".tar.gz",".csv.gz", ".tar",".csv"])
        @error string("outext=",outext, ", not in [.tar.gz, .csv.gz, .tar, .csv]")
    end
end

function check_common_arg(model::Union{AbstractString,AbstractVector},     
    epsf::Union{Nothing,Real},
    epso::Union{Nothing,Real},
    seqerror::Union{Nothing,Real}, workdir::AbstractString,tempdirectory::AbstractString;
    chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing,
    snpsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing)    
    if isa(model,AbstractVector)
        for i in eachindex(model)
            model_i = lowercase(model[i])
            if !in(model_i, ["depmodel","indepmodel","jointmodel"])
                @error string("model[",i,"]=",model_i, " is not in [depmodel, indepmodel, jointmodel]")
            end
        end
    else
        model = lowercase(model)
        if !in(model, ["depmodel","indepmodel","jointmodel"])
            @error string("model=",model, " is not in [depmodel, indepmodel, jointmodel]")
        end
    end
    if !isnothing(epsf) && !(0<=epsf<1)
        @error string("epsf=", epsf,"; founder allelic error probability is not in [0,1)")
    end
    if !isnothing(epso) &&!(0<=epso<1)
        @error string("epso=", epso,"; offspring allelic error probability is not in [0,1)")
    end
    if !isnothing(seqerror) &&!(0<=seqerror<1)
        @error string("seqerror=", seqerror,"; sequencing error probability is not in [0,1)")
    end
    if !isdir(workdir)
        @error string("workdir=", workdir, " is not a directory")
    end
    if !isdir(tempdirectory)
        @error string("tempdirectory=", tempdirectory, " is not a directory")
    end
    if !isnothing(chrsubset)
        if !(eltype(chrsubset) <: Integer)
            @error string("eltype(chrsubset)=",eltype(chrsubset), "; elments of chrusbet are not integers")
        end
        chrmin = min(chrsubset...)
        if chrmin < 1
            @error string("min(chrsubset...)=",chrmin, " is not greater than 0")
        end
    end
    if !isnothing(snpsubset)
        if !(eltype(snpsubset) <: Integer)
            @error string("eltype(snpsubset)=",eltype(snpsubset), "; elments of chrusbet are not Integer")
        end
        snpmin = min(snpsubset...)
        if snpmin < 1
            @error string("min(snpsubset...)=",snpmin, " is not greater than 0")
        end
    end
    nothing
end
