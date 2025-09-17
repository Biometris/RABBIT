
"""
    magiccall(genofile, pedinfo; kwargs...)

single marker genotype call from genofile and pedinfo.

# Positional arguments

`genofile::AbstractString` genotypic data file.

`pedinfo::Union{MagicBase.JuncDist,AbstractString}` specifies pedigree information via a pedigree fille or a string designcode or via a struct juncdist::JuncDist.

# Keyword arguments

`model::Union{AbstractString,AbstractVector}="jointmodel"`:  prior depedence of ancestral prior process
  along the two homologous chromosomes within an offspring. It must be "depmodel",
  "indepmodel", or "jointmodel". 

`likeparam::Union{Nothing,LikeParam}=LikeParam(offspringerror=0.005)`: parameters for genotypic data model. 
  If isinfererror = true, parameters with values being nothing will be inferred.   

`softthreshlikeparam::SoftThreshLikeParam=SoftThreshLikeParam()`: markers with inferred likeparam values > softthreshlikeparam values will be deleted only if they are also outliers. 

`threshlikeparam::Union{Nothing,LikeParam}=ThreshLikeParam()`: markers with inferred likeparam values > threshlikeparam values will be deleted. 
 
`priorlikeparam::Union{Nothing, PriorLikeParam}=PriorLikeParam()`: priors for likelihood parameters. 
 
`isfounderinbred::Bool=true`: if true, founders are inbred, and otherwise they are outbred.

`israndallele::Bool=true`: if true, genotyping error model follows the random allelic model, and otherwise the random genotypic model. 

`israwcall::Bool= false`: if true, perform raw genotype calling. 

`threshcall::Real = 0.9`: genotypes are called if maximum posterior probability > threshcall.

`iscalloffspring::Bool=true`: if true, offspring genotypes are called.

`isdelmultiallelic::Bool=true`: if true, delete markers with >=3 alleles. 

`isdelmonomorphic::Bool=true`: if true, delete monomorphic markers. 

`minmaf::Real = 0.05`: delete makrs with minor allele frequency (MAF) < 0.05. 

`maxmiss::Real = 0.99`: delete makrs with genotype missing frequency > 0.99. 

`isinfererror::Bool = !israwcall`: if true, infer marker specific likelihood parameters that have values of nothing in likeparam. 

`samplesize::Integer = 200`: number of posterior samples of founder genotypes

`burnin::Integer = 20`: number of initial iterations to be discarded. 

`isparallel::Bool=true`: if true, parallel multicore computing over chromosomes.

`outstem::Union{Nothing,AbstractString}="outstem"`: stem of output filenames.

`outext::AbstractString=".vcf.gz"`: extension of output file for imputed geno.

`logfile::Union{Nothing, AbstractString,IO} = outstem*"_magiccall.log"`:
log file or IO for writing log. If it is nothing, no log file.

`workdir::AbstractString=pwd()`: working directory for input and output files.

`verbose::Bool=true`: if true, print details on the stdout.  

"""
function magiccall(genofile::AbstractString,pedinfo::Union{Integer,AbstractString};
    model::AbstractString="jointmodel", 
    isfounderinbred::Bool=true,           
    israndallele::Bool=true, 
    likeparam::Union{Nothing,LikeParam} = LikeParam(offspringerror=0.005), 
    softthreshlikeparam::SoftThreshLikeParam = SoftThreshLikeParam(),
    threshlikeparam::Union{Nothing, ThreshLikeParam} = ThreshLikeParam(),
    priorlikeparam::Union{Nothing, PriorLikeParam} = PriorLikeParam(), 
    tukeyfence::Real=3,    
    israwcall::Bool= false, 
    threshcall::Real = 0.9,
    byfounder::Integer=0,
    iscalloffspring::Bool=true, 
    isdelmultiallelic::Bool=true,
    isdelmonomorphic::Bool=true,    
    minmaf::Real = 0.05, # set monomorphic subpopulation to missing if maf < minmaf
    maxmiss::Real = 0.99,         
    isinfererror::Bool = !israwcall, 
    samplesize::Integer = 200,
    burnin::Integer=20,
    isparallel::Bool=true,    
    outstem::AbstractString = "outstem",
    outext::AbstractString = ".vcf.gz",
    logfile::Union{Nothing, AbstractString,IO} = outstem*"_magiccall.log",
    workdir::AbstractString=pwd(),
    verbose::Bool=true)
    starttime = time()
    formatpriority=["AD","GP","GT"]    
    commentstring="##"
    logio = MagicBase.set_logfile_begin(logfile, workdir, "magiccall"; verbose,delim="-")
    if isa(logfile, AbstractString)
        msg=string("Julia version: ",VERSION)
        MagicBase.printconsole(logio,verbose,msg)
        MagicBase.printpkgst(logio,verbose,"MagicCall")
    end
    MagicBase.check_infiles(genofile,pedinfo; isbreedped=false,
        io = logio, commentstring,workdir,verbose)
    inext = last(MagicBase.split_allext(genofile))
    inext in [".vcf",".vcf.gz"] || @error string("genofile ext must be .vcf or .vcf.gz")
    in_open = inext in [".vcf.gz"] ? GZip.open : open    
    isparallel = isparallel && nworkers()>1        
    msg = string("list of args: \n",
        "genofile = ", genofile, "\n",
        "pedinfo = ", pedinfo, "\n",
        "model = ", model, "\n",
        "likeparam = ", likeparam, "\n",
        "softthreshlikeparam = ", softthreshlikeparam, "\n",		
        "threshlikeparam = ", threshlikeparam, "\n",		
        "priorlikeparam = ", priorlikeparam, "\n",	        
        "tukeyfence = ", tukeyfence, "\n", 
        "isfounderinbred = ", isfounderinbred, "\n",        
        "israwcall = ", israwcall, "\n",        
        "israndallele = ", israndallele, "\n",        
        "byfounder = ", byfounder, "\n",        
        "iscalloffspring = ", iscalloffspring, "\n",                
        "isdelmultiallelic = ", isdelmultiallelic, "\n",
        "isdelmonomorphic = ", isdelmonomorphic, "\n",        
        "minmaf = ", minmaf, "\n",        
        "maxmiss = ", maxmiss, "\n",                
        "threshcall = ", threshcall, "\n",                
        "israwcall = ", israwcall, "\n",        
        "isinfererror = ", isinfererror, "\n",                
        "samplesize = ", samplesize, "\n",             
        "burnin = ", burnin, "\n",             
        "isparallel = ", isparallel, isparallel ? string("(nworker=",nworkers(),")") : "", "\n",
        "outstem = ", outstem, "\n",
        "outext = ", outext, "\n",
        "logfile = ", logfile, "\n",
        "workdir = ", workdir, "\n",
        "verbose = ", verbose)
    printconsole(logio,verbose,msg)    
    if isinfererror       
        if isnothing(threshlikeparam) 
            threshlikeparam = ThreshLikeParam()
            msg = string("reset threshlikeparam = ", threshlikeparam)
            printconsole(logio,verbose,msg)
        end    
        if isnothing(priorlikeparam) 
            priorlikeparam = PriorLikeParam()
            msg = string("reset priorlikeparam = ", priorlikeparam)
            printconsole(logio,verbose,msg)
        end    
        ischangedls, priorlikeparam = reset_priorlikeparam(priorlikeparam)
        if any(ischangedls)
            msg = MagicBase.get_info_likeparam(priorlikeparam; ismultiline=true)
            printconsole(logio,verbose,msg)        
        end    
        if isnothing(likeparam)
            likeparam = LikeParam()
            msg = string("reset likeparam = ", likeparam)
            printconsole(logio,verbose,msg)
        end        
        msg = MagicBase.get_info_likeparam(likeparam; isinfererror, ismultiline=true)
        printconsole(logio,verbose,msg)           
    else
        msg = string("baseerror=", MagicBase.get_likeproperty(likeparam, :baseerror))
        printconsole(logio,verbose,msg)    
    end
    genofile2 = MagicBase.getabsfile(workdir,genofile)
    isfile(genofile2) || error(string("Could not find ",genofile2))    
    lastcomment = MagicBase.findlastcomment(genofile2; commentstring)
    nheader = lastcomment + 1
    outext in [".vcf",".vcf.gz"] || @error string("outfile ext must be in [.vcf,.vcf.gz]")    
    # set parallel
    nsnp = MagicBase.vcf_count_markers(genofile2)     
    msg = string("#markers=", nsnp)
    if isparallel 
        nwork = nworkers()
        seglen = max(10,min(500,div(nsnp, nwork)))
        paragraphls = collect(Iterators.partition(1:nsnp,seglen*nwork))
        msg *= string(" divided into ", length(paragraphls), " blocks")        
        nwork > 1 && (msg *= string(", markers in each block being distributed over ", nwork, " workers"))        
    else
        paragraphls = nothing
    end    
    printconsole(logio,verbose,msg)       
    outstem *= "_magiccall"
    out_open = outext in [".vcf.gz",".csv.gz"] ? GZip.open : open
    calledfile = outstem*"_called"*outext
    calledfile2 = getabsfile(workdir,calledfile)            
    delmarkerfile = outstem*"_delmarker.vcf.gz"
    delmarkerfile2 = getabsfile(workdir,delmarkerfile)
    in_open(genofile2,"r") do inio                          
        out_open(calledfile2, "w") do outio                
            GZip.open(delmarkerfile2,"w") do delio
                magiccall_io(inio, outio, delio, logio, pedinfo,model, 
                    likeparam, threshlikeparam, priorlikeparam, 
                    formatpriority, israndallele,isfounderinbred, byfounder, 
                    isdelmultiallelic, isdelmonomorphic, minmaf,maxmiss, 
                    threshcall, israwcall, isinfererror,
                    iscalloffspring, samplesize, burnin, 
                    paragraphls,commentstring, nheader, workdir, verbose)            
            end
        end        
    end
    outfile = outstem*"_geno"*outext
    outfile2 = getabsfile(workdir,outfile)            
    if isinfererror
        tused = @elapsed outliers = shift_outliers(calledfile2, outfile2, delmarkerfile2; tukeyfence, softthreshlikeparam, threshlikeparam, workdir)
        msg = string("delete #outlier markers = ", length(outliers), ", tused=", round(tused,digits=1), "s")
        printconsole(logio,verbose,msg)        
    else
        mv(calledfile2, outfile2, force=true)
    end
    msg = string("save called_genofile: ", outfile)
    printconsole(logio,verbose,msg)        
    msg = string("save delmarker_genofile: ", delmarkerfile)
    printconsole(logio,verbose,msg)       
    if !israwcall
        for (i,mapfile) in enumerate([outfile2])                
            try 
                figerr = plotmarkererror(mapfile;tukeyfence, workdir)
                if !isnothing(figerr)                        
                    markererrfile = string(outstem, i==1 ? "" : "_delmarker", "_inferred_error.png")
                    try 
                        MagicBase.savefig(figerr,getabsfile(workdir,markererrfile))
                    catch err
                        msg = string("Could not savefig. ",err)
                        @warn msg
                        printconsole(logio, false, "Warning: "*msg)
                    end
                end 
            catch err                
                msg = string(err, ". Could not plot markererror for mapfile= ",mapfile)
                @warn msg
                printconsole(logio, false, "Warning: "*msg)
            end
        end
    end
    
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, logio, tused,"magiccall"; verbose,delim="-")
    return nothing
end

# function threshlike(;
#     foundererror=0.2, offspringerror=0.2, 
#     baseerror=0.05, allelicbias=0.9, 
#     allelicoverdispersion=1.0, allelicdropout=0.05)
#     ThreshLikeParam(foundererror=foundererror, offspringerror=offspringerror, 
#         baseerror=baseerror, allelicbias=allelicbias, 
#         allelicoverdispersion=allelicoverdispersion, allelicdropout=allelicdropout)
# end

function reset_priorlikeparam(priorlikeparam::PriorLikeParam)
    priorls = []
    errnamels = propertynames(priorlikeparam)
    ischangedls = falses(length(errnamels))
    for i in eachindex(ischangedls)
        errname = errnamels[i]
        prior = getproperty(priorlikeparam,errname)
        if isnothing(prior)
            ischangedls[i] = true
            if errname in [:foundererror, :offspringerror]
                prior = Beta(1,19)
            elseif errname in [:baseerror]
                prior = Beta(1,199)
            elseif errname in [:allelicbias]
                prior = Beta(1.01,1.01)
            elseif errname in [:allelicoverdispersion]
                prior = Exponential(0.3)            
            elseif errname in [:allelicdropout]
                prior = Beta(1,19)
            end
        end
        push!(priorls,prior)
    end
    ischangedls, PriorLikeParam(foundererror=priorls[1],offspringerror=priorls[2], peroffspringerror=priorls[3],
        baseerror =priorls[4], allelicbias=priorls[5],allelicoverdispersion=priorls[6], allelicdropout=priorls[7])
end

function magiccall_io(inio::IO,outio::IO,delio::IO,
    logio::Union{Nothing, IO},
    pedinfo::Union{Integer,AbstractString}, 
    model::AbstractString, 
    likeparam::LikeParam,        
    threshlikeparam::ThreshLikeParam,        
    priorlikeparam::PriorLikeParam,        
    formatpriority::AbstractVector,
    israndallele::Bool,
    isfounderinbred::Bool,
    byfounder::Integer,     
    isdelmultiallelic::Bool,
    isdelmonomorphic::Bool,
    minmaf::Real,
    maxmiss::Real,
    threshcall::Real, 
    israwcall::Bool,
    isinfererror::Bool,
    iscalloffspring::Bool, 
    samplesize::Integer, 
    burnin::Integer,
    paragraphls::Union{Nothing,AbstractVector},    
    commentstring::AbstractString,
    nheader::Integer,     
    workdir::AbstractString, verbose::Bool)        
    for _ in 1:nheader-1
        line = readline(inio,keep=true)        
        write(outio, line)
        write(delio, line)
    end
    # parse title row
    titlerow = split(readline(inio,keep=false),"\t")
    magicped, fcols, offcols,newtitlerow = MagicBase.parse_titlerow(titlerow, pedinfo;
        keepvcf=true, commentstring,workdir,logio)        
    write(outio,newtitlerow,"\n")
    write(delio,newtitlerow,"\n")
    isautosome = true    
    prioraa = MagicReconstruct.calprior(magicped,model; isfounderinbred,isautosome)    
    popmakeup = MagicReconstruct.calpopmakeup(magicped, model,prioraa; isfounderinbred,isautosome)
    nstate, nfgl = MagicReconstruct.hmm_nstate_nfgl(popmakeup)    
    # parse rest rows
    msg = string("begin, #founders=",length(fcols), ", #offspring=",length(offcols))
    printconsole(logio,verbose,msg)        
    missingset = [".","./.", ".|."]    
    nmarker = 0
    nmultia = 0
    nmono = 0
    nmaxmiss = 0
    nlargeerr = zeros(Int,7)
    startt0 = time()            
    if isnothing(paragraphls)        
        startt = time()
        while !eof(inio)
            if rem(nmarker, 100) == 0 
                startt = time()
                GC.gc()                
            end
            nmarker += 1                        
            rowstring = readline(inio,keep=false)
            res = magiccall_rowgeno(rowstring,fcols,offcols; israwcall, isdelmultiallelic, isdelmonomorphic,
                israndallele, isfounderinbred,byfounder,model,popmakeup, formatpriority,
                minmaf,maxmiss,threshcall, isinfererror, iscalloffspring, samplesize, burnin,
                likeparam, threshlikeparam, priorlikeparam, missingset, nstate,nfgl)        
            if res[1] == "maxmiss"
                write(delio,res[2],"\n")
                nmaxmiss += 1
            elseif res[1] == "biallelic"
                write(outio,res[2],"\n")
            elseif res[1] == "multiallelic"
                nmultia += 1
                if isdelmultiallelic  
                    write(delio,res[2],"\n")
                else
                    write(outio,res[2],"\n")
                end
            elseif res[1] == "monomorphic"
                nmono += 1
                if isdelmonomorphic  
                    write(delio,res[2],"\n")
                else
                    write(outio,res[2],"\n")
                end
            elseif occursin(r"^largeerror",res[1])
                up_nlargeerror!(nlargeerr,res[1])                
                write(delio,res[2],"\n")            
            else
                @error string("unknown resid=",res[1])                
            end            
            if rem(nmarker, 100) == 0
                nincl = nmarker - nmultia - nmono - sum(nlargeerr) - nmaxmiss
                msg = string("#marker=", nmarker, ", #marker_incl=",nincl)          
                nmultia > 0 && (msg *= string(", #multiallelic=", nmultia))
                nmono > 0 && (msg *= string(", #monomorphic=", nmono))
                sum(nlargeerr) > 0 && (msg *= string(", #largeerror=", nlargeerr))
                nmaxmiss > 0 && (msg *= string(", #>=maxmiss=", nmaxmiss))
                msg *= string(", tused=", round(time()-startt,digits=1),"s")
                printconsole(logio,verbose,msg)
            end
        end  
    else    
        startt = time()            
        for paragraph in paragraphls          
            startt = time()              
            nline = length(paragraph)            
            nmarker += nline            
            multirows = [readline(inio,keep=false) for i in 1:nline]                        
            res = pmap(x-> magiccall_rowgeno(x,fcols,offcols; israwcall, isdelmultiallelic, isdelmonomorphic,
                israndallele,isfounderinbred,byfounder, model, popmakeup, formatpriority, 
                minmaf, maxmiss, threshcall, isinfererror,iscalloffspring, samplesize, burnin,
                likeparam, threshlikeparam, priorlikeparam, missingset, nstate,nfgl), multirows)            
            for i in res
                if i[1] == "maxmiss"
                    write(delio,i[2],"\n")
                    nmaxmiss += 1
                elseif i[1] == "biallelic"
                    write(outio,i[2],"\n")
                elseif i[1] == "multiallelic"
                    nmultia += 1
                    if isdelmultiallelic  
                        write(delio,i[2],"\n")
                    else
                        write(outio,i[2],"\n")
                    end                
                elseif i[1] == "monomorphic"
                    nmono += 1
                    if isdelmonomorphic  
                        write(delio,i[2],"\n")
                    else
                        write(outio,i[2],"\n")
                    end
                elseif occursin(r"largeerror_",i[1])
                    up_nlargeerror!(nlargeerr,i[1])                    
                    write(delio,i[2],"\n")                     
                else
                    @error string("unknown resid=",i[1])                
                end
            end    
            GC.gc()
            @everywhere GC.gc()
            nincl = nmarker - nmultia - nmono - sum(nlargeerr) - nmaxmiss
            msg = string("#marker=", nmarker, ", #marker_incl=",nincl)             
            nmultia > 0 && (msg *= string(", #multiallelic=", nmultia))            
            nmono > 0 && (msg *= string(", #monomorphic=", nmono))
            sum(nlargeerr) > 0 && (msg *= string(", #largeerror=", nlargeerr))
            nmaxmiss > 0 && (msg *= string(", #>=maxmiss=", nmaxmiss))
            msg *= string(", tused=", round(time()-startt,digits=1),"s")
            printconsole(logio,verbose,msg)            
        end
    end
    msg = string("end, #founders=",length(fcols), ", #offspring=",length(offcols))
    nincl = nmarker - nmultia - nmono - sum(nlargeerr) - nmaxmiss
    msg *= string(", #marker=", nmarker, ", #marker_incl=",nincl)          
    nmultia > 0 && (msg *= string(", #multiallelic=",nmultia))    
    nmono > 0 && (msg *= string(", #monomorphic=",nmono))
    sum(nlargeerr) > 0 && (msg *= string(", #largeerror=",nlargeerr))
    nmaxmiss > 0 && (msg *= string(", #>=maxmiss=", nmaxmiss))
    msg *= string(", tused=", round(time()-startt0,digits=1),"s")
    printconsole(logio,verbose,msg)
end

function up_nlargeerror!(nlargeerr::AbstractVector,largeerrid::AbstractString)
    if largeerrid == "largeerror_foundererror"
        nlargeerr[1] += 1
    elseif largeerrid == "largeerror_offspringerror"
        nlargeerr[2] += 1
    elseif largeerrid == "largeerror_baseerror"
        nlargeerr[4] += 1
    elseif largeerrid == "largeerror_allelicbias"
        nlargeerr[5] += 1
    elseif largeerrid == "largeerror_allelicoverdispersion"
        nlargeerr[6] += 1
    elseif largeerrid == "largeerror_allelicdropout"
        nlargeerr[7] += 1
    end
    nlargeerr
end

function magiccall_rowgeno(rowstring::AbstractString,
    fcols::AbstractVector,offcols::AbstractVector; 
    israwcall::Bool,
    isdelmultiallelic::Bool,
    isdelmonomorphic::Bool,
    israndallele::Bool,
    isfounderinbred::Bool,
    byfounder::Integer, 
    model::AbstractString, 
    popmakeup::AbstractDict,    
    formatpriority::AbstractVector,
    minmaf::Real,
    maxmiss::Real,
    threshcall::Real,     
    isinfererror::Bool,
    iscalloffspring::Bool, 
    samplesize::Integer, 
    burnin::Integer,
    likeparam::LikeParam,    
    threshlikeparam::ThreshLikeParam,    
    priorlikeparam::PriorLikeParam,    
    missingset::AbstractVector,    
    nstate::Integer, 
    nfgl::Integer)
    rowgeno = split(rowstring,"\t")
    ismultiallele = length(split(rowgeno[5],",")) > 1 # col5=alternative                     
    ismultiallele && isdelmultiallelic && return ("multiallelic", rowstring)
    baseerror = MagicBase.get_likeproperty(likeparam, :baseerror)
    res = extract_rowgeno!(rowgeno, fcols,offcols,formatpriority,
        missingset,isfounderinbred,isdelmonomorphic,baseerror)    
    res == "monomorphic" && return ("monomorphic",join(rowgeno,"\t"))
    inputformat, fgeno, founderformat, offgeno, offspringformat = res
    isoffmiss = rowgeno[offcols] .== "NA"
    # @info "test" fgeno founderformat offgeno offspringformat maxlog = 10
    # imputation/correction and error estimations    
    if israwcall                        
        fhaplo_GT, fhaplo_GP = infer_fhaplo_rawcall(fgeno,founderformat; baseerror,isfounderinbred,threshcall)                 
        esterrors = infer_error_rawcall(offgeno,offspringformat; popmakeup, model,isinfererror,
            likeparam,priorlikeparam,israndallele)                        
    else        
        fhaplo_GT, fhaplo_GP, esterrors = infer_fhaploerror_singlesite(fgeno,founderformat, offgeno,offspringformat; 
            model, popmakeup, isfounderinbred, israndallele, isinfererror, byfounder, 
            threshcall, samplesize, burnin,likeparam, priorlikeparam)                    
    end                    
    add_error_info!(rowgeno,esterrors)    
    # update rowgeno for founder geno         
    outformat =  copy(inputformat)    
    in(:GP, outformat) || push!(outformat,:GP)
    for i in fcols
        rowgeno[i] == "NA" && (rowgeno[i] ="." )
    end     
    if in(:GT, outformat)
        i_gt = findfirst(isequal(:GT),outformat)
        for i in eachindex(fcols)
            gg = split(rowgeno[fcols[i]],":")
            gg[i_gt] = fhaplo_GT[i]
            rowgeno[fcols[i]] = join(gg,":")
        end
    else
        pushfirst!(outformat,:GT)
        for i in eachindex(fcols)
            rowgeno[fcols[i]] = string(fhaplo_GT[i], ":", rowgeno[fcols[i]])
        end
        # To add GT for offcols
        for i in eachindex(offcols)
            rowgeno[offcols[i]] = string("./.:", rowgeno[offcols[i]])
        end        
        if offspringformat == "GT" 
            @error string("inconsistent among offspringformat=",offspringformat,", FORMAT=",rowgeno[9]) 
        end
    end    
    if in(:GP, outformat)
        i_gp = findfirst(isequal(:GP),inputformat)
        if isnothing(i_gp)
            outformat[end] == :GP || @error string("unexpected outformat=",outformat)
            for i in eachindex(fcols)                
                rowgeno[fcols[i]] = string(rowgeno[fcols[i]], ":", fhaplo_GP[i])
            end           
            # To add GP for offcols
            for i in eachindex(offcols)
                rowgeno[offcols[i]] = string(rowgeno[offcols[i]],":.")
            end      
        else
            for i in eachindex(fcols)
                gg = split(rowgeno[fcols[i]],":")            
                gg[i_gp] = fhaplo_GP[i] 
                rowgeno[fcols[i]] = join(gg,":")
            end
        end
    else
        @error string("unexpected outformat=",outformat, " for rowgeno=",rowgeno)        
    end
    ismono = in(unique(reduce(vcat, split.(fhaplo_GT,"/"))),[["0"],["1"]])
    if isdelmonomorphic && ismono                
        return ("monomorphic",join(rowgeno,"\t"))
    end    
    islargeerror, largeerrorid = get_isdelmarker(esterrors,threshlikeparam)
    if islargeerror                 
        return ("largeerror_"*largeerrorid,join(rowgeno,"\t"))
    end        
    resid = ismultiallele ? "multiallelic" : (ismono ? "monomorphic" : "biallelic")

    # update rowgeno for offspring geno
    if iscalloffspring
        if israwcall
            offpostprob = infer_offpostprob_rawcall(offgeno,offspringformat, esterrors,israndallele)
        else
            offpostprob = infer_offpostprob_rawcall(offgeno,offspringformat, esterrors,israndallele)
            # offpostprob = infer_offpostprob_singlesite(fhaplo_GT,offgeno, offspringformat, popmakeup,esterrors,nstate,nfgl, 
            #     israndallele,isfounderinbred)
        end
        offgeno_GT, offgeno_GP = postprob2vcfgeno(offpostprob; callthreshold=threshcall, digits=4)        
        # filter by missing and minmaf
        noff = length(offgeno_GT)    
        alleles = reduce(vcat, split.(offgeno_GT,"/"))
        deleteat!(alleles, alleles .== ".")
        n1 = sum(alleles .== "0")
        n2 = sum(alleles .== "1")
        missfreq = 1 - (n1+n2)/(2*noff)
        if missfreq > maxmiss
            resid = "maxmiss"
        else            
            if minmaf > 0 && (n1+n2) > 1/minmaf            
                if min(n1,n2)/(n1+n2) < minmaf
                    resid = "monomorphic"
                end
            end
        end    
        # save 
        if in(:GT, outformat)
            i_gt = findfirst(isequal(:GT),outformat)
            for i in eachindex(offcols)
                gg = split(rowgeno[offcols[i]],":")            
                gg[i_gt] = offgeno_GT[i]
                rowgeno[offcols[i]] = join(gg,":")
            end
        else
            @error string("unexpected outformat=",outformat, " for rowgeno=",rowgeno)
        end
        if in(:GP, outformat)
            i_gp = findfirst(isequal(:GP),outformat)            
            for i in eachindex(offcols)
                gg = split(rowgeno[offcols[i]],":")
                gg[i_gp] = offgeno_GP[i]
                rowgeno[offcols[i]] = join(gg,":")
            end        
        else
            @error string("unexpected outformat=",outformat, " for rowgeno=",rowgeno)                               
        end    
    end
    rowgeno[offcols[isoffmiss]] .= "."
    rowgeno[9] = join(outformat,":")
    (resid, join(rowgeno,"\t"))
end

function infer_fhaplo_rawcall(fgeno::AbstractVector, founderformat::String; 
    baseerror::Real=0.001, isfounderinbred::Bool,threshcall::Real)  
    if founderformat == "GT"
        alleledict = Dict("1"=>"0","2"=>"1","N"=>".")
        fhaplo_GT = [ismissing(i) ? "./." : join(replace(split(i,""),alleledict...),"/") for i in fgeno]
        fhaplo_GP = ["." for _ in eachindex(fgeno)]
    elseif founderformat == "GT_haplo"
        haplodict = Dict("1"=>"0|0","2"=>"1|1","N"=>".|.",missing=>".|.")        
        fhaplo_GT = [haplodict[i] for i in fgeno]
        fhaplo_GP = ["." for _ in eachindex(fgeno)]
    elseif founderformat == "GP"
        fhaplo_GT = ["./." for _ in eachindex(fgeno)]
        fhaplo_GP = repace.(fgeno, "N"=>".")
    elseif founderformat =="GP_haplo"
        fhaplo_GP = ["." for _ in eachindex(fgeno)]
    elseif founderformat in ["AD","GP"]
        if founderformat == "AD"            
            if isfounderinbred 
                fhaplo_postprob = MagicBase.genoprobhaplo.(fgeno,baseerror)
            else
                fhaplo_postprob = MagicBase.genoprobdiplo_biallelic.(fgeno, baseerror)
            end
        else
            # GP
            fhaplo_postprob = isfounderinbred ? MagicBase.genoprob2haploprob.(fgeno) : fgeno    
        end                 
        fhaplo_GT, fhaplo_GP = postprob2vcfgeno(fhaplo_postprob; callthreshold=threshcall, digits=2)
    else
        @error string("TODO rawcall with founderformat =",founderformat)            
    end
    fhaplo_GT, fhaplo_GP
end

function postprob2vcfgeno(postproblsls::AbstractVector; callthreshold=0.95,digits=6)        
    res_GT = String[]
    res_GP = String[]
    haplodict = Dict("1"=>"0/0","2"=>"1/1","N"=>"./.",missing=>"./.")        
    alleledict = Dict("1"=>"0","2"=>"1","N"=>".")
    for gp in postproblsls
        if ismissing(gp)
            push!(res_GT, "./.")
            push!(res_GP, ".")        
        else
            ishaplo = length(gp) == 2
            gp2 = round.(gp; digits)
            if ishaplo
                g = MagicBase.callfromprob(gp,callthreshold; ishaplo)
                g2 = haplodict[g]
                push!(res_GT, g2)        
                push!(res_GP, join([gp2[1],0.0, gp2[2]],","))        
            else
                g = MagicBase.callfromprob(gp,callthreshold; ishaplo)
                g2 = join([alleledict[i] for i in split(g,"")],"/")
                push!(res_GT, g2)
                push!(res_GP, join(gp2,","))        
            end
        end
    end
    res_GT, res_GP
end




function extract_rowgeno!(rowgeno::AbstractVector,fcols::AbstractVector,offcols::AbstractVector,
    formatpriority::AbstractVector,
    missingset::AbstractVector,    
    isfounderinbred::Bool,
    isdelmonomorphic::Bool,
    baseerror::Real)    
    ids = Symbol.(formatpriority)
    vals = 1:length(formatpriority)
    prioritytuple = (; zip(ids, vals)...)
    inputformat = Symbol.(split(rowgeno[9],":"))
    formatcode = [get(prioritytuple,i,-1) for i in inputformat] # -1 denotes missing fromat            
    # fgeno
    fgeno = Vector(undef, length(fcols))     
    founderformat = parse_rowgeno!(fgeno,rowgeno[fcols],formatcode, formatpriority, missingset)        
    parse2_internalgeno!(fgeno,founderformat; ishaplo=isfounderinbred)        
    if isfounderinbred        
        threshcall = 0.9
        if founderformat == "AD"            
            fgeno .= MagicBase.callfromprob.(MagicBase.genoprobhaplo.(fgeno, baseerror),threshcall; 
                isphased=false,ishaplo=true)            
        elseif founderformat == "GP"
            fgeno .= MagicBase.callfromprob.(fgeno, threshcall; isphased=false,ishaplo=isfounderinbred)
        elseif founderformat == "GT"
            0
        else 
            @error string("TODO for inbred founderformat =",founderformat)            
        end   
        founderformat = "GT_haplo"
    end    
    if occursin(r"^GT",founderformat)
        ismono = in(unique(split(join(fgeno),"")),[["1"],["2"]])
        if isdelmonomorphic && ismono        
            return "monomorphic"
        end  
    end
    # set offspring geno
    offgeno = Vector(undef, length(offcols))
    offspringformat = parse_rowgeno!(offgeno, rowgeno[offcols],formatcode, formatpriority, missingset)            
    parse2_internalgeno!(offgeno,offspringformat; ishaplo=false)   
    if occursin(r"^GT",offspringformat)
        alleleset = unique(split(join(offgeno),""))
        setdiff!(alleleset,["N"])
        ismono = in(alleleset,[["1"],["2"]])
        if isdelmonomorphic && ismono        
            return "monomorphic"
        end      
    end
    inputformat, fgeno, founderformat, offgeno,offspringformat
end

function get_isdelmarker(esterrors,threshlikeparam)    
    epsf,epso, _,baseerror, allelicbias,allelicoverdispersion,allelicdropout = esterrors
    epsf > threshlikeparam.foundererror && return (true, "foundererror")
    epso > threshlikeparam.offspringerror && return (true, "offspringerror")
    baseerror > threshlikeparam.baseerror && return (true, "baseerror")
    thresh = threshlikeparam.allelicbias
    (1-thresh <= allelicbias <= thresh) || return (true, "allelicbias")
    allelicoverdispersion > threshlikeparam.allelicoverdispersion && return (true, "allelicoverdispersion")
    allelicdropout > threshlikeparam.allelicdropout && return (true, "allelicdropout")
    (false, "")        
end

function add_error_info!(rowgeno::AbstractVector,esterrors)
    epsf,epso, _, baseerror, allelicbias,allelicoverdispersion,allelicdropout = esterrors
    info = string("FOUNDERERROR=",round(epsf,digits = 5))
    info *= string(";OFFSPRINGERROR=",round(epso,digits = 5))
    info *= string(";BASEERROR=",round(baseerror,digits = 5))
    info *= string(";ALLELICBIAS=",round(allelicbias,digits = 5))
    info *= string(";ALLELICOVERDISPERSION=",round(allelicoverdispersion,digits = 5))
    info *= string(";ALLELEDROPOUT=",round(allelicdropout,digits = 5))    
    inputinfo = strip(rowgeno[8])
    if !in(inputinfo, [".",""])
        info = string(inputinfo, ";", info)
    end
    rowgeno[8] = info    
    info
end

function infer_offpostprob_singlesite(fhaplo_GT::AbstractVector,offgeno::AbstractVector, offspringformat::AbstractString, 
    popmakeup::AbstractDict, esterrors::AbstractVector,     
    nstate::Integer,nfgl::Integer,israndallele::Bool,
    isfounderinbred::Bool)
    epsf,epso, _, baseerror, allelicbias,allelicoverdispersion,allelicdropout = esterrors
    if isfounderinbred
        fhaplo = first.(split.(fhaplo_GT,"/"))
    else
        fhaplo = reduce(vcat, split.(fhaplo_GT,"/"))
    end    
    replace!(fhaplo, "0"=>"1","1"=>"2","."=>"N")
    # @info "test" fhaplo fhaplo_GT maxlog=10
    res = Vector(undef, length(offgeno))
    for popid in keys(popmakeup)                
        offls = popmakeup[popid]["offspring"]
        nzstate =  popmakeup[popid]["nzstate"]
        nzorigin =  popmakeup[popid]["nzorigin"]
        ishaploid =  popmakeup[popid]["ishaploid"]		        
        initprob = spzeros(nstate)
        initprob[nzstate] .= popmakeup[popid]["initprob"]       
        ismiss = [(i==[0,0] || i==[0.25,0.5,0.25]) for i in view(offgeno, offls)]
        isnonmiss  = .!ismiss                        
        sitegeno = view(offgeno, offls[isnonmiss])              
        if ishaploid                                         
            if !isempty(sitegeno)
                fderivehaplo = calfderive(fhaplo, nzorigin,ishaploid)
                if offspringformat == "AD"
                    # genofrmat = AD             
                    like = MagicReconstruct.haplolikeGBS(sitegeno, epso,baseerror)          
                elseif offspringformat == "GT"
                    like = MagicReconstruct.haplolike_GT(sitegeno,epso)  
                else        
                    # genoformat = GP (genotpe prob), GL (log10 scaled likeliihood)
                    like = MagicReconstruct.haplolikeGBS(sitegeno, epso)            
                end   
                priorhaplo = MagicReconstruct.haploprior(epsf)                          
                ls = priorhaplo[:,fderivehaplo[nzstate]] * initprob[nzstate]
                post = [ls .* i for i in eachrow(like)]                
                res[offls[isnonmiss]] .= [normalize(i,1) for i in post]   
            end
            if any(ismiss)
                res[offls[ismiss]]  .= [[0.5,0.5] for _ in 1:sum(ismiss)]        
            end
        else
            nfgl == length(fhaplo) || @error string("inconsistent nfgl=",nfgl, ", nfgl2=",length(fhaplo), ", fhaplo=",fhaplo) maxlog=10
            nstate == nfgl^2 || error(string("mismatch between nfgl=",nfgl, " and nstate=",nstate))					       
            if !isempty(sitegeno)                        
                ibdbool = allequal.(nzorigin)
                nonibdbool = .!ibdbool
                initprob_ibd = initprob[ibdbool]
                initprob_nonibd = initprob[nonibdbool]
                fderivediplo = calfderive(fhaplo, nzorigin,ishaploid)
                if offspringformat == "AD"
                    # format = AD
                    like = MagicReconstruct.diplolikeGBS(sitegeno,epso,baseerror,
                        allelicbias,allelicoverdispersion,allelicdropout; israndallele)   
                elseif offspringformat == "GT"
                    like = MagicReconstruct.diplolike_GT(sitegeno,epso; israndallele)
                else        
                    # format = GP    
                    like = MagicReconstruct.diplolikeGBS(sitegeno,epso; israndallele)
                end                                
                priordiplo = MagicReconstruct.diploprior(epsf)                                       
                ls = priordiplo.nonibd[:,fderivediplo[nonibdbool]] * initprob_nonibd
                ls .+= priordiplo.ibd[:,fderivediplo[ibdbool]] * initprob_ibd 
                post = [ls .* i for i in eachrow(like)]
                res[offls[isnonmiss]] .= [begin 
                    v = normalize(i,1)
                    [v[1],v[2]+v[3],v[4]]
                end for i in post]
            end
            if any(ismiss)
                res[offls[ismiss]] .= [[0.25,0.5,0.25] for _ in 1:sum(ismiss)]
            end
        end                
    end        
    res
end

function calfderive(fhaplo::AbstractVector, nzorigin::AbstractVector, ishaploid::Bool)
    # fhaplo contains three possible alleles (integer codes): 0(missing), 1 , 2
    if ishaploid        
        alleleset = ["N","1","2"]
        derrule=Dict(alleleset .=> eachindex(alleleset))
        fderive = [get(derrule,i,-1) for i in fhaplo[nzorigin]]
    else
        # nfgl = length(fhaplo)      
        genoset = ["NN", "N1", "1N", "N2", "2N", "11", "12", "21", "22"] 
        derrule=Dict(genoset .=> eachindex(genoset))
        fderive = [get(derrule,join(fhaplo[i]),-1) for i in nzorigin]        
    end
    if -1 in fderive
        error("there exist unknow founder alleles; unique(fhaplo)=",unique(fhaplo))
    end
    fderive
end

function parse2_internalgeno!(resgeno::AbstractVector,format::AbstractString; ishaplo)    
    if format == "GT"                
        dict=Dict("0"=>"1","1"=>"2", "."=>"N","/"=>"","|"=>"|")
        for i in eachindex(resgeno)
            if resgeno[i] == "NA"
                resgeno[i] = ishaplo ? "N" : "NN"
            else
                g = join([get(dict,j,"2") for j = split(resgeno[i],"")])             
                if ishaplo 
                    g2 = occursin("1",g) ? (occursin("2",g) ? "N" : "1") : (occursin("2",g) ? "2" : "N")
                else
                    g2 = occursin("|",g) ? split(g,"|") : g
                end
                resgeno[i] = g2
            end
        end
    elseif format == "AD"        
        for i in eachindex(resgeno)
            if resgeno[i] in ["NA", "."]
                resgeno[i] = [0,0]
            else
                resgeno[i] = parse.(Int,split(resgeno[i],","))
            end       
        end    
    elseif format == "GP"     
        for i in eachindex(resgeno)
            if resgeno[i] in ["NA", "."]
                resgeno[i] = ishaplo ? Float32[0.5,0.5] : Float32[0.25,0.5,0.25]
            else
                resgeno[i] = parse.(Float32,split(resgeno[i],","))            
            end       
        end   
    else
        error(string("Could not call for rowgeno=",resgeno," with format=",format))  
    end
    resgeno
end

function priorfhaplo_singlesite(fgeno::AbstractVector, fformat::AbstractString;        
    foundererror::Real = 0.005,  baseerror::Real, allelicbias::Real,
    allelicoverdispersion::Real,allelicdropout::Real, israndallele::Bool)            
    if fformat=="GT_haplo"        
        alleleset = ["1","2"]
        fhaploset = [alleleset for _ in eachindex(fgeno)]
        a = foundererror
        wdict =Dict(["1"=>[(1-a),a], "2"=>[a, (1-a)], "N"=>[0.5, 0.5], missing=>[0.5, 0.5]])        
        fhaploweight = [wdict[i] for i in fgeno]
    else
        genoset = [["1","1"],["1","2"],["2","2"]]
        fhaploset = [genoset for _ in eachindex(fgeno)]
        if in(fformat, ["GT", "GT_unphased"])            
            # random allelic error model
            a = foundererror            
            if israndallele         
                pls11 = [(1-a)^2, 2a*(1-a), a^2]
                pls12 = [a*(1-a), (a^2+(1-a)^2), a*(1-a)]
                pls22 = [a^2, 2a*(1-a), (1-a)^2]
            else
                pls11 = [1-a, a/2, a/2]
                pls12 = [a/2, 1-a, a/2]
                pls22 = [a/2, a/2, 1-a]
            end
            pls1N = (pls11 .+ pls12) ./ 2
            pls2N = (pls22 .+ pls12) ./ 2
            plsNN = [0.25, 0.5, 0.25]
            wdict = Dict(["11"=>pls11, "12"=>pls12, "22"=>pls22, "1N"=>pls1N, "2N"=>pls2N, "NN"=>plsNN])
            wdict2 = Dict(["21"=>pls12, "N1"=>pls1N, "N2"=>pls2N])
            merge!(wdict,wdict2)
            gset = unique(fgeno)
            d = setdiff(gset, keys(wdict))
            if !isempty(d)
                @error string("unknown genotypes = ",d)
            end
            filter!(x->in(x.first,gset),wdict)    
            fhaploweight = [wdict[i] for i in fgeno]
        elseif fformat == "AD"
            like = MagicReconstruct.diplolikeGBS(fgeno, foundererror,baseerror,
                        allelicbias,allelicoverdispersion,allelicdropout; israndallele)
            fhaploweight = [ismissing(i) ? [0.25, 0.5, 0.25] : [i[1], i[2]+i[3], i[4]] ./ sum(i) for i in eachrow(like)]        
        elseif fformat == "GP"
            fhaploweight = fgeno
        else
            @error string("unknow kind of founder geno: ",fformat)        
        end
    end    
    fhaploset, fhaploweight
end


function callogl_singlesite!(dataprobls::AbstractVector, fhaplo::AbstractVector, offgeno::AbstractVector;
    popmakeup::AbstractDict,
    popidls=keys(popmakeup),
    epsf::Real,
    epso::Real,
    epso_perind::Union{Nothing,AbstractVector}=nothing, 
    baseerror::Real,
    allelicbias::Real,
    allelicoverdispersion::Real,
    allelicdropout::Real,
    israndallele::Bool,
    offspringformat::AbstractString)     
    issiteGT = occursin(r"^GT",offspringformat)    
    sitefhaplo = reduce(vcat,fhaplo)     
    sitemultiallele = 2
    MagicReconstruct.calsitedataprob_singlephase_nucleo!(dataprobls,sitefhaplo,offgeno, sitemultiallele, popmakeup; 
        epsf, epso,epso_perind, baseerror,allelicbias,
        allelicoverdispersion, allelicdropout,israndallele,issiteGT)
    singlelogl = 0.0
    for popid in popidls
        offls = popmakeup[popid]["offspring"]
        initprob  =  popmakeup[popid]["initprob"]
        singlelogl += sum(log(dot(initprob, dataprobls[i])) for i in offls)
    end
    singlelogl
end


function callogl_singlesite(inputfhaplo::AbstractVector,offgeno::AbstractVector; 
    popmakeup::AbstractDict,
    popidls=keys(popmakeup),
    epsf::Real,
    epso::Real,    
    baseerror::Real,
    allelicbias::Real,
    allelicoverdispersion::Real,
    allelicdropout::Real,
    israndallele::Bool,
    offspringformat::AbstractString)  
    only(callogl_singlesite_multiphase([inputfhaplo], offgeno; 
        popmakeup, popidls, epsf, epso, baseerror,
        allelicbias, allelicoverdispersion, 
        allelicdropout,offspringformat,israndallele)) 
end

function callogl_singlesite_multiphase(inputfhaplols::AbstractVector,offgeno::AbstractVector; 
    popmakeup::AbstractDict,
    popidls=keys(popmakeup),
    epsf::Real,
    epso::Real,    
    baseerror::Real,
    allelicbias::Real,
    allelicoverdispersion::Real,
    allelicdropout::Real,
    israndallele::Bool,
    offspringformat::AbstractString)  
    fhaplols = [reduce(vcat, i) for i in inputfhaplols]
    res = zeros(length(fhaplols))
    for popid in popidls
        offls = popmakeup[popid]["offspring"]        
        nzorigin =  popmakeup[popid]["nzorigin"]
        ishaploid =  popmakeup[popid]["ishaploid"]		                
        initprob = popmakeup[popid]["initprob"]       
        ismiss = [(i==[0,0] || i==[0.25,0.5,0.25]) for i in view(offgeno, offls)]        
        isnonmiss  = .!ismiss                        
        sitegeno = view(offgeno, offls[isnonmiss])           
        isempty(sitegeno) && continue   
        if ishaploid               
            if offspringformat == "AD"
                # genofrmat = AD             
                like = MagicReconstruct.haplolikeGBS(sitegeno, epso,baseerror)            
            elseif offspringformat == "GT"
                    # format = GP    
                like = MagicReconstruct.haplolike_GT(sitegeno,epso)
            else        
                # genoformat = GP (genotpe prob), GL (log10 scaled likeliihood)
                like = MagicReconstruct.haplolikeGBS(sitegeno, epso)            
            end   
            priorhaplo = MagicReconstruct.haploprior(epsf)    
            for phase in eachindex(fhaplols)
                fhaplo = fhaplols[phase]
                fderivehaplo = calfderive(fhaplo, nzorigin,ishaploid)
                ls = priorhaplo[:,fderivehaplo] * initprob
                res[phase] += sum(log(dot(ls, i)) for i in eachrow(like))   
            end
        else
            ibdbool = allequal.(nzorigin)
            nonibdbool = .!ibdbool
            initprob_ibd = initprob[ibdbool]
            initprob_nonibd = initprob[nonibdbool]
            
            if offspringformat == "AD"
                # format = AD
                like = MagicReconstruct.diplolikeGBS(sitegeno,epso,baseerror,
                        allelicbias,allelicoverdispersion,allelicdropout; israndallele)
            elseif offspringformat == "GT"
                    # format = GP    
                like = MagicReconstruct.diplolike_GT(sitegeno,epso; israndallele)
            else        
                # format = GP     
                like = MagicReconstruct.diplolikeGBS(sitegeno,epso; israndallele)
            end                                
            priordiplo = MagicReconstruct.diploprior(epsf)     
            for phase in eachindex(fhaplols)
                fhaplo = fhaplols[phase]              
                fderivediplo = calfderive(fhaplo, nzorigin,ishaploid)
                ls = priordiplo.nonibd[:,fderivediplo[nonibdbool]] * initprob_nonibd
                ls .+= priordiplo.ibd[:,fderivediplo[ibdbool]] * initprob_ibd                 
                res[phase] += sum(log(dot(ls,i)) for i in eachrow(like))                
            end
        end                
    end        
    res
end

# function get_subpop_polymorphic(fgeno::AbstractVector, offgeno::AbstractVector, popmakeup::AbstractDict;
#     minmaf::Real=0.05)
#     res = []
#     for popid in keys(popmakeup)
#         founders = popmakeup[popid]["founder"]
#         offspring = popmakeup[popid]["offspring"]
#         fmono = unique(split(join(fgeno[founders]),"")) in [["0"],["1"]]
#         if fmono
#             push!(res, popid)
#             continue
#         end
#         n1n2 = sum(i .> 0 for i in view(offgeno, offspring))
#         p = /(n1n2...)
#         if p < minmaf || p > 1 - minmaf
#             push!(res, popid)
#         end
#     end
#     setdiff(keys(popmakeup),res)
# end

function get_fixedfounders(magicped::MagicPed,offgeno::AbstractVector,offspringformat::AbstractString)
    p2off = MagicBase.get_founder2offspring(magicped;isindex = true)
    res = []
    for p in keys(p2off)
        gls = view(offgeno, p2off[p])
        gset = unique(gls) 
        # println("p=",p, ",gset=",gset)
        if offspringformat == "AD"
            if unique(gset) == [[0,0]] # no information for inferring founders
                push!(res,p)       
            else
                nreadls = reduce(vcat, gls[[i != [0,0] for i in gls]])
                if sum(nreadls .> 0) < 10 # information is too weak to infer founders
                    push!(res,p)       
                end
            end
        elseif occursin(r"^GT",offspringformat)
            alleleset = setdiff(unique(split(join(gset),"")),["|","N"])
            if alleleset  in [["0"],["1"],[]]
                push!(res,p)
            end
        else
            @warn string("unxecpted offspringformat=",offspringformat)
        end
    end
    dict = Dict(magicped.founderinfo[!,:individual] .=> 1:size(magicped.founderinfo,1))
    [dict[i] for i in res]    
end

function infer_fhaploerror_singlesite(fgeno::AbstractVector, founderformat::AbstractString,
    offgeno::AbstractVector,offspringformat::AbstractString;    
    model::AbstractString,
    popmakeup::AbstractDict,        
    isfounderinbred::Bool, 
    israndallele::Bool,    
    isinfererror::Bool,
    byfounder::Integer,
    threshcall::Real, 
    samplesize::Integer,      
    burnin::Integer,
    likeparam::LikeParam,        
    priorlikeparam::PriorLikeParam)    
    liketargetls, epsf, epso, pereoffspringerror, baseerror, allelicbias,allelicoverdispersion,allelicdropout = MagicBase.extract_likeparam(likeparam)		        
    fhaploset,fhaploweight = priorfhaplo_singlesite(fgeno,founderformat; foundererror=epsf, 
        baseerror, allelicbias,allelicoverdispersion,allelicdropout, israndallele
    )
    fhaplo = map((x,y)->x[rand(Categorical(y))], fhaploset, fhaploweight)        
    setdiff!(liketargetls, ["foundererror", "peroffspringerror"])
    if occursin(r"^GT",offspringformat)
        setdiff!(liketargetls,["baseerror","allelicbias","allelicoverdispersion","allelicdropout"])
    end
    if model == "depmodel"
        setdiff!(liketargetls,["allelicbias","allelicoverdispersion","allelicdropout"])
    end    
    isinfererror2 = isinfererror && !isempty(liketargetls) 
    findexlist = calfindexlist(byfounder,fhaploset, popmakeup; isfounderinbred)        
    fhaplo_history = []
    error_history = []    
    maxit = burnin+samplesize    
    # dataprobls = MagicReconstruct.init_dataprobls_singlephase(popmakeup)    
    for it in 1:maxit            
        if in(founderformat, ["AD"])
            fhaploset,fhaploweight = priorfhaplo_singlesite(fgeno,founderformat; foundererror=epsf, 
                baseerror, allelicbias,allelicoverdispersion,allelicdropout, israndallele
            )
        end
        update_fhaplo_singlesite!(fhaplo, fhaploset, fhaploweight, findexlist, offgeno; 
            popmakeup, epsf,epso, baseerror,allelicbias, allelicoverdispersion, 
            allelicdropout,offspringformat,israndallele
        )                    
        it > burnin && push!(fhaplo_history, isfounderinbred ? copy(fhaplo) : copy.(fhaplo))
        if isinfererror2
            for target in liketargetls
                est = first(infer_error_singlesite(target, offgeno; fhaplo, popmakeup, epsf,epso, baseerror,
                    allelicbias, allelicoverdispersion, allelicdropout, offspringformat, 
                    priorlikeparam,israndallele,temperature = it<=1 ? 0 : 1))   # temperature = 1: posterior sampling, 0 = optimize
                if target == "offspringerror"
                    epso = est
                elseif target == "baseerror"
                    baseerror = est
                elseif target == "allelicbias"
                    allelicbias = est
                elseif target == "allelicoverdispersion"
                    allelicoverdispersion = est
                elseif target == "allelicdropout"
                    allelicdropout = est
                else
                    @error string("unknown target=",target)
                end                 
                # println("it=", it,",target=",target,",est=", est, ", logl=",logl)
            end    
            nowerror = [epsf, epso, 0, baseerror,allelicbias, allelicoverdispersion, allelicdropout]
            it > burnin && push!(error_history, nowerror)                        
        end                
    end
    if isinfererror2
        esterrors = median.(eachrow(reduce(hcat,error_history)))
    else
        esterrors = [epsf, epso, 0, baseerror,allelicbias, allelicoverdispersion, allelicdropout]
    end    
    # average fhaplo    
    fhaplo_postprob = fhaplo_history_2postprob(fhaplo_history; isfounderinbred)
    fhaplo_GT, fhaplo_GP = postprob2vcfgeno(fhaplo_postprob; callthreshold=threshcall, digits=2)    
    fhaplo_GT, fhaplo_GP, esterrors
end    


function fhaplo_history_2postprob(fhaplo_history; isfounderinbred)
    genoset = isfounderinbred ? ["1","2"] : [["1","1"],["1","2"],["2","2"]]    
    n = length(fhaplo_history)
    fhaplo_GP = [[count(isequal(i), hh) for i in genoset] ./ n for hh in eachrow(reduce(hcat, fhaplo_history))]
    fhaplo_GP
end


function update_fhaplo_singlesite!(fhaplo::AbstractVector, fhaploset::AbstractVector, fhaploweight::AbstractVector,
    findexlist::AbstractVector, offgeno::AbstractVector;        
    popmakeup::AbstractDict,
    epsf::Real,
    epso::Real,
    baseerror::Real,
    allelicbias::Real,
    allelicoverdispersion::Real,
    allelicdropout::Real,
    offspringformat::AbstractString, 
    israndallele::Bool)        
    for findex in findexlist
        fhaplols = getfhaplols(findex,fhaplo,fhaploset)
        in(fhaplo, fhaplols) || @error string("current fhaplo is not in the full list!")
        logpls = callogl_singlesite_multiphase(fhaplols, offgeno; popmakeup, epsf,epso, baseerror,
            allelicbias, allelicoverdispersion, allelicdropout,offspringformat,israndallele) 
        logpls .+= [callogpri_fhaplo(fhaplo2, fhaploset, fhaploweight) for fhaplo2 in fhaplols]
        logpls .-= MagicBase.logsumexp(logpls)    
        fhaplo .= fhaplols[rand(Categorical(exp.(logpls)))]
    end
    fhaplo
end

function getfhaplols(findex::AbstractVector,fhaplo::AbstractVector,fhaploset::AbstractVector)
    set1 = vec(collect(Iterators.product(fhaploset[findex]...)))
    fhaplols = [begin
        fhaplo2 = copy(fhaplo)
        fhaplo2[findex] .= seg
        fhaplo2
    end for seg in set1]
    fhaplols
end

function callogpri_fhaplo(fhaplo::AbstractVector, fhaploset::AbstractVector, fhaploweight::AbstractVector)    
    fhaploindex = map((x,y)->findfirst(isequal(x),y), fhaplo, fhaploset)
    if in(nothing,fhaploindex)
        println("fhaploindex=",fhaploindex, ",fhaplo=",fhaplo, 
            ",fhaploset=",fhaploset,",fhaploweight=",fhaploweight)
    end
    logpri = sum(log.(map((x,y)->x[y], fhaploweight, fhaploindex)))
    logpri
end


function calfindexlist(byfounder::Integer, fhaploset::AbstractVector, popmakeup::AbstractDict; isfounderinbred)
    nfounder = length(fhaploset)
    if byfounder == -1 
        findexlist = [1:nfounder]        
    else    
        lenls = length.(fhaploset)
        nphase = sum(lenls .>= 2) > 10 ? 2^10 : reduce(*,lenls)
        if byfounder == 0 && nphase <= 3^4
            findexlist = [1:nfounder]
        else
            byfounder2 = byfounder == 0 ? (isfounderinbred ? 8 : 4) : byfounder        
            findexlist = MagicBase.getfindexlist(byfounder2,lenls, popmakeup;defaultby=2)                                
        end
    end
    findexlist
end


function cal_logpri(epsf::Real, 
    epso::Real,
    baseerror::Real, 
    allelicbias::Real,    
    allelicoverdispersion::Real,    
    allelicdropout::Real,    
    priorlikeparam::PriorLikeParam)
    logpri = 0.0    
    prior_errls = [priorlikeparam.foundererror, priorlikeparam.offspringerror, priorlikeparam.baseerror,
        priorlikeparam.allelicbias, priorlikeparam.allelicoverdispersion,priorlikeparam.allelicdropout]    
    errls = [epsf, epso, baseerror,allelicbias,allelicoverdispersion,allelicdropout]
    for i in eachindex(errls,prior_errls)
        isnothing(prior_errls[i]) && continue
        p =  errls[i]
        logpri += logpdf(prior_errls[i], p)        
    end    
    logpri    
end

function infer_error_singlesite(target::AbstractString, offgeno::AbstractVector;
    fhaplo::Union{Nothing,AbstractVector}=nothing, 
    popmakeup::AbstractDict,
    popidls=keys(popmakeup),
    epsf::Real=0.005, 
    epso::Real,
    baseerror::Real, 
    allelicbias::Real,    
    allelicoverdispersion::Real,    
    allelicdropout::Real,
    offspringformat::AbstractString,
    priorlikeparam::PriorLikeParam,
    israndallele::Bool, 
    temperature::Real=0.0)    
    accuracygoal, precisiongoal,itmax = 4, 4,20
    inputerrors = (epsf, epso, baseerror,allelicbias,allelicoverdispersion,allelicdropout)
    function loglfun(x::Real)
        epsf, epso, baseerror,allelicbias,allelicoverdispersion,allelicdropout = inputerrors
        if target == "offspringerror"
            epso = 1/(1+exp(-x))  # inverse of logit transformation     
        elseif target == "foundererror"
            epsf = 1/(1+exp(-x))                      
        elseif target == "baseerror"
            baseerror = 1/(1+exp(-x))                     
        elseif target == "allelicbias"
            allelicbias = 1/(1+exp(-x))    
        elseif target == "allelicoverdispersion"
            allelicoverdispersion = exp(x)
        elseif target == "allelicdropout"
            allelicdropout = 1/(1+exp(-x))                     
        else
            @error string("unknown genoerror type: ",target)
        end        
        if temperature ≈ 0.0
            logjacobi = 0 # adding logjacobi results too large estimate error rates and worse geneic map for testing Rice_F2 
        else
            if target in ["foundererror", "offspringerror", "baseerror","allelicbias","allelicdropout"]
                # estimation refers to transformed space
                logjacobi = x - 2*log(1+exp(x)) # dp/dx = exp(x)/(1+exp(x))^2                          
            elseif target == "allelicoverdispersion"                
                logjacobi = x   # dp/dx = exp(x) # estimation refers to transformed space       
            end    
        end
        logpri = cal_logpri(epsf, epso, baseerror,allelicbias,allelicoverdispersion,
            allelicdropout,priorlikeparam)
        if isnothing(fhaplo)
            logl = callogl_singlesite(offgeno; popmakeup, popidls, epso, baseerror,
                allelicbias, allelicoverdispersion, allelicdropout,offspringformat,israndallele) 
        else
            logl = callogl_singlesite(fhaplo, offgeno; popmakeup, popidls, epsf,epso, baseerror,
                allelicbias, allelicoverdispersion, allelicdropout,offspringformat,israndallele) 
        end        
        logl + logpri + logjacobi
    end    
    if target == "offspringerror"
        xstart = log(epso/(1-epso))
    elseif target == "foundererror"
        xstart = log(epsf/(1-epsf))
    elseif target == "baseerror"
        xstart = log(baseerror/(1-baseerror))
    elseif target == "allelicbias"
        xstart = log(allelicbias/(1-allelicbias))
    elseif target == "allelicoverdispersion"
        xstart = max(-10,log(allelicoverdispersion))
    elseif target == "allelicdropout"
        xstart = log(allelicdropout/(1-allelicdropout))
    else
        @error string("unknown genoerror type: ",target)
    end       
    if target == "allelicoverdispersion"
        lowbound, upbound = max(-12.0,xstart-10.0), xstart+10.0
    else
        fraction_upbound = 0.99
        lowbound, upbound = log(1e-5), log(fraction_upbound/(1-fraction_upbound))        
    end
    xstart = min(max(lowbound,xstart),upbound)
    if temperature ≈ 0.0
        res= MagicBase.brentMax(loglfun,lowbound,upbound;
            xstart, precisiongoal,accuracygoal,maxiter=itmax)            
        x = res[1]
        logl = res[2]
    else
        constraint = x-> lowbound < x < upbound
        res2 = MagicBase.metroplis1d(loglfun; xstart, constraint, temperature,
            stepsize = log(5.0), nstep = 5)
        x = res2[end][1]        
        logl = missing
    end
    est =  target == "allelicoverdispersion" ? exp(x) : 1.0/(1.0+exp(-x))     
    est, logl
end

function callogl_singlesite(offgeno::AbstractVector;        
    popmakeup::AbstractDict,
    popidls=keys(popmakeup),
    epso::Real,
    baseerror::Real,
    allelicbias::Real,
    allelicoverdispersion::Real,
    allelicdropout::Real,
    offspringformat::AbstractString,
    israndallele::Bool)     
    singlelogl = 0.0
    for popid in popidls
        offls = popmakeup[popid]["offspring"]
        ishaploid  =  popmakeup[popid]["ishaploid"]
        if ishaploid            
            if offspringformat == "AD"
                like = MagicReconstruct.haplolikeGBS(view(offgeno,offls),epso,baseerror)
            elseif offspringformat == "GT"
                    # format = GP    
                like = MagicReconstruct.haplolike_GT(view(offgeno,offls),epso)
            else
                offspringformat == "GP" || @error string("unexpected offspringformat = ", offspringformat) maxlog=10
                like = MagicReconstruct.haplolike(view(offgeno,offls),epso)
            end            
        else
            if offspringformat == "AD"                
                like = MagicReconstruct.diplolikeGBS(view(offgeno,offls),epso,baseerror,allelicbias,allelicoverdispersion,allelicdropout; israndallele)    
            elseif offspringformat == "GT"
                    # format = GP    
                like = MagicReconstruct.diplolike_GT(view(offgeno,offls),epso; israndallele)
            else
                offspringformat == "GP" || @error string("unexpected offspringformat = ", offspringformat) maxlog=10
                like = MagicReconstruct.diplolike(view(offgeno,offls),epso;israndallele)
            end            
        end        
        singlelogl += sum(log.(sum(like,dims=2)))
    end
    singlelogl
end

function infer_error_rawcall(offgeno::AbstractVector,offspringformat::AbstractString;    
    popmakeup::AbstractDict, 
    model::AbstractString,        
    isinfererror::Bool,        
    likeparam::LikeParam,        
    priorlikeparam::PriorLikeParam,
    israndallele::Bool)        
    liketargetls, epsf, epso, epso_perind, baseerror, allelicbias,allelicoverdispersion,allelicdropout = MagicBase.extract_likeparam(likeparam)		    
    setdiff!(liketargetls,["foundererror"])
    if occursin(r"^GT",offspringformat)
        setdiff!(liketargetls,["baseerror","allelicbias","allelicoverdispersion","allelicdropout"])
    end
    if model == "depmodel"
        setdiff!(liketargetls,["allelicbias","allelicoverdispersion","allelicdropout"])
    end        
    olderror = [epsf,epso, epso_perind, baseerror,allelicbias, allelicoverdispersion, allelicdropout]    
    olderror[1:3] .= 0.0
    isinfererror || return olderror        
    logl = oldlogl = -Inf
    maxit = 25
    for it in 1:maxit
        for target in liketargetls
            est,logl = infer_error_singlesite(target, offgeno;  popmakeup, epso, baseerror,
                allelicbias, allelicoverdispersion, allelicdropout, offspringformat, 
                priorlikeparam,israndallele,temperature = 0.0)                    
            if target == "offspringerror"
                epso = est
            elseif target == "baseerror"
                baseerror = est
            elseif target == "allelicbias"
                allelicbias = est
            elseif target == "allelicoverdispersion"
                allelicoverdispersion = est
            elseif target == "allelicdropout"
                allelicdropout = est
            else
                @error string("unknown target=",target)
            end                 
            # println("it=", it,",target=",target,",est=", est, ", logl=",logl)
        end    
        newerror = [epsf, epso, 0, baseerror,allelicbias, allelicoverdispersion, allelicdropout]
        diff = abs.(olderror .- newerror)            
        if (max(diff[1:3]...)<1e-3 && max(diff[4:6]...) < 1e-2) || (logl-oldlogl <= 0.1)
            isinfererror = false
        end
        olderror .= newerror
        oldlogl = logl        
        # @info string("it=", it, " errors=",round.(olderror,digits=5), ",logl=",oldlogl)
        if it == maxit 
           msg = string("reach maxit=", it, " errors=",round.(olderror,digits=5),",logl=",oldlogl)            
           @warn msg
        end        
        isinfererror || break
    end
    epso_perind = 0
    esterrors = [epsf,epso, epso_perind, baseerror,allelicbias, allelicoverdispersion, allelicdropout]
    esterrors
end    

function infer_offpostprob_rawcall(offgeno::AbstractVector, offspringformat::AbstractString,esterrors::AbstractVector,israndallele::Bool)
    _,epso,_, baseerror, allelicbias,allelicoverdispersion,allelicdropout = esterrors        
    if offspringformat == "AD"
        # format = AD
        like0 = MagicReconstruct.diplolikeGBS(offgeno,epso,baseerror,
            allelicbias,allelicoverdispersion,allelicdropout; israndallele)
    elseif offspringformat == "GT"
        like0 = MagicReconstruct.diplolike_GT(offgeno,epso; israndallele)
    else        
        # format = GP
        offspringformat == "GP" || @error string("unexpected offspringformat = ", offspringformat)
        like0 = MagicReconstruct.diplolikeGBS(offgeno,epso; israndallele)
    end         
    like = eachrow(like0)
    ls = unique(length.(skipmissing(like)))
    if ls == [4] || ls == []        
        # uniform prior for 4 possible genotypes
        [ismissing(i) ? [0.25, 0.5, 0.25] : normalize([i[1],i[2]+i[3],i[4]],1) for i in like]
    elseif ls == [2]
        [ismissing(i) ? [0.5,0.5] : normalize(i,1) for i in like]
    else        
        error(string("wrong like size: ",ls))
    end  
end



function parse_rowgeno!(resgeno::AbstractVector,
    rowgeno::AbstractVector,formatcode::AbstractVector,
    formatpriority::AbstractVector, missingset::AbstractVector)
    ncell = length(rowgeno)
    formatls = Vector(undef,ncell)
    if all(formatcode .< 0) #-1 denotes missing fromat
        resgeno .= "."
        commonformat = "GT"
    else
        ThreadsX.foreach(eachindex(resgeno,formatls,rowgeno)) do i
            @inbounds resgeno[i], formatls[i] = parse_vcfcell(rowgeno[i], formatcode,
                formatpriority, missingset)                
        end
        ls = skipmissing(formatls)
        if isempty(ls)
            commonformat = "GT"
        else        
            formatls2 = unique(ls)
            posls = [findfirst(formatpriority .== i[1:2]) for i in formatls2]
            commonformat =  formatls2[argmin(posls)]     
            b =[i !== commonformat for i in formatls]
            resgeno[b] .= "."           
        end        
    end    
    commonformat
end

function parse_vcfcell(cellgeno::AbstractString,formatcode::AbstractVector,
    formatpriority::AbstractVector, missingset::AbstractVector)
    n = length(formatcode)
    geno = split(cellgeno,":")
    ng = length(geno)
    if ng < n
        geno = vcat(geno, repeat(["."],n-ng))
    elseif ng > n
        geno = geno[1:n]
    end
    b = [(geno[i] in missingset) || formatcode[i] < 0 for i in 1:n]  # -1 of formatcode denotes missing fromat
    if all(b)
        genoset = geno[formatcode .>= 1]
        if in("./.", genoset)
            g = "./."
            f = "GT" 
        elseif in(".|.", genoset)
            g = ".|."
            f = "GT" 
        else
            g = "."
            f = missing
        end
    else
        # the less the formatcode element, the higher the priority
        upbound = length(formatpriority)+10
        i = argmin(@. formatcode + upbound * b)             
        g = geno[i]
        f = formatpriority[formatcode[i]]
    end
    g, f    
end
