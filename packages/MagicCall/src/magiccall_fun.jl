
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

`isfounderinbred::Bool=true`: if true, founders are inbred, and otherwise they are outbred.

`likeparameters::LikeParameters=LikeParameters(offspringerror=0.03, peroffspringerror=0.0)`: parameters for genotypic data model. 
  If isinfererror = true, parameters with values being nothing will be inferred. 

`threshlikeparameters::ThreshLikeParameters=ThreshLikeParameters()`: markers with inferred likeparameters values > threshlikeparameters values will be deleted. 

`priorlikeparameters::PriorLikeParameters=PriorLikeParameters(offspringerror=Beta(1.05,9),seqerror=Beta(1.05,9))`: priors for likelihood parameters

`israndallele::Bool=true`: if true, genotyping error model follows the random allelic model, and otherwise the random genotypic model. 

`threshcall::Real = 0.9`: offspring genotypes are call if 
  the maximum posterior probability > threshcall.

`isdelmultiallelic::Bool=true`: if true, delete markers with >=3 alleles. 

`isdelmonomorphic::Bool=true`: if true, delete monomorphic markers. 

`minmaf::Real = 0.05`: delete markers with minor allele frequency < minmaf. 

`maxmiss::Real = 0.99`: delete markers with genotype missing frequency > maxmiss. 

`israwcall::Bool= false`: if true, perform raw genotype calling. 

`isinfererror::Bool = !israwcall`: if true, infer marker specific likelihood parameters that have values of nothing in likeparameters. 

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
    likeparameters::LikeParameters=LikeParameters(offspringerror=0.03, peroffspringerror=0.0),    
    threshlikeparameters::ThreshLikeParameters=ThreshLikeParameters(),    
    priorlikeparameters::PriorLikeParameters=PriorLikeParameters(offspringerror=Beta(1.05,9),seqerror=Beta(1.05,9)),                          
    threshcall::Real = 0.9,
    israwcall::Bool= false, 
    byfounder::Integer=1,
    isdelmultiallelic::Bool=true,
    isdelmonomorphic::Bool=true,    
    minmaf::Real = 0.05, # set monomorphic subpopulation to missing if maf < minmaf
    maxmiss::Real = 0.99,         
    isinfererror::Bool = !israwcall, 
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
        "likeparameters = ", likeparameters, "\n",
        "threshlikeparameters = ", threshlikeparameters, "\n",		
        "priorlikeparameters = ", priorlikeparameters, "\n",	        
        "isfounderinbred = ", isfounderinbred, "\n",        
        "israndallele = ", israndallele, "\n",        
        "byfounder = ", byfounder, "\n",        
        "isdelmultiallelic = ", isdelmultiallelic, "\n",
        "isdelmonomorphic = ", isdelmonomorphic, "\n",        
        "minmaf = ", minmaf, "\n",        
        "maxmiss = ", maxmiss, "\n",        
        "threshcall = ", threshcall, "\n",        
        "israwcall = ", israwcall, "\n",        
        "isinfererror = ", isinfererror, "\n",                
        "isparallel = ", isparallel, isparallel ? string("(nworker=",nworkers(),")") : "", "\n",
        "outstem = ", outstem, "\n",
        "outext = ", outext, "\n",
        "logfile = ", logfile, "\n",
        "workdir = ", workdir, "\n",
        "verbose = ", verbose)
    printconsole(logio,verbose,msg)    
    if !israwcall
        msg = string("likeparameters = (foundererror, offspringerror, peroffspringerror, seqerror, allelebalancemean, allelebalancedisperse, alleledropout)")
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
        seglen = min(100,div(nsnp, nwork))
        paragraphls = collect(Iterators.partition(1:nsnp,seglen*nwork))
        msg *= string(" divided into ", length(paragraphls), " blocks")
        nwork > 1 && (msg *= string("; each block distributed over ", nwork, " workers"))        
    else
        paragraphls = nothing
    end    
    printconsole(logio,verbose,msg)   
    outstem *= "_magiccall"
    in_open(genofile2,"r") do inio                
        out_open = outext in [".vcf.gz",".csv.gz"] ? GZip.open : open
        outfile = outstem*"_geno"*outext
        outfile2 = getabsfile(workdir,outfile)            
        delmarkerfile = outstem*"_delmarker.vcf.gz"
        delmarkerfile2 = getabsfile(workdir,delmarkerfile)            
        out_open(outfile2, "w") do outio                
            GZip.open(delmarkerfile2,"w") do delio
                magiccall_io(inio, outio, delio, logio, pedinfo,model, 
                    likeparameters, threshlikeparameters, priorlikeparameters, 
                    formatpriority, israndallele,isfounderinbred, byfounder, 
                    isdelmultiallelic, isdelmonomorphic, minmaf,maxmiss, 
                    threshcall, israwcall, isinfererror,
                    paragraphls,commentstring, nheader, workdir, verbose)            
            end
        end
        msg = string("save called_genofile: ", outfile)
        printconsole(logio,verbose,msg)        
        msg = string("save delmarker_genofile: ", delmarkerfile)
        printconsole(logio,verbose,msg)       
        if !israwcall
            for (i,mapfile) in enumerate([outfile2])                
                try 
                    figerr = plotmarkererror(mapfile;tukeyfence=3.0, workdir)
                    if !isnothing(figerr)                        
                        markererrfile = string(outstem, i==1 ? "" : "_delmarker", "_inferred_error.png")
                        try 
                            MagicBase.savefig(figerr,getabsfile(workdir,markererrfile))
                        catch err
                            msg = string("Could not savefig. ",err)
                            @warn msg
                            printconsole(logio, verbose, "Warning: "*msg)
                        end
                    end 
                catch err                
                    msg = string(err, ". Could not plot markererror for mapfile= ",mapfile)
                    @warn msg
                    printconsole(logio, verbose, "Warning: "*msg)
                end
            end
        end
    end
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, logio, tused,"magiccall"; verbose,delim="-")
    return nothing
end

function magiccall_io(inio::IO,outio::IO,delio::IO,
    logio::Union{Nothing, IO},
    pedinfo::Union{Integer,AbstractString}, 
    model::AbstractString, 
    likeparameters::LikeParameters,        
    threshlikeparameters::ThreshLikeParameters,        
    priorlikeparameters::PriorLikeParameters,        
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
            res = magiccall_rowgeno(rowstring,fcols,offcols,magicped; israwcall, isdelmultiallelic, isdelmonomorphic,
                israndallele, isfounderinbred,byfounder,model,popmakeup, formatpriority,
                minmaf,maxmiss,threshcall, isinfererror,
                likeparameters, threshlikeparameters, priorlikeparameters, missingset, nstate,nfgl)        
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
            res = pmap(x-> magiccall_rowgeno(x,fcols,offcols,magicped; israwcall, isdelmultiallelic, isdelmonomorphic,
                israndallele,isfounderinbred,byfounder, model, popmakeup, formatpriority, 
                minmaf, maxmiss, threshcall, isinfererror,
                likeparameters, threshlikeparameters, priorlikeparameters, missingset, nstate,nfgl), multirows)            
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
    elseif largeerrid == "largeerror_seqerror"
        nlargeerr[4] += 1
    elseif largeerrid == "largeerror_allelebalancemean"
        nlargeerr[5] += 1
    elseif largeerrid == "largeerror_allelebalancedisperse"
        nlargeerr[6] += 1
    elseif largeerrid == "largeerror_alleledropout"
        nlargeerr[7] += 1
    end
    nlargeerr
end

function magiccall_rowgeno(rowstring::AbstractString,
    fcols::AbstractVector,offcols::AbstractVector,
    magicped::MagicPed;
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
    likeparameters::LikeParameters,    
    threshlikeparameters::ThreshLikeParameters,    
    priorlikeparameters::PriorLikeParameters,    
    missingset::AbstractVector,    
    nstate::Integer, 
    nfgl::Integer)
    rowgeno = split(rowstring,"\t")
    ismultiallele = length(split(rowgeno[5],",")) > 1 # col5=alternative                     
    ismultiallele && isdelmultiallelic && return ("multiallelic", rowstring)
    res = transform_rowgeno!(rowgeno, fcols,offcols,likeparameters,formatpriority,
        missingset,isfounderinbred,isdelmonomorphic)    
    res == "monomorphic" && return ("monomorphic",join(rowgeno,"\t"))
    inputformat, fgeno, offgeno, offspringformat = res
    isoffmiss = rowgeno[offcols] .== "NA"
    # imputation/correction and error estimations    
    if israwcall
        fhaplo = fgeno
        esterrors = MagicCall.infer_singlesite_rawcall(offgeno,offspringformat; popmakeup, model,isinfererror,
            likeparameters,priorlikeparameters,israndallele)
    else
        fixedfounders = get_fixedfounders(magicped,offgeno,offspringformat)     
        fhaplo, esterrors = infer_singlesite(fgeno,fixedfounders, offgeno,offspringformat; 
            model, popmakeup, israndallele,isfounderinbred, isinfererror, byfounder, 
            likeparameters, priorlikeparameters)    
    end
    newfgeno = join.(fhaplo)
    calledgeno2vcf!(newfgeno; ishaplo = isfounderinbred)       
    add_error_info!(rowgeno,esterrors)    
    ismono = in(unique(fhaplo),[["1"],["2"]])
    if isdelmonomorphic && ismono        
        if !israwcall
            rowgeno[8] *= string(";FOUNDERHAPLO=",replace(join(newfgeno),"/"=>""))
        end
        return ("monomorphic",join(rowgeno,"\t"))
    end    
    islargeerror, largeerrorid = get_isdelmarker(esterrors,threshlikeparameters)
    if islargeerror         
        if !israwcall
            rowgeno[8] *= string(";FOUNDERHAPLO=",replace(join(newfgeno),"/"=>""))
        end        
        return ("largeerror_"*largeerrorid,join(rowgeno,"\t"))
    end        
    resid = ismultiallele ? "multiallelic" : (ismono ? "monomorphic" : "biallelic")
    # update rowgeno for founder geno         
    outformat =  inputformat
    for i in fcols
        rowgeno[i] == "NA" && (rowgeno[i] ="." )
    end     
    if in(:GT, outformat)
        i_gt = findfirst(isequal(:GT),outformat)
        for i in eachindex(fcols)
            gg = split(rowgeno[fcols[i]],":")
            gg[i_gt] = newfgeno[i]
            rowgeno[fcols[i]] = join(gg,":")
        end
        offspringformat == "GT" && return (resid,join(rowgeno,"\t"))
    else
        pushfirst!(outformat,:GT)
        for i in eachindex(fcols)
            rowgeno[fcols[i]] = string(newfgeno[i], ":", rowgeno[fcols[i]])
        end
        for i in eachindex(offcols)
            rowgeno[offcols[i]] = string("./.:", rowgeno[offcols[i]])
        end
        # To add GT for offcols
        if offspringformat == "GT" 
            @error string("inconsistent among offspringformat=",offspringformat,", FORMAT=",rowgeno[9]) 
        end
    end    
    # update rowgeno for offspring geno
    if israwcall
        offpostprob = singlesite_offpostprob_raw(offgeno,offspringformat, esterrors,israndallele)
    else
        offpostprob = singlesite_offpostprob(fhaplo,offgeno, offspringformat, popmakeup,esterrors,nstate,nfgl, israndallele)
    end
    newoffgeno = MagicBase.callfromprob.(offpostprob, threshcall; isphased=false,ishalfcall=true)        
    
    geno_mono2miss!(newoffgeno, offpostprob, popmakeup; minmaf = minmaf)
    b = [occursin("N",i) for i in newoffgeno]    
    noff = length(newoffgeno)
    if maxmiss < 1 && all(b)
        resid = "maxmiss"
    else
        alleles = split(join(newoffgeno),"")
        deleteat!(alleles, alleles .== "N")
        n1 = sum(alleles .== "1")
        n2 = sum(alleles .== "2")
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
    end
    calledgeno2vcf!(newoffgeno; ishaplo = false)            
    if in(:GT, outformat)
        i_gt = findfirst(isequal(:GT),outformat)
        for i in eachindex(offcols)
            gg = split(rowgeno[offcols[i]],":")            
            gg[i_gt] = newoffgeno[i]
            rowgeno[offcols[i]] = join(gg,":")
        end
    else
        @error string("unexpected outformat=",outformat, " for rowgeno=",rowgeno)
    end
    if in(:GP, outformat)
        i_gp = findfirst(isequal(:GP),outformat)
        for i in eachindex(offcols)
            gg = split(rowgeno[offcols[i]],":")
            gg[i_gp] = join(round.(offpostprob[i],digits=5),",")
            rowgeno[offcols[i]] = join(gg,":")
        end
    else
        push!(outformat,:GP)
        for i in eachindex(fcols)
            rowgeno[fcols[i]] = string(rowgeno[fcols[i]],":.")
        end
        for i in eachindex(offcols)
            rowgeno[offcols[i]] = string(rowgeno[offcols[i]], ":", join(round.(offpostprob[i],digits=5),","))
        end        
    end    
    rowgeno[offcols[isoffmiss]] .= "."
    rowgeno[9] = join(outformat,":")
    (resid, join(rowgeno,"\t"))
end

function transform_rowgeno!(rowgeno::AbstractVector,fcols::AbstractVector,offcols::AbstractVector,
    likeparameters::LikeParameters,
    formatpriority::AbstractVector,
    missingset::AbstractVector,    
    isfounderinbred::Bool,
    isdelmonomorphic::Bool)
    keepvcf = true    
    ids = Symbol.(formatpriority)
    vals = 1:length(formatpriority)
    prioritytuple = (; zip(ids, vals)...)
    inputformat = Symbol.(split(rowgeno[9],":"))
    formatcode = [get(prioritytuple,i,-1) for i in inputformat] # -1 denotes missing fromat            
    # fgeno
    fgeno = Vector(undef, length(fcols))     
    founderformat = MagicBase.parse_rowgeno!(fgeno,rowgeno[fcols],
        formatcode, formatpriority, missingset,keepvcf)        
    ismono = in(unique(split(replace(join(fgeno),"/"=>"", "|"=>""),"")),[["0"],["1"]])
    if isdelmonomorphic && ismono        
        return "monomorphic"
    end  
    parse2_internalgeno!(fgeno,founderformat; ishaplo=isfounderinbred)    
    threshcall = 0.95 # fixed threshold for founders
    if founderformat == "AD"
        seqerror = MagicBase.get_seqerror(likeparameters)
        if isfounderinbred
            fgeno .= MagicBase.callfromprob.(MagicBase.genoprobhaplo.(fgeno, seqerror),threshcall; 
                isphased=false,ishaplo=true)
        else
            fgeno .= MagicBase.callfromprob.(MagicBase.genoprobdiplo.(fgeno, seqerror),threshcall; 
                isphased=false,ishaplo=false)
        end
    elseif founderformat == "GP"
        fgeno = MagicBase.callfromprob.(fgeno, threshcall; isphased=false,ishaplo=isfounderinbred)
    end    
    ismono = in(unique(split(join(fgeno),"")),[["1"],["2"]])
    if isdelmonomorphic && ismono        
        return "monomorphic"
    end      
    # set offspring geno
    offgeno = Vector(undef, length(offcols))
    offspringformat = MagicBase.parse_rowgeno!(offgeno, rowgeno[offcols],
        formatcode, formatpriority, missingset,keepvcf)            
    parse2_internalgeno!(offgeno,offspringformat; ishaplo=false)   
    alleleset = unique(split(join(offgeno),""))
    setdiff!(alleleset,["N"])
    ismono = in(alleleset,[["1"],["2"]])
    if isdelmonomorphic && ismono        
        return "monomorphic"
    end      
    inputformat, fgeno,offgeno,offspringformat
end

function get_isdelmarker(esterrors,threshlikeparameters)    
    epsf,epso, _,seqerror, allelebalancemean,allelebalancedisperse,alleledropout = esterrors
    epsf > threshlikeparameters.foundererror && return (true, "foundererror")
    epso > threshlikeparameters.offspringerror && return (true, "offspringerror")
    seqerror > threshlikeparameters.seqerror && return (true, "seqerror")
    thresh = threshlikeparameters.allelebalancemean
    (1-thresh <= allelebalancemean <= thresh) || return (true, "allelebalancemean")
    allelebalancedisperse > threshlikeparameters.allelebalancedisperse && return (true, "allelebalancedisperse")
    alleledropout > threshlikeparameters.alleledropout && return (true, "alleledropout")
    (false, "")        
end

function add_error_info!(rowgeno::AbstractVector,esterrors)
    epsf,epso, _, seqerror, allelebalancemean,allelebalancedisperse,alleledropout = esterrors
    info = string("FOUNDERERROR=",round(epsf,digits = 5))
    info *= string(";OFFSPRINGERROR=",round(epso,digits = 5))
    info *= string(";SEQERROR=",round(seqerror,digits = 5))
    info *= string(";ALLELEBALANCEMEAN=",round(allelebalancemean,digits = 5))
    info *= string(";ALLELEBALANCEDISPERSE=",round(allelebalancedisperse,digits = 5))
    info *= string(";ALLELEDROPOUT=",round(alleledropout,digits = 5))    
    inputinfo = strip(rowgeno[8])
    if !in(inputinfo, [".",""])
        info = string(inputinfo, ";", info)
    end
    rowgeno[8] = info    
    info
end

function geno_mono2miss!(calledgeno::AbstractVector,offpostprob::AbstractVector, popmakeup::AbstractDict; minmaf=0.05)    
    for popid in keys(popmakeup)        
        offls = popmakeup[popid]["offspring"]
        alleles = split(join(view(calledgeno,offls)),"")
        deleteat!(alleles, alleles .== "N")
        n1 = sum(alleles .== "1")
        n2 = sum(alleles .== "2")
        ishaploid = popmakeup[popid]["ishaploid"]
        missprob = ishaploid ? [0.5, 0.5] : [0.25,0.5,0.25]
        if n1 + n2 < 10
            calledgeno[offls] .= "NN"
            for i in offls
                offpostprob[i] .= missprob
            end
        else
            if n1/(n1+n2) < minmaf
                calledgeno[offls] .= "NN"
                for i in offls
                    offpostprob[i] .= missprob
                end
            end
        end
    end
    calledgeno
end
function calledgeno2vcf!(calledgeno::AbstractVector; ishaplo)
    if ishaplo
        dict = Dict("1"=>"0/0","2"=>"1/1","N"=>"./.")        
    else
        dict = Dict("11"=>"0/0","12"=>"0/1","22"=>"1/1","1N"=>"0/.","2N"=>"1/.","NN"=>"./.")
    end 
    for i in eachindex(calledgeno)
        calledgeno[i] = dict[calledgeno[i]] 
    end   
    calledgeno
end


function singlesite_offpostprob(inputfhaplo::AbstractVector,offgeno::AbstractVector, offspringformat::AbstractString, 
    popmakeup::AbstractDict, esterrors::Tuple,     
    nstate::Integer,nfgl::Integer,israndallele::Bool)    
    epsf,epso, _, seqerror, allelebalancemean,allelebalancedisperse,alleledropout = esterrors
    fhaplo = reduce(vcat, inputfhaplo)
    res = Vector(undef, length(offgeno))
    for popid in keys(popmakeup)        
        offls = popmakeup[popid]["offspring"]
        nzstate =  popmakeup[popid]["nzstate"]
        ishaploid =  popmakeup[popid]["ishaploid"]		        
        initprob = spzeros(nstate)
        initprob[nzstate] .= popmakeup[popid]["initprob"]       
        ismiss = [(i==[0,0] || i==[0.25,0.5,0.25]) for i in view(offgeno, offls)]
        isnonmiss  = .!ismiss                        
        sitegeno = view(offgeno, offls[isnonmiss])              
        if ishaploid                                         
            if !isempty(sitegeno)
                fderivehaplo = MagicReconstruct.calfderive(reshape(fhaplo,1,:); nzstate,ishaploid)[1,:]		                        
                if offspringformat == "AD"
                    # format = AD
                    like = MagicReconstruct.haplolikeGBS(sitegeno,epso,seqerror)
                else
                    # format = GP or GT
                    like = MagicReconstruct.haplolikeGBS(sitegeno,epso)
                end       
                priorhaplo = MagicReconstruct.haploprior(epsf)                          
                ls = priorhaplo[:,fderivehaplo[nzstate]] * initprob[nzstate]
                post = ls' .* like
                res[offls[isnonmiss]] .= [normalize(i,1) for i in eachrow(post)]   
            end
            if any(ismiss)
                res[offls[ismiss]]  .= [[0.5,0.5] for _ in 1:sum(ismiss)]        
            end
        else
            nfgl == length(fhaplo) || @error string("inconsistent nfgl=",nfgl, ", nfgl2=",length(fhaplo))
            nstate == nfgl^2 || error(string("mismatch between nfgl=",nfgl, " and nstate=",nstate))					       
            if !isempty(sitegeno)                        
                zzstate = setdiff(1:nstate, nzstate)					
                ibdbool= falses(nfgl^2)
                ibdbool[[(i-1)nfgl+i for i=1:nfgl]] .= true
                nonibdbool = .!ibdbool
                ibdbool[zzstate] .= false
                nonibdbool[zzstate] .= false		
                initprob_ibd = initprob[ibdbool]
                initprob_nonibd = initprob[nonibdbool]
                fderivediplo = MagicReconstruct.calfderive(reshape(fhaplo,1,:); nzstate,ishaploid)[1,:]		
                if eltype(first(sitegeno)) <: Integer
                    # format = AD
                    like = MagicReconstruct.diplolikeGBS(sitegeno,epso,seqerror,
                        allelebalancemean,allelebalancedisperse,alleledropout; israndallele)
                else
                    # format = GP
                    like = MagicReconstruct.diplolikeGBS(sitegeno,epso; israndallele)
                end                 
                priordiplo = MagicReconstruct.diploprior(epsf)                       
                ls = priordiplo.nonibd[:,fderivediplo[nonibdbool]] * initprob_nonibd
                ls .+= priordiplo.ibd[:,fderivediplo[ibdbool]] * initprob_ibd 
                post = ls' .* like
                res[offls[isnonmiss]] .= [begin 
                    v = normalize(i,1)
                    [v[1],v[2]+v[3],v[4]]
                end for i in eachrow(post)]
            end
            if any(ismiss)
                res[offls[ismiss]] .= [[0.25,0.5,0.25] for _ in 1:sum(ismiss)]
            end
        end                
    end        
    res
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
            if resgeno[i] in ["NA","."]
                resgeno[i] = [0,0]
            else
                resgeno[i] = parse.(Int,split(resgeno[i],","))
            end       
        end    
    elseif format == "GP"     
        for i in eachindex(resgeno)
            if resgeno[i] == "NA"
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

function priorfhaplo_singlesite(fgeno::AbstractVector, fformat::AbstractString,
    fixedfounderes::AbstractVector;     
    genoerror::Real = 0.05,
    allowmissing::Bool=true,
    allowcorrect::Bool=true)
    # random genotype error model
    if fformat=="GT_haplo"
        if allowcorrect
            if allowmissing
                dict = Dict(["1"=>["1","2","N"],"2"=>["1","2","N"], "N"=>["1","2","N"]])        
                wdict =Dict(["1"=>[1-genoerror,genoerror/2, genoerror/2],
                    "2"=>[genoerror/2, 1-genoerror, genoerror/2], 
                    "N"=>[0.5-genoerror/2, 0.5-genoerror/2, genoerror]])         
            else
                dict = Dict(["1"=>["1","2"],"2"=>["1","2"], "N"=>["1","2"]])        
                wdict =Dict(["1"=>[1-genoerror, genoerror],"2"=>[genoerror,1-genoerror], "N"=>[0.5,0.5]])         
            end
        else
            if allowmissing
                dict = Dict(["1"=>["1","N"],"2"=>["2","N"], "N"=>["1","2","N"]])       
                # genoerror denotes prob of missing 
                wdict =Dict(["1"=>[1-genoerror,genoerror],
                    "2"=>[1-genoerror, genoerror], 
                    "N"=>[0.5-genoerror/2, 0.5-genoerror/2, genoerror]])         
            else
                dict = Dict(["1"=>["1"],"2"=>["2"], "N"=>["1","2"]])        
                wdict =Dict(["1"=>[1],"2"=>[1], "N"=>[0.5,0.5]])         
            end
        end
    elseif fformat=="GT_unphased"
        if allowcorrect
            if allowmissing
                # At most one allelic error, and at most one becomes missing
                dict = Dict(["11"=>[["1","1"],["1","2"],["1","N"]],
                    "12"=>[["1","2"],["1","1"],["2","2"],["1","N"],["2","N"]],
                    "22"=>[["2","2"],["1","2"],["2","N"]],
                    "NN"=>[["1","1"],["1","2"],["2","2"],["1","N"],["2","N"],["N","N"]]])            
                wdict = Dict(["11"=>[1-genoerror,genoerror/2, genoerror/2],
                    "12"=>[1-genoerror, genoerror/4, genoerror/4, genoerror/4, genoerror/4],
                    "22"=>[1-genoerror, genoerror/2, genoerror/2],
                    "NN"=>[0.25-genoerror/3,0.5-genoerror/3,0.25-genoerror/3,genoerror/3,genoerror/3,genoerror/3]])
            else
                # At most one allelic error; unphased/unordered genotypes for single site analysis
                dict = Dict(["11"=>[["1","1"],["1","2"]],
                    "12"=>[["1","2"],["1","1"],["2","2"]],
                    "22"=>[["2","2"],["1","2"]],
                    "NN"=>[["1","1"],["1","2"],["2","2"]]])            
                wdict = Dict(["11"=>[1-genoerror,genoerror],
                    "12"=>[1-genoerror, genoerror/2, genoerror/2],
                    "22"=>[1-genoerror, genoerror],
                    "NN"=>[0.25,0.5,0.25]])
            end
        else
            if allowmissing                
                # genoerror denotes prob of missing 
                dict = Dict(["11"=>[["1","1"],["1","N"]],
                    "12"=>[["1","2"],["1","N"],["2","N"]],
                    "22"=>[["2","2"],["2","N"]],
                    "NN"=>[["1","1"],["1","2"],["2","2"],["1","N"],["2","N"],["N","N"]]])            
                wdict = Dict(["11"=>[1-genoerror,genoerror],
                    "12"=>[1-genoerror, genoerror/2, genoerror/2],
                    "22"=>[1-genoerror, genoerror],
                    "NN"=>[0.25-genoerror/3,0.5-genoerror/3,0.25-genoerror/3,genoerror/3,genoerror/3,genoerror/3]])
            else                
                dict = Dict(["11"=>[["1","1"]],
                    "12"=>[["1","2"]],
                    "22"=>[["2","2"]],
                    "NN"=>[["1","1"],["1","2"],["2","2"]]])            
                wdict = Dict(["11"=>[1],
                    "12"=>[1],
                    "22"=>[1],
                    "NN"=>[0.25,0.5,0.25]])
            end
        end
        gset = unique(fgeno)
        d = setdiff(gset, keys(dict))
        if !isempty(d)
            if allowmissing
                dict2 = Dict(["1N"=>[["1","1"],["1","2"],["2","2"],["1","N"],["2","N"],["N","N"]], 
                    "2N"=>[["1","1"],["1","2"],["2","2"],["1","N"],["2","N"],["N","N"]]])
                pls1N = [0.5-genoerror/2,0.5-genoerror/2,genoerror/4,genoerror/4,genoerror/4,genoerror/4]
                pls2N = [genoerror/4,0.5-genoerror/2,0.5-genoerror/2,genoerror/4,genoerror/4,genoerror/4]
                wdict2 = Dict(["1N"=>pls1N,"2N"=>pls2N])
            else
                dict2 = Dict(["1N"=>[["1","1"],["1","2"],["2","2"]],
                    "2N"=>[["1","1"],["1","2"],["2","2"]]])
                wdict2 = Dict(["1N"=>[1-genoerror,genoerror/2,genoerror/2],
                    "2N"=>[genoerror/2, genoerror/2, 1-genoerror]])
            end
            merge!(dict,dict2)
            merge!(wdict,wdict2)
            d = setdiff(gset, keys(dict))
            isempty(d) || @error string("unexpected genotypes: ",d)
            filter!(x->in(x.first,gset),dict)
        end
    else        
        @error string("unknow kind of founder geno: ",fformat)        
    end    
    fhaploset = [if in(i, fixedfounderes) 
        fformat=="GT_haplo" ? split(fgeno[i],"") : [split(fgeno[i],"")]        
    else 
        get(dict, fgeno[i], missing) 
    end for i in eachindex(fgeno)]    
    b = ismissing.(fhaploset)
    any(b) && @error string("unknown genotypes = ",fgeno[b])    
    fhaploweight = [in(i, fixedfounderes) ? [1.0] : get(wdict, fgeno[i], missing) for i in eachindex(fgeno)]    
    b = ismissing.(fhaploweight)
    any(b) && @error string("unknown genotypes = ",fgeno[b])
    # if !all(length.(fhaploset) .== length.(fhaploweight))
    #     println("fhaploset=",fhaploset, ", length=",length.(fhaploset))
    #     println("fhaploweight=",fhaploweight,", length=",length.(fhaploweight))
    # end
    fhaploset, fhaploweight
end

# function calfhaploset_singlesite(fgeno::AbstractVector, fformat::AbstractString,vcfformat; 
#     isfounderinbred::Bool,
#     allowmissing::Bool=true)
#     if fformat=="GT_haplo"
#         if allowmissing
#             dict = Dict(["1"=>["1","N"],"2"=>["2","N"], "N"=>["1","2","N"]])        
#         else
#             dict = Dict(["1"=>["1"],"2"=>["2"], "N"=>["1","2"]])        
#         end
#     elseif fformat=="GT_unphased"
#         if allowmissing
#             dict = Dict(["11"=>[["1","1"],["1","N"]],"12"=>[["1","2"],["1","N"],["2","N"]],
#                 "22"=>[["2","2"],["2","N"]],
#                 "NN"=>[["1","1"],["1","2"],["2","2"],["1","N"],["2","N"],["N","N"]]])
#         else
#             dict = Dict(["11"=>[["1","1"]],"12"=>[["1","2"]],"22"=>[["2","2"]],"NN"=>[["1","1"],["1","2"],["2","2"]]])
#         end
#         gset = unique(fgeno)
#         d = setdiff(gset, keys(dict))
#         if !isempty(d)
#             if allowmissing
#                 dict2 = Dict(["1N"=>[["1","1"],["1","2"],["1","N"]],"N1"=>[["1","1"],["1","2"],["1","N"]],                    
#                     "21"=>[["1","2"],["1","N"],["2","N"]],
#                     "2N"=>[["2","2"],["1","2"],["2","N"]],"N2"=>[["2","2"],["1","2"],["2","N"]]])
#             else
#                 dict2 = Dict(["1N"=>[["1","1"],["1","2"]],"N1"=>[["1","1"],["1","2"]],
#                     "21"=>[["1","2"]],
#                     "2N"=>[["2","2"],["1","2"]],"N2"=>[["2","2"],["1","2"]]])
#             end
#             merge!(dict,dict2)
#             d = setdiff(gset, keys(dict))
#             isempty(d) || @error string("unexpected genotypes: ",d)
#             filter!(x->in(x.first,gset),dict)
#         end
#     else        
#         @error string("unknow kind of founder geno: ",fformat)        
#     end    
#     fhaploset = [get(dict, i, missing) for i in fgeno]    
#     b = ismissing.(fhaploset)
#     any(b) && @error string("unknown genotypes = ",fgeno[b])
#     # account for loss of heterozygosity due to allele balance bias
#     if vcfformat == "AD" && !isfounderinbred        
#         if allowmissing
#             fhaploset = [i==[["1","1"],["1","N"]] ? [["1","1"],["1","2"],["1","N"]] : (i==[["2","2"],["2","N"]] ? [["2","2"],["1","2"],["2","N"]] : i) for i in fhaploset]
#         else
#             fhaploset = [i==[["1","1"]] ? [["1","1"],["1","2"]] : (i==[["2","2"]] ? [["2","2"],["1","2"]] : i) for i in fhaploset]
#         end
#     end
#     fhaploset
# end


function callogl_singlesite!(dataprobls::AbstractVector, fhaplo::AbstractVector, offgeno::AbstractVector;
    popmakeup::AbstractDict,
    popidls=keys(popmakeup),
    epsf::Real,
    epso::Real,
    epso_perind::Union{Nothing,AbstractVector}=nothing, 
    seqerror::Real,
    allelebalancemean::Real,
    allelebalancedisperse::Real,
    alleledropout::Real,
    israndallele::Bool,
    offspringformat::AbstractString)     
    issiteGT = occursin(r"^GT",offspringformat)    
    sitefhaplo = reduce(vcat,fhaplo)    
    MagicReconstruct.calsitedataprob_singlephase!(dataprobls,sitefhaplo,offgeno,popmakeup; epsf,epso,epso_perind, seqerror,allelebalancemean,
        allelebalancedisperse,alleledropout,israndallele,issiteGT)
    singlelogl = 0.0
    for popid in popidls
        offls = popmakeup[popid]["offspring"]
        initprob  =  popmakeup[popid]["initprob"]
        singlelogl += sum(log(dot(initprob, dataprobls[i])) for i in offls)
    end
    singlelogl
end




function callogl_singlesite2(fhaplo::AbstractVector, offgeno::AbstractVector;
    popmakeup::AbstractDict,
    epsf::Real,
    epso::Real,
    epso_perind::Union{Nothing,AbstractVector}=nothing, 
    seqerror::Real,
    allelebalancemean::Real,
    allelebalancedisperse::Real,
    alleledropout::Real,
    offspringformat::AbstractString,
    israndallele::Bool)     
    issiteGT = occursin(r"^GT",offspringformat)
    isoffphased = false
    fderive, offcode = MagicReconstruct.precompute_chr(reshape(reduce(vcat,fhaplo),1,:),
        reshape(offgeno,1,:),popmakeup, isoffphased, [issiteGT])
    singlelogl = 0.0
    nsnp = size(offcode,1)
    for popid in keys(popmakeup)
        ishaploid = popmakeup[popid]["ishaploid"]
        offls = popmakeup[popid]["offspring"]
        nzstate = popmakeup[popid]["nzstate"]
        initprob  =  popmakeup[popid]["initprob"]
        dataprobseq = [zeros(MagicReconstruct._float_like,length(nzstate))  for _ in 1:nsnp]
        for off = offls
            obsseq = offcode[:,off]                        
            MagicReconstruct.caldataprobseq!(dataprobseq,obsseq,epsf,epso,epso_perind[off], seqerror,
                allelebalancemean,allelebalancedisperse,alleledropout,
                fderive,nzstate, isoffphased, israndallele, [issiteGT],ishaploid)
            singlelogl += log(dot(first(dataprobseq),initprob))         
        end
    end
    singlelogl
end


function get_subpop_polymorphic(fgeno::AbstractVector, offgeno::AbstractVector, popmakeup::AbstractDict;
    minmaf::Real=0.05)
    res = []
    for popid in keys(popmakeup)
        founders = popmakeup[popid]["founder"]
        offspring = popmakeup[popid]["offspring"]
        fmono = unique(split(join(fgeno[founders]),"")) in [["1"],["2"]]
        if fmono
            push!(res, popid)
            continue
        end
        n1n2 = sum(i .> 0 for i in view(offgeno, offspring))
        p = /(n1n2...)
        if p < minmaf || p > 1 - minmaf
            push!(res, popid)
        end
    end
    setdiff(keys(popmakeup),res)
end

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
            if alleleset  in [["1"],["2"],[]]
                push!(res,p)
            end
        else
            @warn string("unxecpted offspringformat=",offspringformat)
        end
    end
    dict = Dict(magicped.founderinfo[!,:individual] .=> 1:size(magicped.founderinfo,1))
    [dict[i] for i in res]    
end

function infer_singlesite(fgeno::AbstractVector, fixedfounders::AbstractVector,
    offgeno::AbstractVector,offspringformat::AbstractString;    
    model::AbstractString,
    popmakeup::AbstractDict,
    israndallele::Bool,
    isfounderinbred::Bool, 
    isinfererror::Bool,
    byfounder::Integer,     
    likeparameters::LikeParameters,        
    priorlikeparameters::PriorLikeParameters)    
    fformat = isfounderinbred ? "GT_haplo" : "GT_unphased"        
    fhaploset,fhaploweight = priorfhaplo_singlesite(fgeno,fformat,fixedfounders; 
        genoerror=0.05,allowmissing=true,allowcorrect=true)        
    fhaplo = map((x,y)->x[rand(Categorical(y))], fhaploset, fhaploweight)    
    liketargetls, epsf, epso, pereoffspringerror, seqerror, allelebalancemean,allelebalancedisperse,alleledropout = MagicBase.extract_likeparameters(likeparameters)		    
    setdiff!(liketargetls, ["peroffspringerror"])
    if occursin(r"^GT",offspringformat)
        setdiff!(liketargetls,["seqerror","allelebalancemean","allelebalancedisperse","alleledropout"])
    end
    if model == "depmodel"
        setdiff!(liketargetls,["allelebalancemean","allelebalancedisperse","alleledropout"])
    end
    isinfererror2 = isinfererror && !isempty(liketargetls) 
    findexlist = calfindexlist(byfounder,fhaploset, popmakeup)    
    oldfhaplo = deepcopy(fhaplo)    
    pereoffspringerror = 0
    olderror = [epsf,epso, pereoffspringerror, seqerror,allelebalancemean, allelebalancedisperse, alleledropout]    
    isinferfhaplo = any(length.(fhaploset) .> 1)    
    oldlogl = -Inf
    maxit = 25
    dataprobls = MagicReconstruct.init_dataprobls_singlephase(popmakeup)
    for it in 1:maxit
        !isinferfhaplo  && !isinfererror2 && break 
        if isinferfhaplo                    
            logl = inferfhaplo_singlesite!(dataprobls,fhaplo, fhaploset, fhaploweight, findexlist, offgeno; 
                popmakeup, epsf,epso, seqerror,allelebalancemean, allelebalancedisperse, 
                alleledropout,offspringformat,priorlikeparameters,israndallele)                    
            if oldfhaplo == fhaplo 
                isinferfhaplo = false
            end
            oldfhaplo .= fhaplo
            oldlogl = logl
        end        
        if isinfererror2
            for target in liketargetls
                est,logl = infererr_singlesite!(dataprobls,target, offgeno; fhaplo, popmakeup, epsf,epso, seqerror,
                    allelebalancemean, allelebalancedisperse, alleledropout, offspringformat, 
                    priorlikeparameters,israndallele)                    
                if target == "foundererror"
                    epsf = est
                elseif target == "offspringerror"
                    epso = est
                elseif target == "seqerror"
                    seqerror = est
                elseif target == "allelebalancemean"
                    allelebalancemean = est
                elseif target == "allelebalancedisperse"
                    allelebalancedisperse = est
                elseif target == "alleledropout"
                    alleledropout = est
                else
                    @error string("unknown target=",target)
                end                 
                # println("it=", it,",target=",target,",est=", est, ", logl=",logl)
            end    
            newerror = [epsf, epso, 0, seqerror,allelebalancemean, allelebalancedisperse, alleledropout]
            diff = abs.(olderror .- newerror)            
            if (max(diff[1:3]...)<1e-3 && max(diff[4:6]...) < 1e-2) || (logl-oldlogl <= 0.1)
                isinfererror2=false
            end
            olderror .= newerror
            oldlogl = logl
        end        
        # @info string("it=", it, " errors=",round.(olderror,digits=5),
        #     ",fhaplo=",join(fhaplo), ",logl=",oldlogl, ",t=",round.([tuse1,tuse2],digits=1),"s")            
        if it == maxit 
           msg = string("reach maxit=", it, " errors=",round.(olderror,digits=5),",fhaplo=",join(fhaplo), ",logl=",oldlogl)            
           @warn msg
        end        
    end
    pereoffspringerror = 0
    esterrors = (epsf,epso, pereoffspringerror, seqerror,allelebalancemean, allelebalancedisperse, alleledropout)
    fhaplo, esterrors
end    

function inferfhaplo_singlesite!(dataprobls,fhaplo::AbstractVector, fhaploset::AbstractVector, fhaploweight::AbstractVector,
    findexlist::AbstractVector, offgeno::AbstractVector;        
    popmakeup::AbstractDict,
    epsf::Real,
    epso::Real,
    seqerror::Real,
    allelebalancemean::Real,
    allelebalancedisperse::Real,
    alleledropout::Real,
    offspringformat::AbstractString, 
    priorlikeparameters::PriorLikeParameters,
    israndallele::Bool)        
    logl = callogl_singlesite!(dataprobls,fhaplo, offgeno; popmakeup, epsf,epso, seqerror,
        allelebalancemean, allelebalancedisperse, alleledropout,offspringformat,israndallele)   
    logl += callogpri_fhaplo(fhaplo, fhaploset, fhaploweight)              
    for findex in findexlist
        fhaplols = getfhaplols(findex,fhaplo,fhaploset)
        for fhaplo2 in fhaplols                         
            logl2 = callogl_singlesite!(dataprobls,fhaplo2, offgeno; popmakeup, epsf,epso, seqerror,
                allelebalancemean, allelebalancedisperse, alleledropout,offspringformat,israndallele) 
            logl2 += callogpri_fhaplo(fhaplo2, fhaploset, fhaploweight)              
            if logl2 > logl
                logl = logl2        
                fhaplo .= fhaplo2            
            end
            # println("findex=",findex, ",logl=",logl,", fhaplo=",join(fhaplo,","))
        end   
    end 
    logpri = cal_logpri(epsf, epso, seqerror,allelebalancemean,allelebalancedisperse,
        alleledropout,priorlikeparameters)    
    logl + logpri
end

function getfhaplols(findex::AbstractVector,fhaplo::AbstractVector,fhaploset::AbstractVector)
    set1 = vec(collect(Iterators.product(fhaploset[findex]...)))
    fhaplols = [begin
        fhaplo2 = copy(fhaplo)
        fhaplo2[findex] .= seg
        fhaplo2
    end for seg in set1]
    setdiff!(fhaplols,[fhaplo])
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


function calfindexlist(byfounder::Integer, fhaploset::AbstractVector, popmakeup::AbstractDict)
    nfounder = length(fhaploset)
    if byfounder == -1 
        findexlist = [1:nfounder]        
    else
        lenls = length.(fhaploset)
        nphase = sum(lenls .>= 2) > 10 ? 2^10 : reduce(*,lenls)
        if nphase <= 64
            findexlist = [1:nfounder]            
        else
            fmissls = length.(fhaploset)
            findexlist = MagicBase.getfindexlist(byfounder,fmissls, popmakeup;defaultby=2)            
            b = [reduce(*,lenls[i]) <= 1 for i in findexlist]
            deleteat!(findexlist,b)
        end
    end
    findexlist
end


function cal_logpri(epsf::Real, 
    epso::Real,
    seqerror::Real, 
    allelebalancemean::Real,    
    allelebalancedisperse::Real,    
    alleledropout::Real,    
    priorlikeparameters::PriorLikeParameters)
    logpri = 0.0    
    prior_errls = [priorlikeparameters.foundererror, priorlikeparameters.offspringerror, priorlikeparameters.seqerror,
        priorlikeparameters.allelebalancemean, priorlikeparameters.allelebalancedisperse,priorlikeparameters.alleledropout]    
    errls = [epsf, epso, seqerror,allelebalancemean,allelebalancedisperse,alleledropout]
    for i in eachindex(errls,prior_errls)
        p =  errls[i]
        logpri += logpdf(prior_errls[i], p)        
    end    
    logpri    
end

function infererr_singlesite!(dataprobls::AbstractVector,target::AbstractString, offgeno::AbstractVector;
    fhaplo::Union{Nothing,AbstractVector}=nothing, 
    popmakeup::AbstractDict,
    popidls=keys(popmakeup),
    epsf::Real=0.005, 
    epso::Real,
    seqerror::Real, 
    allelebalancemean::Real,    
    allelebalancedisperse::Real,    
    alleledropout::Real,
    offspringformat::AbstractString,
    priorlikeparameters::PriorLikeParameters,
    israndallele::Bool)    
    accuracygoal, precisiongoal,itmax = 4, 4,50
    inputerrors = (epsf, epso, seqerror,allelebalancemean,allelebalancedisperse,alleledropout)
    function loglfun(x::Real)
        epsf, epso, seqerror,allelebalancemean,allelebalancedisperse,alleledropout = inputerrors
        if target == "offspringerror"
            epso = 1/(1+exp(-x))  # inverse of logit transformation     
        elseif target == "foundererror"
            epsf = 1/(1+exp(-x))                      
        elseif target == "seqerror"
            seqerror = 1/(1+exp(-x))                     
        elseif target == "allelebalancemean"
            allelebalancemean = 1/(1+exp(-x))    
        elseif target == "allelebalancedisperse"
            allelebalancedisperse = exp(x)
        elseif target == "alleledropout"
            alleledropout = 1/(1+exp(-x))                     
        else
            @error string("unknown genoerror type: ",target)
        end        
        logjacobi = 0 # adding logjacobi results too large estimate error rates and worse geneic map for testing Rice_F2 
        # if target in ["foundererror", "offspringerror", "seqerror","allelebalancemean"]
        #     # estimation refers to transformed space
        #     logjacobi = x - 2*log(1+exp(x)) # dp/dx = exp(x)/(1+exp(x))^2              
        # elseif target  == "alleledropout"
        #     logjacobi = 0 
        # elseif target == "allelebalancedisperse"                
        #     logjacobi = x   # dp/dx = exp(x) # estimation refers to transformed space       
        # end    
        logpri = cal_logpri(epsf, epso, seqerror,allelebalancemean,allelebalancedisperse,
            alleledropout,priorlikeparameters)
        if isnothing(fhaplo)
            logl = callogl_singlesite(offgeno; popmakeup, popidls, epso, seqerror,
                allelebalancemean, allelebalancedisperse, alleledropout,offspringformat,israndallele) 
        else
            logl = callogl_singlesite!(dataprobls,fhaplo, offgeno; popmakeup, popidls, epsf,epso, seqerror,
                allelebalancemean, allelebalancedisperse, alleledropout,offspringformat,israndallele) 
        end        
        logl + logpri + logjacobi
    end    
    if target == "offspringerror"
        xstart = log(epso/(1-epso))
    elseif target == "foundererror"
        xstart = log(epsf/(1-epsf))
    elseif target == "seqerror"
        xstart = log(seqerror/(1-seqerror))
    elseif target == "allelebalancemean"
        xstart = log(allelebalancemean/(1-allelebalancemean))
    elseif target == "allelebalancedisperse"
        xstart = max(-10,log(allelebalancedisperse))
    elseif target == "alleledropout"
        xstart = log(alleledropout/(1-alleledropout))
    else
        @error string("unknown genoerror type: ",target)
    end       
    if target == "allelebalancedisperse"
        lowbound, upbound = max(-12.0,xstart-10.0), xstart+10.0
    else
        fraction_upbound = 0.99
        lowbound, upbound = log(1e-5), log(fraction_upbound/(1-fraction_upbound))        
    end
    xstart = min(max(lowbound,xstart),upbound)
    res= MagicBase.brentMax(loglfun,lowbound,upbound;
        xstart, precisiongoal,accuracygoal,maxiter=itmax)            
    x = res[1]
    est =  target == "allelebalancedisperse" ? exp(x) : 1.0/(1.0+exp(-x)) 
    logl = res[2]
    est,logl
end

function callogl_singlesite(offgeno::AbstractVector;        
    popmakeup::AbstractDict,
    popidls=keys(popmakeup),
    epso::Real,
    seqerror::Real,
    allelebalancemean::Real,
    allelebalancedisperse::Real,
    alleledropout::Real,
    offspringformat::AbstractString,
    israndallele::Bool)     
    singlelogl = 0.0
    for popid in popidls
        offls = popmakeup[popid]["offspring"]
        ishaploid  =  popmakeup[popid]["ishaploid"]
        if ishaploid            
            if offspringformat == "AD"
                like = MagicReconstruct.haplolikeGBS(view(offgeno,offls),epso,seqerror)
            else
                like = MagicReconstruct.haplolike(view(offgeno,offls),epso)
            end            
        else
            if offspringformat == "AD"                
                like = MagicReconstruct.diplolikeGBS(view(offgeno,offls),epso,seqerror,allelebalancemean,allelebalancedisperse,alleledropout; israndallele)    
            else
                like = MagicReconstruct.diplolike(view(offgeno,offls),epso;israndallele)
            end
            
        end        
        singlelogl += sum(log.(sum(like,dims=2)[:,:]))
    end
    singlelogl
end

function infer_singlesite_rawcall(offgeno::AbstractVector,offspringformat::AbstractString;    
    popmakeup::AbstractDict, 
    model::AbstractString,        
    isinfererror::Bool,        
    likeparameters::LikeParameters,        
    priorlikeparameters::PriorLikeParameters,
    israndallele::Bool)        
    liketargetls, epsf, epso, epso_perind, seqerror, allelebalancemean,allelebalancedisperse,alleledropout = MagicBase.extract_likeparameters(likeparameters)		    
    setdiff!(liketargetls,["foundererror"])
    if occursin(r"^GT",offspringformat)
        setdiff!(liketargetls,["seqerror","allelebalancemean","allelebalancedisperse","alleledropout"])
    end
    if model == "depmodel"
        setdiff!(liketargetls,["allelebalancemean","allelebalancedisperse","alleledropout"])
    end        
    olderror = [epsf,epso, 0, seqerror,allelebalancemean, allelebalancedisperse, alleledropout]    
    isinfererror || return olderror    
    logl = oldlogl = -Inf
    maxit = 25
    dataprobls = []
    for it in 1:maxit
        for target in liketargetls
            est,logl = infererr_singlesite!(dataprobls,target, offgeno;  popmakeup, epso, seqerror,
                allelebalancemean, allelebalancedisperse, alleledropout, offspringformat, priorlikeparameters,israndallele)                    
            if target == "offspringerror"
                epso = est
            elseif target == "seqerror"
                seqerror = est
            elseif target == "allelebalancemean"
                allelebalancemean = est
            elseif target == "allelebalancedisperse"
                allelebalancedisperse = est
            elseif target == "alleledropout"
                alleledropout = est
            else
                @error string("unknown target=",target)
            end                 
            # println("it=", it,",target=",target,",est=", est, ", logl=",logl)
        end    
        newerror = [epsf, epso, 0, seqerror,allelebalancemean, allelebalancedisperse, alleledropout]
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
    peroffspringerror = 0
    esterrors = [epsf,epso, peroffspringerror, seqerror,allelebalancemean, allelebalancedisperse, alleledropout]
    esterrors
end    

function singlesite_offpostprob_raw(offgeno::AbstractVector, offspringformat::AbstractString,esterrors::AbstractVector,israndallele::Bool)
    _,epso,_, seqerror, allelebalancemean,allelebalancedisperse,alleledropout = esterrors    
    if offspringformat == "AD"
        like = MagicReconstruct.diplolikeGBS(offgeno,epso,seqerror,allelebalancemean,allelebalancedisperse,alleledropout; israndallele)    
    else
        like = MagicReconstruct.diplolike(offgeno,epso;israndallele)
    end 
    if size(like,2) == 4
        like[:,2] .+= like[:,3]        
        # uniform prior for 4 possible genotypes
        [normalize(i,1) for i in eachrow(view(like,:,[1,2,4]))]
    elseif size(like,2) == 2
        [normalize(i,1) for i in eachrow(like)]
    else
        error(string("wrong size of like: ",size(like)))
    end  
end
