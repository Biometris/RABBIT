function filter_marker(genofile::AbstractString,pedinfo::Union{Integer,AbstractString};
    model::AbstractString="jointmodel",
    epso::Real = 0.005,     
    isfounderinbred::Bool=true,
    threshcall::Real = 0.9,
    formatpriority::AbstractVector=["AD","GT"],
    isphysmap::Bool=false,
    recomrate::Real=1.0,    
    commentstring::AbstractString="##",    
    minmonotest::Integer = 20,
    mono2miss::Union{Nothing,Bool} = true,	        
    isdelinconsistent::Bool = false,    
    minmaf::Real=0.05,            
    missfilter::Function=(fmiss,omiss)-> fmiss<=1.0 && omiss<=1.0,        
    isparallel::Bool=false,
    outstem::AbstractString= "outstem",
    logfile::Union{AbstractString,IO} = outstem*"_filter_marker.log",
    workdir::AbstractString=pwd(),
    verbose::Bool=true)
    marker_mono_prob >= 1.0 && return magicgeno
    starttime = time()
    if !isfounderinbred
        @error string("TODO: filter_marker for outbred founders!")
        return
    end
    logio = MagicBase.set_logfile_begin(logfile, workdir, "filter_marker";
        verbose,delim="=")
    MagicReconstruct.info_file_arg(genofile, pedinfo, formatpriority,isphysmap, recomrate,
        commentstring,workdir, logio, verbose)
    magicgeno = formmagicgeno(genofile, pedinfo;
        isfounderinbred, formatpriority, isphysmap, recomrate,
        commentstring, missingstring="NA", workdir)
    outstem *= "_filter_marker"    
    filter_marker!(magicgeno; model, isfounderinbred, threshcall,
        minmonotest,epso,mono2miss,
        isdelinconsistent,minmaf, missfilter, 
        isparallel, outstem, logfile=logio, workdir, verbose)
    outfile = outstem*".vcf.gz"
    savegenodata(outfile,magicgeno; commentstring,workdir)
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, logio, tused,"filter_marker";
        verbose,delim="=")
    magicgeno
end

function filter_marker!(magicgeno::MagicGeno;
    model::AbstractString="jointmodel",
    epso::Real = 0.005,         
    isfounderinbred::Bool=true,    
    threshcall::Real = 0.9,
    minmonotest::Integer = 20,	    
    mono2miss::Union{Nothing,Bool} = true,	    
    isdelinconsistent::Bool = false,    
    minmaf::Real=0.05,            
    missfilter::Function=(fmiss,omiss)-> fmiss<=1.0 && omiss<=1.0,    
    isparallel::Bool=false,
    outstem::AbstractString= "outstem",
    logfile::Union{AbstractString,IO} = outstem*"_filter_marker.log",
    workdir::AbstractString=pwd(),
    verbose::Bool=true)
    starttime = time()
    marker_mono_correct = true
    logio = MagicBase.set_logfile_begin(logfile, workdir, "filter_marker!"; verbose)    
    isparallel = isparallel && nworkers()>1
    msg = string("list of options: \n",
        "model = ", model, "\n",
        "epso = ", epso, "\n",        
        "isfounderinbred = ", isfounderinbred, "\n",        
        "threshcall =", threshcall, "\n",
        "minmonotest = ", minmonotest, "\n",                
        "mono2miss = ", mono2miss, "\n",	
        "isdelinconsistent = ", isdelinconsistent, "\n",		
        "minmaf = ", minmaf, "\n",        
        "missfilter = ", first(code_lowered(missfilter)), "\n",                
        "isparallel = ", isparallel, isparallel ? string("(nworkers=",nworkers(),")") : "", "\n",        
        "workdir = ",workdir,"\n",
        "outstem = ", outstem,"\n",
        "logfile = ",logfile,"\n",
        "verbose = ",verbose)
    MagicBase.printconsole(logio,verbose,msg)    
    MagicBase.reset_juncdist!(magicgeno.magicped,model;io=logio,verbose,isfounderinbred)
    model = MagicBase.reset_model(magicgeno.magicped,model;io=logio,verbose)
    founderformat,_= MagicBase.setunphasedgeno!(magicgeno)
    MagicBase.info_magicgeno(magicgeno;io=logio,verbose)    
    MagicBase.info_missing(magicgeno;io=logio,verbose)
    if founderformat in ["GT_haplo","GT_unphased","GT_phased"]
        msg = "founders genoformat must be GT_haplo, GT_unphased, or GT_phased"
        printconsole(logio,false,"ERROR: "*msg)
        @error msg
    end    
    printconsole(logio,verbose,"-----monomorphic test for each subpopulation at each marker")
    outfiles = test_monomorphic!(magicgeno; model, epso,isfounderinbred, threshcall, 
        minmonotest,marker_mono_correct,mono2miss,
        minmaf, 
        isparallel, outstem, logio, workdir, verbose)
    MagicBase.info_missing(magicgeno;io=logio,verbose)    
    msg = string("-----filter sequentially for marker_consistent, missfilter, and maf >= ",minmaf)
    printconsole(logio,verbose,msg)
    snpsumfile = first(outfiles)    
    plot_missing_maf(snpsumfile; mono2miss, isdelinconsistent,minmaf,
        logio, workdir, verbose,outstem)    
    filter_marker!(magicgeno,snpsumfile; 
        isdelinconsistent, missfilter, minmaf, 
        isparallel, logio, workdir, verbose)
    MagicBase.info_missing(magicgeno;io=logio,verbose)
    MagicBase.info_magicgeno(magicgeno;io=logio,verbose)    
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, logio, tused,"filter_marker!"; verbose)
    outfiles
end

function filter_marker!(magicgeno::MagicGeno,snpsumfile::AbstractString;    
    isdelinconsistent::Bool = false,
    missfilter::Function=(fmiss,omiss)-> fmiss<=1.0 && omiss<=1.0,    
    minmaf::Real=0.05,                
    isparallel::Bool=false,    
    commentstring::AbstractString = "##", 
    logio::Union{Nothing,IO} = nothing,
    workdir::AbstractString=pwd(),    
    verbose::Bool=true)
    startt = time()
    nchr = length(magicgeno.markermap)    
    snpsumdf = CSV.read(getabsfile(workdir,snpsumfile),DataFrame; comment=commentstring, missingstring=["NA","missing"])
    # println(snpsumdf[1:20,:])
    snpsum_gdf = groupby(snpsumdf,:linkagegroup)    
    msg = "calculate offspring missing and MAF with mono2miss=true"
    printconsole(logio,verbose,msg)    
    if isparallel && nprocs()>1
        magicgenols =[submagicgeno(magicgeno,chrsubset=[chr]) for chr=1:nchr]
        gdfls = [snpsum_gdf[[chr]] for chr = 1:nchr]
        magicgeno.markermap = nothing
        magicgeno.foundergeno = nothing
        magicgeno.offspringgeno = nothing
        GC.gc()
        resls = pmap((x,y)-> filter_marker_chr!(x, y,1; isdelinconsistent,
            missfilter,minmaf, logio=nothing,verbose), magicgenols,gdfls)        
        if !isnothing(logio)
            for chr=1:nchr
                iobuffer = resls[chr][2]
                write(logio,String(take!(iobuffer)))
                flush(logio)
                close(iobuffer)
            end
        end
        magicgenols .= first.(resls)
        magicgeno.markermap = [i.markermap[1] for i in magicgenols]
        magicgeno.foundergeno = [i.foundergeno[1] for i in magicgenols]
        magicgeno.offspringgeno = [i.offspringgeno[1] for i in magicgenols]
        magicgenols =  nothing	
        ndells = sum(last.(resls))
        GC.gc()
    else
        ndells = zeros(Int,3)
        for chr in 1:nchr
            ndells .+= last(filter_marker_chr!(magicgeno,snpsum_gdf, chr; isdelinconsistent,
                missfilter, minmaf, logio,verbose))
        end
    end
    msg = string("all_chrs, #inconsistent=",ndells[1],
        ", #missing=",ndells[2], ", #(maf<",minmaf,")=",ndells[3],
        ", #del=",sum(ndells)," of ", size(snpsumdf,1))
    msg *= string(", t=",round(time()-startt,digits=1),"s")
    MagicBase.printconsole(logio,verbose,msg)
end

function filter_marker_chr!(magicgeno::MagicGeno,snpsum_gdf::GroupedDataFrame, chr::Integer;   
    isdelinconsistent::Bool = false,
    minmaf::Real=0.05,            
    missfilter::Function=(fmiss,omiss)-> fmiss<=1.0 && omiss<=1.0,    
    logio::Union{Nothing,IO} = nothing,
    verbose::Bool=true)
    startt = time()
    isnothing(logio) && (logio=IOBuffer(append=true))
    check_chr_marker(magicgeno,snpsum_gdf,chr;logio,verbose)       
    chrid = string(magicgeno.markermap[chr][1,:linkagegroup])
    nsnp = size(magicgeno.markermap[chr],1)
    msg = string("chr=",chrid)
    sumdf = snpsum_gdf[chr]
    # filter snp with inconsistent monomorphic populations
    if isdelinconsistent
        inconsistls = sumdf[!,:marker_inconsistent]
        ndel_inconsist = sum(inconsistls)
    else
        inconsistls = falses(nsnp)
        ndel_inconsist = 0
    end
    msg *= string(", #inconsistent=",ndel_inconsist)
    # filter missing    
    fmissls = sumdf[!,:founder_miss]
    offmissls = sumdf[!,:offspring_miss_mono2miss]            
    missdells = map((f,o)->!missfilter(f,o), fmissls, offmissls)
    ndel_miss = sum(missdells .&& .!inconsistls)
    msg *= string(", #missing=",ndel_miss)
    # filter MAF
    maf = [isnan(i) ? i : (i>0.5 ? 1.0-i : i) for i in sumdf[!,:freqA2_mono2miss]]
    mafdells = [isnan(i) ? false : i < minmaf for i in maf]
    ndel_maf = sum(mafdells .&& .!missdells .&& .!inconsistls)
    msg *= string(", #(maf<",minmaf, ")=",ndel_maf)
    keepls = .!(mafdells .|| missdells .|| inconsistls)
    if !all(keepls)
        magicgeno.foundergeno[chr] = magicgeno.foundergeno[chr][keepls,:]
        magicgeno.offspringgeno[chr] = magicgeno.offspringgeno[chr][keepls,:]
        magicgeno.markermap[chr] = magicgeno.markermap[chr][keepls,:]
    end
    ndells = [ndel_inconsist,ndel_miss,ndel_maf]
    ndel = nsnp - sum(keepls)
    ndel == sum(ndells) || @warn string("inconsistent #del_markers")
    msg *= string(", #del=",ndel, " of ",nsnp)
    msg *= string(", t=",round(time()-startt,digits=1),"s")
    MagicBase.printconsole(logio,verbose,msg)
    magicgeno,logio, ndells
end

function check_chr_marker(magicgeno::MagicGeno,snpsum_gdf::GroupedDataFrame,chr::Integer;
    logio::Union{Nothing,IO}=nothing,verbose::Bool=true)    
    nchr = length(magicgeno.markermap)
    if nchr != length(snpsum_gdf) 
        msg = "inconsistent number of chromosomes"
        printconsole(logio,false,msg)
        verbose && @error msg
    end
    res = true
    chrid  = string(magicgeno.markermap[chr][1,:linkagegroup])
    chrid2 = string(snpsum_gdf[chr][1,:linkagegroup]) 
    if chrid != chrid2
        msg = string("inconsistent [chrid,chrid2] =", [chrid,chrid2])
        printconsole(logio,false,msg)
        verbose && @error msg
        res = false
    end
    markers = string.(magicgeno.markermap[chr][!,:marker]); 
    markers2 = string.(snpsum_gdf[chr][!,:marker])
    if issetequal(markers,markers2)
        b = markers .!= markers2
        if any(b)
            msg = string("inconsistent marker ordering in chr=",chrid)   
            msg *= string(",magicgeno marekrs: ",markers[b], ", snpsum markers: ",markers2[b])
            printconsole(logio,false,msg)
            verbose && @error msg
            res = false
        end
    else
        msg = string("inconsistent markers in chr=",chrid)     
        d12 = setdiff(markers,markers2)
        isempty(d12) || (msg *= string(", markers in magicgeno but not in snpsum: ",join(d12,",")))
        d21 = setdiff(markers2,markers)
        isempty(d21) || (msg *= string(", markers in snpsum but not in magicgeno: ",join(d21,",")))
        printconsole(logio,false,msg)
        verbose && @error msg
        res = false
    end
    res
end

function test_monomorphic!(magicgeno::MagicGeno;
    model::AbstractString="jointmodel",
    epso::Real = 0.005,     
    isfounderinbred::Bool=true,    
    threshcall::Real = 0.9,
    minmonotest::Integer = 20,	        
    marker_mono_correct::Bool = true,
    mono2miss::Union{Nothing,Bool} = true,	    
    minmaf::Real=0.05,            
    isparallel::Bool=false,
    outstem::AbstractString= "outstem",
    logio::Union{Nothing,IO} = nothing,
    workdir::AbstractString=pwd(),
    verbose::Bool=true)                
    magicprior= MagicReconstruct.calmagicprior(magicgeno,model;isfounderinbred)    
    popmakeup = MagicReconstruct.calpopmakeup(magicgeno,1, model, magicprior;isfounderinbred)                    
    marker_mono_prob = 0.99
    info_freqa2_minsize(popmakeup;marker_mono_prob,epso, minmonotest,io=logio,verbose)
    nchr = length(magicgeno.markermap)    
    outfiles = [getabsfile(workdir,outstem*i) for i in ["_marker_summary.csv", "_founder_change.csv"]]  
    open(outfiles[1],"w") do outio
        sumcols = ["linkagegroup", "marker"]
        popidls = collect(keys(popmakeup)) 
        append!(sumcols, [string(popid,"_A1|A2") for popid in popidls])
        append!(sumcols, ["marker_inconsistent", "founder_miss", "offspring_miss", 
            "offspring_miss_mono2excl", "offspring_miss_mono2miss", 
            "freqA2", "freqA2_mono2miss","freqnonmiss_mono2miss"])
        msg = "##col_1, linkagegroup, linkage group ID\n"
        msg *= "##col_2, marker, marker ID\n"
        for i in eachindex(popidls)
            msg *= string("##col_",i+2, ", ", string(popidls[i],"_A1|A2"), ", n1|n2 where n1 = number of 1st allele A1 && n2 = number of 2nd allele A2 for subpoulation = ",popidls[i], " && \"=>0|0\" denotes the genotypes in the subpopulation being set to missing\n")
        end
        descripls = [
            "marker_inconsistent, if the marker has inconsistent founder genotypes based on the dominant alleles in each monomorphic subpoulation", 
            "founder_miss, fraction of missing founder genotypes after filtering && it is used by missfilter", 
            "offspring_miss, fraction of missing offspring genotypes after filtering && offspring_miss is the same as offspring_miss_mono2miss if mono2miss = true", 
            "offspring_miss_mono2excl, fraction of missing offspring genotypes after excluding monomorphic subpoulations", 
            "offspring_miss_mono2miss, fraction of missing offspring genotypes after setting genotypes in monomorphic subpoulations to missing && it is used by missfilter", 
            "freqA2, frequency of 2nd allele A2 after filtering && freqA2 is the same as freqA2_mono2miss if mono2miss = true", 
            "freqA2_mono2miss, frequency of 2nd allele A2 after setting genotypes in monomorphic subpoulations to missing",             
            "freqnonmiss_mono2miss, frequency of nonmissing offspring genotypes (excluding halfcalled 1N and 2N) after setting genotypes in monomorphic subpoulations to missing"
        ]
        for i in eachindex(descripls)
            msg *= string("##col_", i+length(popidls)+2,",", descripls[i],"\n")
        end
        write(outio, msg)
        write(outio, join(sumcols,","), "\n")                                        
    end
    colnames = Symbol.(["linkagegroup","marker","individual","oldgeno","newgeno",
        "subpop","n1","n2","chr_index","marker_index","ind_index","marker_inconsistent"])    
    open(outfiles[2],"w") do outio
        descripls = ["linkagegroup, linkage group ID",
            "marker, marker ID",
            "individual, founder ID",
            "oldgeno, founder genotype before filtering && allele code N = missing && allele code 1 = reference allele && allele code 2 = alternative allele",
            "newgeno, founder genotype after filtering",
            "subpop, subpopulation ID",
            "n1, number of 1st allele (refernece alelle)",
            "n2, number of 2nd allele (alternative alelle)",
            "chr_index, index of the chromosome",
            "marker_index, index of the marker",
            "ind_index, index of the founder",
            "marker_inconsistent, if the marker has inconsistent founder genotypes based on the dominant alleles in each monomorphic subpoulation", 
        ]
        msg = ""
        for i in eachindex(descripls)
            msg *= string("##col_", i, ", ", descripls[i], "\n")
        end
        write(outio, msg)
        write(outio,join(colnames,","),"\n")                
    end    
    msg = string("#marker_inconsist = #markers with inconsistent monomorphic subpopulations\n")
    msg *= string("#impute_f = #imputed founder genotypes out of missing founder genotypes\n")
    marker_mono_correct && (msg *= string("#correct_f = #corrected founder genotypes out of observed founder genotypes\n"))    
    printconsole(logio,verbose,msg)  
    startt = time()
    if isparallel && nprocs()>1
        magicgenols =[submagicgeno(magicgeno,chrsubset=[chr]) for chr=1:nchr]
        magicgeno.markermap = nothing
        magicgeno.foundergeno = nothing
        magicgeno.offspringgeno = nothing
        GC.gc()
        resls = pmap(x-> test_monomorphic_chr!(x, 1, model,magicprior; 
            epso,isfounderinbred, threshcall, minmaf, marker_mono_prob, 
            minmonotest,marker_mono_correct,mono2miss,
            logio=nothing,outstem,workdir,verbose), magicgenols)
        magicgenols .= first.(resls)
        magicgeno.markermap = [i.markermap[1] for i in magicgenols]
        magicgeno.foundergeno = [i.foundergeno[1] for i in magicgenols]
        magicgeno.offspringgeno = [i.offspringgeno[1] for i in magicgenols]
        changecounts = [i[2] for i in resls]
        if !isnothing(logio)
            for chr=1:nchr
                iobuffer = resls[chr][3]
                write(logio,String(take!(iobuffer)))
                flush(logio)
                close(iobuffer)
            end
        end
        chroutfilesls= [i[4] for i in resls]
        for chroutfiles in chroutfilesls
            for i in eachindex(chroutfiles)
                open(outfiles[i],"a") do outio
                    if isfile(chroutfiles[i])
                        write(outio,read(chroutfiles[i]))
                        flush(outio)
                    end
                end
            end 
            rm.(chroutfiles;force=true)
        end 
    else
        changecounts = Vector(undef, nchr)
        for chr in 1:nchr
            _, changecounts[chr],_, chroutfiles = test_monomorphic_chr!(magicgeno, chr, model,magicprior; 
                epso,isfounderinbred, threshcall, minmaf, marker_mono_prob,
                minmonotest,marker_mono_correct,mono2miss,
                logio,outstem,workdir,verbose
            )
            for i in eachindex(chroutfiles)
                open(outfiles[i],"a") do outio
                    if isfile(chroutfiles[i])
                        write(outio,read(chroutfiles[i]))
                        flush(outio)
                    end
                end
            end  
            rm.(chroutfiles;force=true)
        end
    end
    mono_inconsist, f_impute, f_correct, off_hetero_2missing = sum(changecounts)
    off_hetero = round(/(off_hetero_2missing...),digits=5)
    if off_hetero > 0
        msg = "heterozygous offspring genotypes for subpopulations with depmodel are set to missing. \n"
        msg *= string("offspring frac_heterozygous out of observed genotypes = ",off_hetero)
        MagicBase.printconsole(logio,verbose,msg)
    end
    f_change = [round(/(i...),digits=5) for i in [mono_inconsist,f_impute, f_correct]]
    msg = string("all_chrs, marker_inconsistent=", f_change[1],
        ", impute_founder=", f_change[2])
	if marker_mono_correct
		msg *= string(", correct_founder=", f_change[3])
	end
    msg *= string(", t=",round(time()-startt,digits=1),"s") 
    MagicBase.printconsole(logio,verbose,msg)           
    outfiles
end

function test_monomorphic_chr!(magicgeno::MagicGeno, chr::Integer,
    model::AbstractString, magicprior::NamedTuple;
    epso::Real = 0.005,
    isfounderinbred::Bool=true,
    threshcall::Real = 0.9,
    minmaf::Real=0.05,            
    marker_mono_prob::Real = 0.99,
    minmonotest::Integer = 20,	        
    marker_mono_correct::Bool = true,
    mono2miss::Union{Nothing,Bool} = true,	    
    logio::Union{Nothing,IO}=nothing,
    outstem::AbstractString= "outstem",
    workdir::AbstractString=pwd(),
    verbose::Bool=true)    
    startt = time()
    isnothing(logio) && (logio=IOBuffer(append=true))
    chrid = magicgeno.markermap[chr][1,:linkagegroup]
    chrid = ismissing(chrid) ? "NA" : chrid
    colnames = Symbol.(["linkagegroup","marker","individual","oldgeno","newgeno",
        "subpop","n1","n2","chr_index","marker_index","ind_index","marker_inconsistent"])
    outfiles = [getabsfile(workdir,outstem*"_"*chrid*i) for i in ["_marker_summary.csv", "_founder_change.csv"]]  
    # colnames of outfile1: chromosome,marker,each_subpopid, marker_inconsistent,founder_miss,offspring_miss,offspring_miss_mono2excl,offspring_miss_mono2miss,freqA2,freqA2_mono2miss
    # colnames of outfile2: colnames    
    popmakeup = MagicReconstruct.calpopmakeup(magicgeno,chr, model, magicprior;isfounderinbred)
    off_hetero_2missing = hetero2missing!(magicgeno,chr,popmakeup) # for pop with ishaploid = true    
    founderid = magicgeno.magicped.founderinfo[!,:individual]    
    foundergeno = magicgeno.foundergeno[chr]
    offgeno = magicgeno.offspringgeno[chr]
    offformatls = magicgeno.markermap[chr][!,:offspringformat]
    pop2offindices,offspringidls = get_pop2offindices(magicgeno.magicped)  
    # calledgeno    
    calledgeno = call_offgeno(magicgeno,model,chr; callthreshold=threshcall)	    
    popinfols = get_popinfols(popmakeup; epso, marker_mono_prob, minmonotest)   
    # cal #obs and #missing    
    nf_miss = sum(foundergeno .== "N") + sum(foundergeno .== "NN")
    nf_obs = length(foundergeno) - nf_miss    
    nsnp = size(calledgeno,1)             
    nf_impute = 0
    nf_correct= 0    
    nmarker_inconsist = 0    
    open(outfiles[1],"w") do io_allelefreq
        open(outfiles[2],"w") do io_fchange            
            popsizels = [length(val["offspring"]) for val in values(popmakeup)]    
            foundersls = [val["founder"] for val in values(popmakeup)]    
            totfounders = unique(reduce(vcat,foundersls))
            totpopsize = sum(popsizels)          
            # println("totsize = ", totpopsize,", popsizels=", popsizels)
            founderchanges = Vector()                
            for snp in 1:nsnp
                nnonmiss_mono2miss = 0
                marker_inconsistent = false
                snpid = magicgeno.markermap[chr][snp,:marker]
                founderchanges = Vector()                    
                n1n2dict = OrderedDict{String,Vector{Vector{Int}}}()
                for popinfo in popinfols                    
                    popid, ishaploid, founders, offspring, freq_a2, min_popsize = popinfo                    
                    offseq = calledgeno[snp,offspring]
                    if ishaploid
                        n1 = 2*sum(offseq .== "11")
                        n2 = 2*sum(offseq .== "22")
                        n11 = n1; 
                        n22 = n2
                        n12 = 0 
                    else
                        n11 = sum(offseq .== "11")
                        n12 = sum(offseq .== "12")
                        n22 = sum(offseq .== "22")
                        n1 = 2*n11+n12+sum(offseq .== "1N")
                        n2 = 2*n22+n12+sum(offseq .== "2N")                        
                    end            
                    push!(n1n2dict, popid=>[[n1,n2]])                    
                    nnonmiss_mono2miss += n11 + n12 + n22
                    if length(popinfols)==1 || (1-minmaf >= n1/(n1+n2) >= minmaf)
                        continue
                    else    
                        if n1+n2 < 2*min_popsize                        
                            if !isnothing(mono2miss) && mono2miss 
                                misscode = MagicBase.get_missingcode(offformatls[snp])   
                                for i in offspring
                                    offgeno[snp,i] = misscode
                                end
                                if 0 in [n1,n2]                                
                                    n1n2dict[popid] = [[0,0]]
                                else
                                    n1n2dict[popid] = [[n1,n2],[0,0]]
                                end                      
                                nnonmiss_mono2miss -= n11 + n12 + n22 # removing nonmiss from monomorphic subpopulations
                            end                            
                            continue
                        end
                    end
                    prob_post = cal_post_prob(n1,n2,freq_a2)                        
                    ischange = true
                    if last(prob_post) > marker_mono_prob                        
                        newgeno = isfounderinbred ? "2" : "22"
                    elseif first(prob_post) > marker_mono_prob                        
                        newgeno = isfounderinbred ? "1" : "11"
                    else
                        ischange = false
                    end        
                    if ischange                            
                        # TODO: for outbred founders
                        # booking monomorphic
                        fgeno_snp = foundergeno[snp,founders]    
                        for a in unique(fgeno_snp)
                            for i in founders[fgeno_snp .== a]
                                # false denotes default marker_inconsistent
                                push!(founderchanges, [chrid, snpid, founderid[i],
                                    a, newgeno,popid, n1,n2,chr,snp,i,false]) 
                            end
                        end       
                        nnonmiss_mono2miss -= n11 + n12 + n22 # removing nonmiss from monomorphic subpopulations
                    end
                end                                
                if !isempty(founderchanges)  
                    founderchanges2 = permutedims(reduce(hcat,founderchanges))
                    fchangedf = DataFrame(founderchanges2,colnames)                    
                    MagicFilter.set_mono_inconsistent!(fchangedf)                    
                    inconsistls = Vector{Bool}(fchangedf[!,:marker_inconsistent])  
                    marker_inconsistent = any(inconsistls) 
                    # set founder genotypes to missing if they have marker_inconsistent imputations or corrections                            
                    if marker_inconsistent 
                        nmarker_inconsist += 1 
                        fchangedf_inconsist = fchangedf[inconsistls,:]
                        inconsistent_founder_2miss!(foundergeno,fchangedf_inconsist)                            
                    end
                    # founder imputation or correction
                    fchangedf_consist = fchangedf[.!inconsistls,:]
                    if !isempty(fchangedf_consist)
                        nimpute, ncorrect = change_founder!(foundergeno,fchangedf_consist,isfounderinbred,marker_mono_correct)                        
                        nf_impute += nimpute
                        nf_correct += ncorrect
                    end                    
                    # change offspring genotypes in monomorphic subpopulations
                    # offchanges[i] = [chrid, snpid, subpop]
                    if !isnothing(mono2miss)
                        offchanges = change_offspring!(offgeno,calledgeno,fchangedf,offformatls,pop2offindices,mono2miss)                                                                         
                        if mono2miss                                        
                            for  changes in offchanges
                                pop = changes[3] 
                                if 0 in first(n1n2dict[pop])     
                                    # not book change without corrections
                                    n1n2dict[pop] = [[0,0]]
                                else                                
                                    push!(n1n2dict[pop], [0,0])                                
                                end
                            end                        
                        else
                            for changes in offchanges
                                pop = changes[3] 
                                oldn1,oldn2 = first(n1n2dict[pop])
                                if oldn1*oldn2 > 0
                                    # approximate for genotype = AD, where reads of wrong allele were removed
                                    oldn1 > oldn2 && push!(n1n2dict[pop], [oldn1,0])
                                    oldn2 > oldn1 && push!(n1n2dict[pop], [0,oldn2])
                                end
                            end
                        end
                    end
                    b = fchangedf[!,:oldgeno] .== fchangedf[!,:newgeno] .&& .!fchangedf[!,:marker_inconsistent]
                    deleteat!(fchangedf,b)
                    CSV.write(io_fchange,fchangedf;append=true,missingstring="NA")
                end                    
                # no excl
                fgenosnp = foundergeno[snp,totfounders]
                f_miss = round((sum(fgenosnp .== "N") + sum(fgenosnp .== "NN"))/length(fgenosnp),digits=4)
                n1n2ls = last.(values(n1n2dict))
                n1n2 = sum(n1n2ls)
                off_miss = round(1.0 - sum(n1n2)/(2*totpopsize),digits=4)
                freqA2 = round(n1n2[2]/sum(n1n2),digits=4)
                # mono2excl or mono2miss                    
                b = [!in(0,i) for i in n1n2ls] # b is a list of indicators of polymorphic subpoplations
                if any(b)                    
                    n1n2 = sum(n1n2ls[b])
                    off_miss_mono2excl = round(1.0 - sum(n1n2)/(2*sum(popsizels[b])),digits=4)
                    off_miss_mono2miss = round(1.0 - sum(n1n2)/(2*sum(popsizels)),digits=4)
                    freqA2_mono2miss = round(n1n2[2]/sum(n1n2),digits=4)  # same as freqA2_mono2miss
                else
                    off_miss_mono2excl = NaN
                    off_miss_mono2miss = 1.0
                    freqA2_mono2miss = NaN
                end
                n1n2msg = [join(join.(i,"|"),"=>") for i in values(n1n2dict)]
                freqnonmiss_mono2miss = round(nnonmiss_mono2miss/totpopsize,digits=4)
                snpsumls = (marker_inconsistent,f_miss,off_miss,off_miss_mono2excl,off_miss_mono2miss,freqA2,freqA2_mono2miss,freqnonmiss_mono2miss)
                freqmsg = string(join(n1n2msg,","), ",", join(snpsumls,","))
                write(io_allelefreq,join([chrid, snpid],","), ",", freqmsg, "\n")                     
            end
        end        
    end    
    changecounts = [[nmarker_inconsist,nsnp],[nf_impute,nf_miss],[nf_correct,nf_obs], off_hetero_2missing]
    msg = string("chr=", chrid,", #marker_inconsist=", join(changecounts[1],"|"))
    msg *= string(", #impute_f=", join(changecounts[2],"|"))
	marker_mono_correct && (msg *= string(", #correct_f=", join(changecounts[3],"|")))	
    msg *= string(", t=",round(time()-startt,digits=1),"s")
    MagicBase.printconsole(logio,verbose,msg)    
    magicgeno, changecounts,logio,outfiles
end

function hetero2missing!(magicgeno::MagicGeno, chr::Integer,
    popmakeup::AbstractDict)
    offgeno = magicgeno.offspringgeno[chr]
    formatls = magicgeno.markermap[chr][!,:offspringformat]
    nobs = length(offgeno) - sum(MagicBase.count_missing(offgeno,formatls))
    nmiss = 0
    for (popid, makeup) in popmakeup
        offspring = popmakeup[popid]["offspring"]
        ishaploid = popmakeup[popid]["ishaploid"]
        if ishaploid
            b = formatls .== "GT_unphased"
            if any(b)
                subgeno = view(offgeno,b,offspring)
                for i in eachindex(subgeno)
                    if subgeno[i] in ["12","21"]
                        subgeno[i] = "NN"
                        nmiss += 1
                    end
                end
            end
            b = formatls .== "AD"
            if any(b)
                subgeno = view(offgeno,b,offspring)
                for i in eachindex(subgeno)
                    if prod(subgeno[i]) > 0
                        subgeno[i] .= [0,0]
                        nmiss += 1
                    end
                end
            end
            b = formatls .== "GP"
            if any(b)
                subgeno = view(offgeno,b,offspring)
                for i in eachindex(subgeno)
                    g = subgeno[i]
                    if length(g) == 3
                        g[2] = 0.0
                        s = sum(g)
                        if s ≈ 0
                            subgeno[i] .= [0.5, 0.0, 0.5]
                            nmiss += 1
                        else
                            subgeno[i] .= g ./ s
                        end
                    elseif length(g) == 4
                        g[[2,3]] .= 0.0
                        s = sum(g)
                        if s ≈ 0
                            subgeno[i] .= [0.5, 0.0, 0.0, 0.5]
                            nmiss += 1
                        else
                            subgeno[i] .= g ./ s
                        end
                    end
                end
            end
        end
    end
    [nmiss,nobs]
end

function get_popinfols(popmakeup::AbstractDict; 
    epso::Real,marker_mono_prob::Real,minmonotest::Integer)    
    res = [begin 
        isfounderinbred = popmakeup[popid]["isfounderinbred"]
        founders = popmakeup[popid]["founder"]
        fgls = isfounderinbred ? founders : reduce(vcat, [[2*p-1,2*p] for p in founders])
        offspring = popmakeup[popid]["offspring"]        
        initprob = popmakeup[popid]["initprob"]
        nzorigin = popmakeup[popid]["nzorigin"]
        ishaploid = popmakeup[popid]["ishaploid"]
        freq_a2 = possible_freq_a2(ishaploid,nzorigin, fgls,initprob,epso)
        min_popsize = cal_min_popsize(freq_a2, marker_mono_prob; minmonotest)    
        (popid, ishaploid, founders, offspring, freq_a2, min_popsize)
    end for popid in keys(popmakeup)]
    res
end

# update inconsist status among  monomorphic subpos at the marker.
# assume fchangedf resulting from a single marker
function set_mono_inconsistent!(fchangedf::DataFrame)
    gdf  = groupby(fchangedf,:ind_index)
    len = [length(unique(i[!,:newgeno])) for i in gdf]
    b = len .> 1    
    if any(b)        
        for df in  gdf[b]
            df[!,:marker_inconsistent] .= true            
        end
    end    
    fchangedf
end

function inconsistent_founder_2miss!(foundergeno::AbstractMatrix,fchangedf_inconsistent::DataFrame)
    gchangedf = groupby(fchangedf_inconsistent,[:chr_index,:marker_index,:ind_index])    
    for df in gchangedf
        snp = df[1,:marker_index]
        ind = df[1,:ind_index]
        oldg = df[1,:oldgeno]
        foundergeno[snp,ind] = length(oldg)==1 ? "N" : "NN"        
    end	
    return nothing
end

function change_founder!(foundergeno::AbstractMatrix,fchangedf::DataFrame,    
    isfounderinbred::Bool,
	marker_mono_correct::Bool)
    nimpute = 0
    ncorrect = 0    
    gchangedf = groupby(fchangedf,[:chr_index,:marker_index,:ind_index])    
    for df in gchangedf
        snp = df[1,:marker_index]
        ind = df[1,:ind_index]
        oldg = df[1,:oldgeno]        
        allequal(df[!,:newgeno]) || @warn string("unequal changes of fouder genotypess: ",unique(df[!,:newgeno]))
        newg = df[1,:newgeno] # newg in ["11","22"] or ["1","2"] (inbredfounder)
        oldg == newg && continue
        if isfounderinbred 
            if oldg in ["N"]
                # impute inbred founder
                nimpute += 1
                foundergeno[snp,ind] = newg        
            else
                oldg in ["1","2"] || @error string("unknown founder genotype: ",oldg)
                # correct inbred founder
                if marker_mono_correct && oldg != newg        
                    ncorrect += 1
                    foundergeno[snp,ind] = newg
                end
            end
        else
            if oldg in ["NN"]
                # impute inbred founder
                nimpute += 1
                foundergeno[snp,ind] = newg        
            elseif oldg in ["1N","N1"]
                if newg == "22"
                    if marker_mono_correct  
                        ncorrect += 1
                        foundergeno[snp,ind] = newg
                    end
                else
                    newg == "11"  || @error string("unknown new founder genotype: ",newg)
                    nimpute += 1
                    foundergeno[snp,ind] = newg  
                end
            elseif oldg in ["2N","N2"]
                if newg == "11"
                    if marker_mono_correct  
                        ncorrect += 1
                        foundergeno[snp,ind] = newg
                    end
                else
                    newg == "22"  || @error string("unknown new founder genotype: ",newg)
                    nimpute += 1
                    foundergeno[snp,ind] = newg  
                end
            elseif oldg in ["11","12","21", "22"]                
                # correct inbred founder
                if marker_mono_correct && oldg != newg        
                    ncorrect += 1
                    foundergeno[snp,ind] = newg
                end
            else
                @error string("unknown founder genotype: ",oldg)
            end
        end
    end	
    nimpute, ncorrect
end

function info_freqa2_minsize(popmakeup::AbstractDict;
    marker_mono_prob::Real,epso::Real,
    minmonotest::Integer,
    io::Union{IO,Nothing}, verbose::Bool)
    pop_freqa2_minsize = Dict([begin    
        isfounderinbred = popmakeup[popid]["isfounderinbred"]
        founders = popmakeup[popid]["founder"]
        fgls = isfounderinbred ? founders : reduce(vcat, [[2*p-1,2*p] for p in founders])
        offspring = popmakeup[popid]["offspring"]        
        initprob = popmakeup[popid]["initprob"]
        nzorigin = popmakeup[popid]["nzorigin"]
        ishaploid = popmakeup[popid]["ishaploid"]
        freq_a2 = possible_freq_a2(ishaploid,nzorigin, fgls,initprob,epso)        
        min_popsize = cal_min_popsize(freq_a2, marker_mono_prob; minmonotest)
        popid=>(freq_a2,min_popsize)
    end for popid in keys(popmakeup)])
    popidls = collect(keys(pop_freqa2_minsize))
    valuels = values(pop_freqa2_minsize)
    for val in unique(valuels)
        b = [val == i  for i in valuels]
        msg = string("sub-populations ", popidls[b], ": ",
            "  possible a2-freq = ", round.(val[1],digits=4),
            "; minimum sub-popsize = ", val[2], " for monomorphic test")
        MagicBase.printconsole(io,verbose,msg)
    end
end

function possible_freq_a2(ishaploid::Bool, nzorigin::AbstractVector,
    fgls::AbstractVector,initprob::AbstractVector,epso::Real)
    nf = length(fgls)
    if allequal(initprob) || nf > 20
        freq_a2 = [i/nf for i in 0:nf]
    else
        if ishaploid
            hh = Iterators.product(repeat([0:1],nf)...)
            freq_a2 = unique([dot(i, initprob) for i in hh])
        else
            dict= Dict(fgls .=> 1:nf)
            nzorigin2 = [[dict[i] for i in j] for j in nzorigin]
            hh = Iterators.product(repeat([0:1],nf)...)
            freq_a2 = unique([0.5*dot([sum(h[j]) for j in nzorigin2], initprob) for h in hh])
        end
    end        
    freq_a2 = [i*(1-epso)+(1-i)*epso for i in freq_a2]    
    sort(freq_a2)
end


function cal_post_prob(n1::Integer, n2::Integer,freq_a2::AbstractVector)
    # logl = [begin
    #     if f_a2 <= 0.5
    #         alpha = 1.0
    #         beta = alpha/f_a2 - alpha
    #     else
    #         beta = 1.0
    #         alpha = beta/(1-f_a2) - beta
    #     end
    #     logbeta(n2+alpha,n1+beta) - logbeta(alpha,beta)
    # end for f_a2  in freq_a2]
    logl = [n1*log(1-i)+n2*log(i) for i in freq_a2]
    prob_post = normalize(exp.(logl .- max(logl...)),1)
    prob_post
end

function cal_min_popsize(freq_a2::AbstractVector,marker_mono_prob::Real;
    minmonotest::Integer=20)
    function f(n2,freq_a2, marker_mono_prob)
        prob_post = cal_post_prob(n1,n2,freq_a2)
        d = last(prob_post) - marker_mono_prob
        # println("n2=", n2,";d=", d, "; prob_post=",prob_post)
        d>0 ? -d : d
    end
    n1 = 0
    n2 = minmonotest - 1
    fn = f(n2,freq_a2, marker_mono_prob)
    step = 1
    while n2 <= 1000
        n2 +=step
        fn_new = f(n2,freq_a2, marker_mono_prob)
        if fn_new>fn
            fn = fn_new
        else
            break
        end
    end
    max(minmonotest,n1+n2)
end

function get_pop2offindices(magicped::MagicPed)
    offls = magicped.offspringinfo[!,:individual]
    offdict = Dict(offls .=> 1:length(offls))
    pop2offs = MagicBase.get_subpop2offspring(magicped)
    pop2offindices = Dict([pop => [get(offdict,i,nothing) for i in offs] for (pop,offs) in pop2offs])
    pop2offindices,offls
end

function change_offspring!(offgeno::AbstractMatrix,    
    calledgeno::AbstractMatrix,
    fchangedf::AbstractDataFrame,
    offformatls::AbstractVector,
    pop2offindices::AbstractDict,
    mono2miss::Bool)
    offchanges = []
    gdf = groupby(fchangedf,[:chr_index,:marker_index,:subpop])
    for df in gdf
        chr = df[1,:chr_index]
        snp = df[1,:marker_index]
        subpop = df[1,:subpop]
        offspring = pop2offindices[subpop]
        gformat = offformatls[snp]        
        misscode = MagicBase.get_missingcode(gformat)   
        if mono2miss
            # set all offspring genotypes to missing                           
            for i in offspring
                offgeno[snp,i] = misscode
            end
            push!(offchanges,[chr,snp,subpop])
        else
            # correct inconsistent offspring genotypes
            n1 = df[1,:n1]
            n2 = df[1,:n2]
            if n1*n2 > 0
                wrongallele = n1 > n2 ? "2" : "1"                
                b = [occursin(wrongallele,i) for i in calledgeno[snp,offspring]]
                offspring2 = offspring[b]
                if occursin("GT", gformat)
                    offgeno[snp,offspring2] .= replace.(offgeno[snp,offspring2],wrongallele => "N")
                elseif occursin("AD", gformat)
                    # AD genotype is a vector [r1,r2], the number of reads for each allele
                    wrongindex = wrongallele == "1" ? 1 : 2
                    for i in offspring2
                        offgeno[snp,i][wrongindex] = 0
                    end
                else
                    for i in offspring2
                        offgeno[snp,i] = misscode
                    end
                end
                push!(offchanges,[chr,snp,subpop])
            end    
        end    
    end
    offchanges
end

function plot_missing_maf(snpsumfile::AbstractString;
    mono2miss::Union{Nothing,Bool}=true,
    isdelinconsistent::Bool = false,
    minmaf::Real=0.05,
    commentstring::AbstractString = "##", 
    outstem::AbstractString= "outstem",    
    workdir::AbstractString=pwd(),
    logio::Union{Nothing,IO} = nothing,
    verbose::Bool=true)    
    snpsumdf = CSV.read(getabsfile(workdir,snpsumfile),DataFrame; comment=commentstring, missingstring=["NA","missing"])
    try 
        MagicFilter.plot_missdist(snpsumdf; mono2miss,isdelinconsistent)
        outfile = outstem*"_marker_missing.png"
        savefig(getabsfile(workdir,outfile))   
        msg = string("plot missing distribution in ", outfile)  
        printconsole(logio,verbose,msg)
    catch err
        msg = string(err,". Could not plot marker missing distribution")
        verbose && @warn msg
        printconsole(logio,false,"Warning: "*msg)
    end
    try 
        MagicFilter.plot_mafdist(snpsumdf; mono2miss, isdelinconsistent, annotate_MAF = minmaf)
        outfile = outstem*"_marker_MAF.png"
        savefig(getabsfile(workdir,outfile))    
        msg = string("plot maf distribution in ", outfile)  
        printconsole(logio,verbose,msg) 
    catch err
        msg = string(err, ". Could not plot marker maf distribution")
        verbose && @warn msg
        printconsole(logio,false,"Warning: "*msg)
    end
    return 
end


function plot_missdist(snpsumdf::AbstractDataFrame;        
    mono2miss::Union{Nothing,Bool}=true,
    isdelinconsistent::Bool = false,
    annotate_foundermiss::Union{Nothing,Real}=nothing,
    annotate_offspringmiss::Union{Nothing,Real}=nothing)
    col = :founder_miss
    if isdelinconsistent 
        b = .!snpsumdf[!,:marker_inconsistent]
        fmiss = snpsumdf[b,col]
    else
        fmiss = snpsumdf[!,col]            
    end
    b = isnan.(fmiss)
    any(b) && deleteat!(fmiss,b)
    if isempty(fmiss)
        fig_founder = nothing
    else
        fcdf = ecdf(fmiss)        
        title = "Empirical cumulative distribution for markers in founder"    
        fig_founder = plot(x -> fcdf(x), -0.0001, 1.0001;            
            xlabel = "Fraction of missing genotypes in founders per marker",
            ylabel = "Cumulative fraction of markers",
            legend=false,            
            title,
        )
        if !isnothing(annotate_foundermiss)
            x = round(annotate_foundermiss,digits=5)
            y = round(fcdf(x),digits=5)
            scatter!(g1,[x],[y],
                series_annotations = text((x,y), :bottom)
            )
        end
    end    
    if !isnothing(mono2miss) && mono2miss 
       colls = [:offspring_miss_mono2miss]    
    else
        colls = [:offspring_miss]
    end
    # titlels = string.(colls)
    titlels =["Empirical cumulative distribution for markers in offspring"]
    gls = [begin      
        title = titlels[c]        
        col = colls[c]       
        if isdelinconsistent 
            b = .!snpsumdf[!,:marker_inconsistent]
            offmiss = snpsumdf[b,col]
        else
            offmiss = snpsumdf[!,col]            
        end        
        b = isnan.(offmiss)
        any(b) && deleteat!(offmiss,b)
        if isempty(offmiss)
            fig = nothing
        else
            offcdf = ecdf(offmiss)
            fig = plot(x -> offcdf(x), -0.0001, 1.0001;                
                xlabel = "Fraction of missing genotypes in offspring per marker",
                ylabel = "cumulative fraction of markers",
                legend=false,                
                title,
            )
            if !isnothing(annotate_offspringmiss)
                x = round(annotate_offspringmiss,digits=5)
                y = round(offcdf(x),digits=5)
                scatter!(g2,[x],[y],
                    series_annotations = text((x,y), :bottom)
                )
            end
        end
        fig
    end for c in eachindex(colls)]
    pushfirst!(gls,fig_founder)
    b = isnothing.(gls)
    any(b) && deleteat!(gls,b)
    ng = length(gls)
    plot(gls...; 
        layout =  (ng,1),   
        linewidth = 2, 
        size = (1000,800), 
        left_margin=20Plots.mm,        
        right_margin=20Plots.mm,                        
        bottom_margin=15Plots.mm,                
    )
end

function plot_mafdist(snpsumdf::AbstractDataFrame;    
    mono2miss::Union{Nothing,Bool}=true,
    isdelinconsistent::Bool = false,
    annotate_MAF::Union{Nothing,Real}=nothing)
    col = (!isnothing(mono2miss) && mono2miss)  ? :freqA2_mono2miss : :freqA2        
    if isdelinconsistent 
        b = .!snpsumdf[!,:marker_inconsistent]
        freq = snpsumdf[b,col]
    else
        freq = snpsumdf[!,col]            
    end
    freq = freq[.!isnan.(freq)]
    if isempty(freq)
        g1 = nothing
    else
        freq = [i>0.5 ? 1-i : i for i in freq]
        fcdf = ecdf(freq)        
        g1 = plot(x -> fcdf(x), -0.0001, 0.5;
            xlabel = "Minor allele frequency (MAF)",
            ylabel = "Cumulative fraction of markers",
            legend=false,
            linewidth = 2, 
            size = (1000,600),         
            left_margin=20Plots.mm,
            bottom_margin=15Plots.mm,             
            title = "Empirical cumulative distribution of MAF",
        )
        if !isnothing(annotate_MAF)
            x = round(annotate_MAF,digits=5)
            y = round(fcdf(x),digits=5)
            scatter!(g1,[x],[y],
                series_annotations = text((x,y), :bottom)
            )
        end
        g1
    end    
end
