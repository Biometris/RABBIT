
function imputeaccuracy(truegenofile::AbstractString,
    imputed_genofile::AbstractString,
    pedinfo::Union{Integer,MagicBase.JuncDist, AbstractString};
    isfounderinbred::Bool=true,
    isphysmap_true::Bool=false, 
    isphysmap_imputed::Bool=false, 
    alignfounder::Bool=true,
    formatpriority::AbstractVector=["GT","GP", "AD"],
    threshcall::Real=0.95,
    workdir::AbstractString=pwd(),
    outstem::AbstractString="outstem",
    io::Union{Nothing,IO}= nothing, 
    verbose::Bool=true)
    isdir(workdir) || @error string(workdir, " is not a directory")
    truegenofile2 = getabsfile(workdir,truegenofile)
    # ext = last(splitext(truegenofile2))
    # MagicBase.alignchromosome! will set the chromosomes of truegeno if they are missing
    truegeno = formmagicgeno(truegenofile2,pedinfo;
        isfounderinbred,formatpriority,isphysmap = isphysmap_true, recomrate = 1, 
    )    
    MagicBase.rawgenoprob!(truegeno; targets = ["founders","offspring"], baseerror = 0.001, isfounderinbred)
    MagicBase.rawgenocall!(truegeno; targets = ["founders","offspring"], callthreshold = threshcall, isfounderinbred)
    magicgeno = formmagicgeno(imputed_genofile,pedinfo; 
        isfounderinbred,formatpriority,isphysmap=isphysmap_imputed,workdir)
    MagicBase.rawgenocall!(magicgeno; targets = ["founders","offspring"], callthreshold = threshcall, isfounderinbred)
    facc = imputeaccuracy!(truegeno,magicgeno; 
        isfounderinbred, alignfounder, targetfounder=true)
    outfiles = [outstem*i for i in ["_founderacc.csv","_offspringacc.csv"]]
    open(getabsfile(workdir,outfiles[1]),"w") do outio
        descripls = ["founder, founder ID",
            "noffspring, number of offspring for the founder",
            "nmarker, number of markers",
            "miss_afterimpute, fraction of missing genotypes for the founder after imputation", 
            "ntrue, number of groud true genotypes for the founder", 
            "nonimpute, fraction of ground true genotypes that are not imputed", 
            "correctimpute, fraction of ground true genotypes that are correctly imputedd"
        ]
        for i in eachindex(descripls)
            write(outio, string("##col_",i,", ", descripls[i],"\n"))
        end        
        CSV.write(outio,facc; append=true, header=true)
    end
    printconsole(io, verbose, string("founder accuracy: ", facc))
    msg = string("save founder imputation accuracy in ", outfiles[1])
    printconsole(io, verbose, msg)    
    offacc = imputeaccuracy!(truegeno,magicgeno;isfounderinbred, targetfounder=false)
    open(joinpath(workdir,outfiles[2]), "w") do outio
        descripls = ["founder, founder ID",
            "subpop, subpopulation ID",
            "subpopsize, size of the subpopulation",
            "nmaker, number of markers", 
            "miss_afterimpute, fraction of missing genotypes for the subpopulation after imputation", 
            "ntrue, number of ground true genotypes for the subpopulation", 
            "nonimpute, fraction of ground true genotypes that are not imputed", 
            "correctimpute, fraction of ground true genotypes that are correctly imputedd"
        ]
        for i in eachindex(descripls)
            write(outio, string("##col_",i,", ", descripls[i]),"\n")
        end
        CSV.write(outio,offacc; append=true, header=true)
    end
    printconsole(io, verbose, string("offspring accuracy: ", offacc))
    msg = string("save founder imputation accuracy in ", outfiles[2])
    printconsole(io, verbose, msg)        
    facc, offacc
end

function imputeaccuracy!(truegeno::MagicGeno, magicgeno::MagicGeno;
    isfounderinbred::Bool,
    alignfounder::Bool=true,
    targetfounder::Bool=true)
    alignoffspring!(truegeno,magicgeno) 
    alignchromosome!(truegeno,magicgeno)     
    alignmarker!(truegeno,magicgeno)     
    alignfounder && alignfounder!(truegeno,magicgeno)
    alignfounder_phase!(truegeno,magicgeno)
    if targetfounder
        accuracy_mask_impute_founder!(truegeno,magicgeno; isfounderinbred)
    else
        accuracy_mask_impute_offspring!(truegeno,magicgeno)
    end
end

function accuracy_mask_impute_founder!(truegeno::MagicGeno, magicgeno::MagicGeno; isfounderinbred)
    founder_progeny = MagicBase.get_founder2offspring(magicgeno.magicped)
    ppidls = magicgeno.magicped.founderinfo[!,:individual]
    noffspringls = [length(founder_progeny[i]) for i in ppidls]
    nchr = length(magicgeno.markermap)
    nfounder =  size(magicgeno.magicped.founderinfo,1)   
    estformat = magicgeno.markermap[1][1,:founderformat]
    in(estformat,["GT_haplo","GT_phased","GT_unphased"]) || @error string("unexpected founder genoformat: ",estformat)
    res = sum([begin        
        chrtrue = truegeno.foundergeno[chr]
        chrest = magicgeno.foundergeno[chr]
        # isfounderinbred = !(typeof(chrtrue[1]) <: AbstractVector)
        # align markers
        snpidls = truegeno.markermap[chr][!,:marker]
        estsnpidls = magicgeno.markermap[chr][!,:marker]
        if snpidls  != estsnpidls
            @error "to first perform alignmarker!"            
        end
        chrres = [begin
            nsnp = size(chrest,1)
            if isfounderinbred
                nmissing = sum(chrest[:,f] .== "N")
                isnonmiss = chrtrue[:,f] .!= "N"
            else
                nmissing = sum([i == "NN" for i in join.(chrest[:,f])])
                isnonmiss = [!occursin("N",i) for i in join.(chrtrue[:,f])]
            end
            trueg = join.(chrtrue[:,f][isnonmiss])
            estg = join.(chrest[:,f][isnonmiss])
            nmask = length(trueg)
            if nmask == 0                
                nnonimpute = 0
                ncorrectimpute = 0
            else                
                if isfounderinbred
                    nnonimpute = sum(estg .== "N")
                    ncorrectimpute = sum(estg .== trueg)
                else     
                    # unphased
                    nnonimpute = sum([occursin("N",i)  for i in estg])
                    trueg2 = replace(trueg,"21"=>"12")
                    estg2 = replace(estg,"21"=>"12")
                    ncorrectimpute = sum(estg2 .== trueg2)
                
                end                
            end
            [nsnp, nmissing, nmask, nnonimpute, ncorrectimpute]
        end for f in 1:nfounder]
        chrres
    end for chr in 1:nchr])
    resdf  = DataFrame(reduce(hcat,res)', [:nmarker, :nmissing, :nmask,:nnonimpute, :ncorrectimpute])
    freq_nonimpute = resdf[!,:nnonimpute] ./ resdf[!,:nmask]
    freq_correctimpute = resdf[!,:ncorrectimpute] ./ (resdf[!,:nmask] - resdf[!,:nnonimpute])
    # nmissing = #missing genotypes in magicgeno after imputing
    # nmasked = #genotypes that are  known in truegeno but missing in magicgeno
    # accuracy_impute = among imputed masked genotypes, fraction of corectly imputed genotypes
    # freq_nonimpute = among masked genotypes, fraction of non-imputed genotypes
    DataFrame("founder"=>ppidls,
        "noffspring"=>noffspringls,
        "nmarker" =>resdf[!,:nmarker],
        "miss_afterimpute" => resdf[!,:nmissing] ./ resdf[!,:nmarker],
        "ntrue"=> resdf[!,:nmask],
        "nonimpute"=>freq_nonimpute,
        "correctimpute"=>freq_correctimpute)
end

function accuracy_mask_impute_offspring!(truegeno::MagicGeno, magicgeno::MagicGeno)    
    subpop_offs = MagicBase.get_subpop2offspring(magicgeno.magicped)
    offls  = magicgeno.magicped.offspringinfo[!,:individual]
    offdict = Dict(offls .=> 1:length(offls))
    subpopls = unique(magicgeno.magicped.offspringinfo[!,:member])
    offlsls = [[get(offdict,i,nothing) for i in subpop_offs[subpop]] for subpop in subpopls]
    nchr = length(magicgeno.markermap)
    res = sum([begin        
        # remove phasing
        chrformat = copy(magicgeno.markermap[chr][!,:offspringformat])
        # println("chr=",chr, ",est offformat=",unique(chrformat))
        if in("GT_phased",chrformat)            
            chrformat = [i == "GT_phased" ? "GT_unphased" : i for i in chrformat]
            chrest= replace([join(i) for i=magicgeno.offspringgeno[chr]],"21"=>"12")
        else
            chrest = replace(magicgeno.offspringgeno[chr],"21"=>"12")        
        end
        chrformat = copy(truegeno.markermap[chr][!,:offspringformat])
        # println("chr=",chr, ",true offformat=",unique(chrformat))
        chrformat = [i == "GT_phased" ? "GT_unphased" : i for i in chrformat]
        chrtrue= replace([join(i) for i=truegeno.offspringgeno[chr]],"21"=>"12")        
        # align markers
        snpidls = truegeno.markermap[chr][!,:marker]
        estsnpidls = magicgeno.markermap[chr][!,:marker]
        if snpidls  != estsnpidls
            @error "to first perform alignmarker!"            
        end  
        # println("true gset=",unique(chrtrue),",est gset=",unique(chrest))      
        [begin
            nsnp = size(chrest,1)
            missingcodes = ["NN","1N","2N","N1","N2"]
            nmissing = sum([in(i, missingcodes) for i in chrest[:,offls]])
            ismiss = truegeno_ismiss(chrtrue[:,offls],chrformat)
            isnonmiss = .!ismiss
            trueg = chrtrue[:,offls][isnonmiss]
            estg = chrest[:,offls][isnonmiss]
            nmask = length(trueg)
            nnonimpute = isempty(estg) ? 0 : sum([in(i, missingcodes) for i in estg])
            ncorrectimpute =  isempty(estg) ? 0 : sum(estg .== trueg)
            [nsnp, nmissing, nmask, nnonimpute, ncorrectimpute]
        end for offls in offlsls]
    end for chr in 1:nchr])
    resdf  = DataFrame(reduce(hcat,res)', [:nmarker, :nmissing, :nmask,:nnonimpute, :ncorrectimpute])
    freq_nonimpute = resdf[!,:nnonimpute] ./ resdf[!,:nmask]
    freq_correctimpute = resdf[!,:ncorrectimpute] ./ (resdf[!,:nmask] - resdf[!,:nnonimpute])
    subpopsizels = length.(offlsls)
    # nmissing = #missing genotypes in magicgeno after imputing
    # nmasked = #genotypes that are  known in truegeno but missing in magicgeno
    # accuracy_impute = among imputed masked genotypes, fraction of corectly imputed genotypes
    # freq_nonimpute = among masked genotypes, fraction of non-imputed genotypes
    DataFrame("subpop"=>subpopls,
        "subpopsize"=>subpopsizels,
        "nmarker" =>resdf[!,:nmarker],
        "miss_afterimpute" => resdf[!,:nmissing] ./ (resdf[!,:nmarker] .* subpopsizels),
        "ntrue"=> resdf[!,:nmask],
        "nonimpute"=>freq_nonimpute,
        "correctimpute"=>freq_correctimpute)
end
