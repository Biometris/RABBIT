
function magicaccuracy(truefile::AbstractString,estfile::AbstractString,
    pedinfo::Union{Integer,AbstractString};  
    isfounderinbred::Bool=true,  
    commentstring::AbstractString="##",
    workdir::AbstractString=pwd(),
    isphysmap_true::Bool=false, 
    alignfounder::Bool=true,
    chrsubset::Union{Nothing,AbstractVector}=nothing,
    minprob::Real=0)
    truegeno= formmagicgeno(truefile,pedinfo; isfounderinbred, 
        isphysmap = isphysmap_true, recomrate = 1.0, 
        formatpriority=["GT"], workdir
    )
    estfile2=getabsfile(workdir,estfile)
    ext = last(splitext(estfile2))
    myopen = ext == ".gz" ? GZip.open : open
    line1 = myopen(estfile2,"r") do io
       readline(io)
    end
    if occursin("RABBIT,designinfo",line1)
        # TODO: robust test if ancestry file
       magicest = readmagicancestry(estfile2)
       offkind = only(unique(vcat([unique(i[!,:offspringformat]) for i=truegeno.markermap]...)))
       if offkind != "discretefgl"
           msg = "genotypes in truefile are not discretefgl"
           error(msg)
       end
    else
        magicest = formmagicgeno(estfile2,pedinfo; isfounderinbred, formatpriority=["GT"], commentstring)
        offkind = unique(vcat([unique(i[!,:offspringformat]) for i=truegeno.markermap]...))
        if !issubset(offkind, ["GT_phased","GT_unphased"])
            msg = string("offspring formats: ", offkind, " are in GT format")
            error(msg)
        end
    end
    magicaccuracy!(truegeno, magicest; isfounderinbred,alignfounder,chrsubset,minprob)
end

function magicaccuracy!(truegeno::MagicGeno,magicest::Union{MagicGeno,MagicAncestry};
    isfounderinbred::Bool=true,  
    alignfounder::Bool=true,
    chrsubset::Union{Nothing,AbstractVector}=nothing,
    minprob::Real=0)        
    alignoffspring!(truegeno,magicest) # revise offspring of truegeno to be consistent with magicest, aumming offspring in magicest must be in trugeno
    alignchromosome!(truegeno,magicest)     
    nchr = length(truegeno.markermap)
    chrsubset = isnothing(chrsubset) ? (1:nchr) : intersect(chrsubset,1:nchr)
    truesnps = reduce(vcat, [i[!,:marker] for i in truegeno.markermap[chrsubset]])
    estsnps =  reduce(vcat,[i[!,:marker] for i in magicest.markermap[chrsubset]])
    # nsnptrue: #markers in truegeno
    # nungroup: #markers in truegeno but not in magicest
    # nunknown: #markers in magicest but not in truegeno
    # nmisgroup: #markers in different chromosomes between magicest and truegeno
    nsnptrue = length(truesnps)
    nungroup = length(setdiff(truesnps, estsnps))
    unknown = setdiff(estsnps, truesnps)
    nunknown = length(unknown)
    nmisgroup = sum([length(setdiff(magicest.markermap[chr][!,:marker], 
        truegeno.markermap[chr][!,:marker],unknown)) for chr in chrsubset])
    snpsum = (nsnptrue=nsnptrue, nungroup=nungroup,nunknown=nunknown,nmisgroup=nmisgroup)
    # align markers 
    alignmarker!(truegeno,magicest) # afterward, same sets of markers
    # align founders
    if isa(magicest, MagicGeno)
        MagicBase.rawgenoprob!(magicest; targets = ["founders"], seqerror=0.001, isfounderinbred)
        MagicBase.rawgenocall!(magicest; targets = ["founders"], callthreshold = 0.95, isfounderinbred)
    end
    alignfounder && alignfounder!(truegeno,magicest)
    alignfounder_phase!(truegeno,magicest)    
    parentacc = calparentacc(truegeno, magicest; chrsubset)
    if typeof(magicest) <: MagicGeno
        offacc = caloffspringacc(truegeno, magicest; chrsubset)
    else
        offacc = caloffspringacc(truegeno, magicest; chrsubset,minprob)
    end    
    Dict("marker"=>snpsum, "founder"=>parentacc, "offspring"=>offacc)
end



function caloffspringacc(truegeno::MagicGeno,magicancestry::MagicAncestry;
    chrsubset::AbstractVector,minprob::Real=0)    
    isviterbi = !isnothing(magicancestry.viterbipath)    
    if isviterbi
        caloffspringacc_viterbi(truegeno, magicancestry; chrsubset)
    else
        caloffspringacc_fw(truegeno, magicancestry; chrsubset,minprob)
    end    
end

function caloffspringacc_viterbi(truegeno::MagicGeno,magicancestry::MagicAncestry;
    chrsubset::Union{Nothing,AbstractVector}=nothing)
    isviterbi = !isnothing(magicancestry.viterbipath)    
    isviterbi || error("viterbipath is nothing")
    is_path_diplo,path_fgls = MagicBase.get_path_fgls(magicancestry)
    isnothing(chrsubset)
    if is_path_diplo
        # ignore phase information
        nfgl = length(magicancestry.statespace["haplotype"])
        stateindex2 = MagicBase.prior_diploindex(nfgl)
        stateindex2 == path_fgls || @error "inconsistent viterbi states"
        viterbipath = magicancestry.viterbipath 
    else            
        stateindex = path_fgls
        stateindex2 = [[i,i] for i in stateindex]
        viterbipath = magicancestry.viterbipath 
    end
    if isnothing(chrsubset)
        chrsubset = 1:length(magicancestry.markermap)
    end
    res = [begin 
        # check markers
        snpidls = truegeno.markermap[chr][!,:marker]
        estsnpidls = magicancestry.markermap[chr][!,:marker]
        if snpidls  != estsnpidls
            @error "to first perform aligntruegeno!"
        end
        chrmap = magicancestry.markermap[chr]
        chrlen = chrmap[end,:poscm]-chrmap[1,:poscm]
        chrid = chrmap[1,:linkagegroup]        
        ischrx = chrid in ["chrx"]
        ischrx && @error string("TODO: chrx")
        # if ischrx
        #     chrtrue = [i[2] == "Y" ? [i[1],i[1]] : i for i=chrtrue]  # TODO for chrest        
        # end
        chrtrue = truegeno.offspringgeno[chr]
        chrest = stateindex2[viterbipath[chr] ]        
        # calculate phasing error
        nswitcherr,nheteroerr = sum(MagicBase.calphaseacc(chrtrue[:,off],chrest[:,off]) for off in 1:size(chrtrue,2))
        # calculate calling error
        nmiss = 0
        bnotequal = map((x,y) -> sort(x) != sort(y), chrtrue,chrest)
        ncallerr = isempty(bnotequal) ? 0 : (sum(bnotequal) - nmiss)        
        ncall = length(chrtrue) -  nmiss
        avgnrecom = sum(MagicBase.calnumrecom([chrest]))/size(chrest,2)        
        assignerr = NaN
        [ncallerr, nmiss, ncall, assignerr, nswitcherr,nheteroerr,avgnrecom,chrlen]        
    end for chr in chrsubset]
    ressum = sum(res)    
    callerr = round(ressum[1]/ressum[3],digits=6)
    callmiss=round(ressum[2]/(ressum[2]+ressum[3]),digits=6)    
    assignerr = NaN
    switcherr = round(/(ressum[5]...),digits=6)
    heteroerr = round(/(ressum[6]...),digits=6)
    avgnrecom = round(ressum[7],digits=2) # in cM
    genlen =round(ressum[8],digits=2) # average #recombation per offspring
    (ncallerr=ressum[1], ncalledorgin=ressum[3], callerr=callerr,callmiss=callmiss,assignerr=assignerr,
        phase_nswitcherr=ressum[5][1], phase_switcherr=switcherr, 
        phase_nheteroerr=ressum[6][1],phase_heteroerr=heteroerr,
        avgnrecom=avgnrecom,genomelen=genlen)
end

function caloffspringacc_fw(truegeno::MagicGeno,magicancestry::MagicAncestry;
    chrsubset::Union{Nothing,AbstractVector}=nothing,
    minprob = 0.7)    
    isdiploprob = !isnothing(magicancestry.diploprob)
    if isdiploprob    
        nfgl = length(magicancestry.statespace["haploindex"])
        haploindex = magicancestry.statespace["haploindex"]    
        diploindex = MagicBase.prior_diploindex(nfgl)            
        genoindex = MagicBase.prior_genoindex(nfgl)            
        genoindexdict = Dict(genoindex .=> 1:length(genoindex))
    else
        isnothing(magicancestry.haploprob) &&  @error string("condprob is nothing in magicancestry")    
        haploindex = magicancestry.statespace["haploindex"]    
    end
    if isnothing(chrsubset)
        chrsubset = 1:length(magicancestry.markermap)
    end
    res = [begin 
        # check markers
        snpidls = truegeno.markermap[chr][!,:marker]
        estsnpidls = magicancestry.markermap[chr][!,:marker]
        if snpidls  != estsnpidls
            @error "to first perform aligntruegeno!"
        end
        chrmap = magicancestry.markermap[chr]
        chrlen = chrmap[end,:poscm]-chrmap[1,:poscm]
        chrid = chrmap[1,:linkagegroup]        
        ischrx = chrid in ["chrx"]
        ischrx && @error string("TODO: chrx")
        # if ischrx
        #     chrtrue = [i[2] == "Y" ? [i[1],i[1]] : i for i=chrtrue]  # TODO for chrest        
        # end
        chrtrue = truegeno.offspringgeno[chr]        
        if isdiploprob
            # calculate assignerr
            nsnp,noff = size(chrtrue)
            chrgenoprob = magicancestry.genoprob[chr]
            chrtrue2 = [genoindexdict[sort(i,rev=true)] for i in chrtrue]
            assignp = [chrgenoprob[off][m,chrtrue2[m,off]] for m=1:nsnp, off=1:noff]
            assignerr = 1-sum(assignp)/length(assignp)
            # calcualte phasing error
            chrcall = first(MagicBase.ancestrycall([magicancestry.diploprob[chr]],minprob=minprob))
            chrcall2 = [ismissing(i) ? [-1,-1] : diploindex[i] for i in chrcall]
            nswitcherr,nheteroerr = sum(MagicBase.calphaseacc(chrtrue[:,off],chrcall2[:,off]) for off in 1:size(chrtrue,2))
            # calculate callerr
            chrcall = first(MagicBase.ancestrycall([magicancestry.genoprob[chr]],minprob=minprob))            
            avgnrecom = sum(MagicBase.calnumrecom([chrcall]))/noff
            b = ismissing.(chrcall)
            nmiss = sum(b)
            chrcall2 = [ismissing(i) ? [-1,-1] : genoindex[i] for i in chrcall]
            bnotequal = map((x,y) -> sort(x) != sort(y), view(chrtrue,.!b),view(chrcall2,.!b))
            ncallerr = isempty(bnotequal) ? 0 : sum(bnotequal)
            ncall = length(bnotequal)
        else
            # calculate assignerr
            nsnp,noff = size(chrtrue)
            chrtrue2 = [allequal(i) ? first(i) : missing for i in chrtrue]
            chrhaploprob = magicancestry.haploprob[chr]                                    
            assignp = [ismissing(chrtrue2[m,off]) ? missing : chrhaploprob[off][m,chrtrue2[m,off]] for m=1:nsnp, off=1:noff]
            isgood = .![ismissing(i) || isnan(i) for i in assignp]
            assignerr = 1-sum(assignp[isgood])/sum(isgood)            
            # calcualte phasing error
            nswitcherr,nheteroerr = [0,0], [0,0]
            # calculate callerr
            chrcall = first(MagicBase.ancestrycall([chrhaploprob],minprob=minprob))            
            avgnrecom = sum(MagicBase.calnumrecom([chrcall]))/noff
            b = ismissing.(chrcall) .|| ismissing.(chrtrue2)
            nmiss = sum(b)            
            bnotequal = map(!=, view(chrtrue2,.!b),view(chrcall,.!b))
            ncallerr = isempty(bnotequal) ? 0 : sum(bnotequal)
            ncall = length(bnotequal)
        end
        [ncallerr, nmiss, ncall, assignerr, nswitcherr,nheteroerr,avgnrecom,chrlen]    
    end for chr in chrsubset]    
    ressum = sum(res)    
    callerr = round(ressum[1]/ressum[3],digits=6)
    callmiss=round(ressum[2]/(ressum[2]+ressum[3]),digits=6)    
    assignerr = round(ressum[4]/length(res), digits=6)
    switcherr = round(/(ressum[5]...),digits=6)
    heteroerr = round(/(ressum[6]...),digits=6)
    avgnrecom = round(ressum[7],digits=2) # in cM
    genlen =round(ressum[8],digits=2) # average #recombation per offspring
    (ncallerr=ressum[1], ncalledorigin=ressum[3], callerr=callerr,callmiss=callmiss,assignerr=assignerr,
        phase_nswitcherr=ressum[5][1], phase_switcherr=switcherr, 
        phase_nheteroerr=ressum[6][1],phase_heteroerr=heteroerr,
        avgnrecom=avgnrecom,genomelen=genlen)
end

function caloffspringacc(truegeno::MagicGeno,magicest::MagicGeno;
    chrsubset::AbstractVector)
    check_trueformat(truegeno; isfouderformat=false)    
    # check only GT format
    res = [begin
        # check markers
        snpidls = truegeno.markermap[chr][!,:marker]
        estsnpidls = magicest.markermap[chr][!,:marker]
        if snpidls  != estsnpidls
            @error "to first perform aligntruegeno!"
        end        
        chrest = magicest.offspringgeno[chr]
        chrtrue = truegeno.offspringgeno[chr]
        true_chrformat = truegeno.markermap[chr][!,:offspringformat]
        est_chrformat = magicest.markermap[chr][!,:offspringformat]
        phasedsnps = findall(@. true_chrformat == "GT_phased" && est_chrformat == "GT_phased")
        if isempty(phasedsnps)
            nswitcherr,nheteroerr = [0, 0], [0,0]
        else
            nswitcherr,nheteroerr = sum(MagicBase.calphaseacc(chrtrue[phasedsnps,off],chrest[phasedsnps,off]) for off in 1:size(chrtrue,2))
        end
        if in("GT_phased",est_chrformat)                        
            est_chrformat = [i == "GT_phased" ? "GT_unphased" : i for i in est_chrformat]
            chrest = replace(join.(chrest),"21"=>"12")            
        else
            chrest = replace(chrest,"21"=>"12")
        end
        true_chrformat = [i == "GT_phased" ? "GT_unphased" : i for i in true_chrformat]
        chrtrue = replace(join.(chrtrue),"21"=>"12")
        # calculate error rate
        ismiss = truegeno_ismiss(chrtrue,true_chrformat)
        isnonmiss = .!ismiss # nonmiss bool for truegeno
        bnotequal = map(!=, chrtrue[isnonmiss],chrest[isnonmiss])
        nmiss = sum([in(i, ["NN","1N","2N","N1","N2"]) for i=chrest[isnonmiss][bnotequal]])
        ngenoerr = isempty(bnotequal) ? 0 : (sum(bnotequal) - nmiss)        
        noffspringgeno = sum(isnonmiss) -  nmiss
        # println("[chr, ngenoerr, nmiss, noffspringgeno, nswitcherr,nheteroerr]=",[chr,ngenoerr, nmiss, noffspringgeno, nswitcherr,nheteroerr])
        [ngenoerr, nmiss, noffspringgeno, nswitcherr,nheteroerr]        
    end for chr=chrsubset]    
    ressum = sum(res)    
    offspringerr = round(ressum[1]/ressum[3],digits=6)
    offspringmiss=round(ressum[2]/(ressum[2]+ressum[3]),digits=6)    
    switcherr = round(/(ressum[4]...),digits=6)
    heteroerr = round(/(ressum[5]...),digits=6)
    (ngenoerr=ressum[1],ncalledgeno=ressum[3], 
        genoerr=offspringerr,genomiss=offspringmiss,
        phase_nswitcherr=ressum[4][1], phase_switcherr=switcherr, 
        phase_nheteroerr=ressum[5][1],phase_heteroerr=heteroerr)
end

function calparentacc(truegeno::MagicGeno,magicest::Union{MagicGeno,MagicAncestry};
    chrsubset::AbstractVector)
    check_trueformat(truegeno; isfouderformat=true)
    # TODO:  phased parental genotypes, nphaseerr    
    res = [begin
        # check markers
        snpidls = truegeno.markermap[chr][!,:marker]
        estsnpidls = magicest.markermap[chr][!,:marker]
        if snpidls  != estsnpidls
            @error "to first perform aligntruegeno!"
        end 
        chrest = magicest.foundergeno[chr]
        chrtrue = truegeno.foundergeno[chr]
        true_chrformat = truegeno.markermap[chr][!,:founderformat]
        est_chrformat = magicest.markermap[chr][!,:founderformat]
        phasedsnps = findall(@. true_chrformat == "GT_phased" && est_chrformat == "GT_phased")
        if isempty(phasedsnps)
            nswitcherr,nheteroerr = [0, 0], [0,0]
        else
            nswitcherr,nheteroerr = sum(MagicBase.calphaseacc(chrtrue[phasedsnps,off],chrest[phasedsnps,off]) for off in 1:size(chrtrue,2))
        end
        if in("GT_phased",est_chrformat)                        
            est_chrformat = [i == "GT_phased" ? "GT_unphased" : i for i in est_chrformat]
            chrest = replace(join.(chrest),"21"=>"12")
        else
            chrest = replace(chrest,"21"=>"12")
        end
        true_chrformat = [i == "GT_phased" ? "GT_unphased" : i for i in true_chrformat]
        chrtrue = replace(join.(chrtrue),"21"=>"12")
        # calculate error rate
        ismiss = truegeno_ismiss(chrtrue,true_chrformat)

        isnonmiss = .!ismiss # nonmiss bool for truegeno
        # bnotequal = map(!=, chrtrue[isnonmiss],chrest[isnonmiss])
        estformat = est_chrformat[1]
        if estformat == "GT_haplo"
            bnotequal = map(!=, view(chrtrue,isnonmiss),view(chrest,isnonmiss))                
        elseif estformat == "GT_phased"
            chrtrue_nonmiss = view(chrtrue,isnonmiss)
            chrest_nonmiss = view(chrest,isnonmiss)
            bnotequal = map(!=, chrtrue_nonmiss,chrest_nonmiss)            
        elseif estformat == "GT_unphased"
            chrtrue_nonmiss2 = replace(join.(view(chrtrue,isnonmiss)),"21"=>"12")
            chrest_nonmiss2 = replace(view(chrest,isnonmiss),"21"=>"12")
            bnotequal = map(!=, chrtrue_nonmiss2,chrest_nonmiss2)            
        else
            @error string("unexpected founder genoformat: ",estformat)
        end
        nmiss = sum([in(i, ["N", "NN","1N","2N","N1","N2"]) for i=chrest[isnonmiss][bnotequal]])
        ngenoerr = isempty(bnotequal) ? 0 : (sum(bnotequal) - nmiss)
        nparentgeno = sum(isnonmiss) -  nmiss    
        # println("[chr, ngenoerr, nmiss, noffspringgeno, nswitcherr,nheteroerr]=",[chr,ngenoerr, nmiss, noffspringgeno, nswitcherr,nheteroerr])
        [ngenoerr, nmiss, nparentgeno, nswitcherr,nheteroerr]      
    end for chr=chrsubset]
    ressum = sum(res)    
    foundererr = round(ressum[1]/ressum[3], digits=6)
    foundermiss=round(ressum[2]/(ressum[3]+ressum[2]),digits=6)       
    switcherr = round(/(ressum[4]...),digits=6)
    heteroerr = round(/(ressum[5]...),digits=6) 
    (ngenoerr=ressum[1],ncalledgeno=ressum[3], 
        genoerr=foundererr, genomiss=foundermiss,
        phase_nswitcherr=ressum[4][1], phase_switcherr=switcherr, 
        phase_nheteroerr=ressum[5][1], phase_heteroerr=heteroerr)
end


function calphaseacc(true_phasedgeno::AbstractVector,est_phasedgeno::AbstractVector)
    nsnp = length(true_phasedgeno)
    nsnp == length(est_phasedgeno) || @error string("inconsistent number of markers")
    alleletype = eltype(first(true_phasedgeno))
    if alleletype <: Integer
        ishetero = @. allunique(true_phasedgeno) && allunique(est_phasedgeno)
    elseif alleletype <: AbstractString
        ishetero = [issetequal(["1","2"],true_phasedgeno[i]) && issetequal(["1","2"],est_phasedgeno[i]) for i in 1:nsnp]    
    else
        @error string("unexpected allele data type = ",alleletype, " for truegeno=",true_phasedgeno) maxlog=20
    end
    nhetero = sum(ishetero)
    nhetero == 0 && return [[0,0],[0,0]]
    truehaplo = true_phasedgeno[ishetero]
    esthaplo = est_phasedgeno[ishetero]
    nheteroerr = sum(esthaplo .!= truehaplo)
    esthaplo2 = reverse.(esthaplo)
    nheteroerr2 = sum(esthaplo2 .!= truehaplo)
    if nheteroerr > nheteroerr2
        nheteroerr = nheteroerr2
        esthaplo = esthaplo2
    end    
    segls = MagicBase.splitindex(tuple.(truehaplo,esthaplo))
    posls = first.(segls)
    trueh = truehaplo[posls]
    esth = deepcopy(esthaplo[posls])
    nswitcherr = 0
    nseg = length(posls)
    for i in 1:nseg
        if trueh[i] != esth[i]
            # esth[i:nseg] .= reverse.(esth[i:nseg])
            reverse!.(view(esth,i:nseg))  # reverse! will modify esthaplo (and est_phasedgeno and magicest) if deepcopy  is not used            
            nswitcherr += 1           
        end
    end
    if alleletype <: AbstractString
        trueh == esth || @error string("unexpected mismatch after switching")    
    end
    [[nswitcherr,nhetero-1],[nheteroerr,nhetero]]
end


function check_trueformat(truegeno::MagicGeno;
    isfouderformat::Bool=true)
    formatset = ["GT_haplo","GT_unphased","GT_phased"]
    nchr = length(truegeno.markermap)
    formatcol = isfouderformat ? :founderformat : :offspringformat
    res = trues(nchr)
    for chr in 1:nchr
        f = unique(truegeno.markermap[chr][!,formatcol])
        d = setdiff(f, formatset)
        if !isempty(d)
            @error string("non ", " GT formats = ",
                d, " in", isfouderformat ? "founders" : "offspring", " in chr=",chr)
            res[chr] = false
        end
    end
    all(res)
end




function truegeno_ismiss(chrgeno::AbstractMatrix, chrformat::AbstractVector)
    ismissls = falses(size(chrgeno)...)
    formatset = unique(chrformat)
    for format in formatset
        ii = chrformat .== format
        subgeno = view(chrgeno, ii,:)
        if format == "GT_haplo"
            g = setdiff(unique(subgeno),["1","2"])
            ismissN = occursin.("N",g)
        elseif format == "GT_phased"            
            g = setdiff(unique(subgeno),[["1","1"],["1","2"],["2","1"],["2","2"]])
            ismissN = [in("N",i) for i in g]
        else
            format == "GT_unphased" || @error string("unexpected format: ",format)
            g = setdiff(unique(subgeno),["11","12","21","22"])
            ismissN = occursin.("N",g)
        end        
        all(ismissN) || @error string("unknow genotypes: ", g[.!ismissN])
        if !isempty(g)
            b = [in(i,g) for i in subgeno]
            ismissls[ii,:] .= b
        end
    end
    ismissls
end

function alignoffspring!(truegeno::MagicGeno, magicest::Union{MagicGeno,MagicAncestry})
    # order individuals of truegeno, delete individauls that not in magicest
    n=size(truegeno.magicped.founderinfo,1)
    dict = Dict(truegeno.magicped.founderinfo[!,:individual] .=> 1:n)
    indls = [get(dict, i, nothing) for i=magicest.magicped.founderinfo[!,:individual]]
    if nothing in indls
        b = isnothing.(indls)
        error(string("founders in truegeno but not in magicest: ", indls[b]))
    else
        if indls  != 1:n
            truegeno.magicped.founderinfo = truegeno.magicped.founderinfo[indls,:]
            truegeno.foundergeno = [i[:,indls] for i = truegeno.foundergeno]
        end
    end
    magicest.magicped.founderinfo[!,:individual] == truegeno.magicped.founderinfo[!,:individual] || @error "inconsisent founders"
    n=size(truegeno.magicped.offspringinfo,1)
    dict = Dict(truegeno.magicped.offspringinfo[!,:individual] .=> 1:n)
    indls = [get(dict, i, nothing) for i=magicest.magicped.offspringinfo[!,:individual]]
    if nothing in indls
        b = isnothing.(indls)
        error(string("offspring in truegeno but not in magicest: ", indls[b]))
    else
        if indls  != 1:n
            truegeno.magicped.offspringinfo = truegeno.magicped.offspringinfo[indls,:]
            truegeno.offspringgeno = [i[:,indls] for i = truegeno.offspringgeno]
        end
    end
    magicest.magicped.offspringinfo[!,:individual] == truegeno.magicped.offspringinfo[!,:individual] || @error "inconsisent offspring"
    truegeno
end

function alignchromosome!(truegeno::MagicGeno, magicest::Union{MagicGeno,MagicAncestry})
    # align chromosome groups
    truelg = findtruelg(truegeno.markermap,magicest.markermap)
    expectlg = [[i] for i in 1:length(truegeno.markermap)]
    if truelg != expectlg
        truegeno.markermap = [reduce(vcat,truegeno.markermap[i]) for i in truelg]
        truegeno.foundergeno = [reduce(vcat,truegeno.foundergeno[i]) for i in truelg]
        truegeno.offspringgeno = [reduce(vcat,truegeno.offspringgeno[i]) for i in truelg]
    end    
    truegeno
end


function findtruelg(truemarkermap, estmarkermap; minfreq=0.2, isphysmap = [false, false], verbose=true)
    # true chromosome IDs for each stimated linakge group
    snpinter = [length(intersect(i[!,:marker],j[!,:marker]))
        for i in truemarkermap, j in estmarkermap]                
    # @info "" snpinter
    truelg = [begin
        pos = findall(snpinter[:,chr] .> 1)
        ninterls = snpinter[pos,chr]
        nsum = sum(ninterls)
        pos2 = pos[[n/nsum >= minfreq && n > 5 for n in ninterls]]                
        pos2 = union(pos2,argmax(snpinter[:,chr]))        
        pos2
    end for chr in 1:size(snpinter,2)]        
    chrcol = [i ? "physchrom" : "linkagegroup" for i in isphysmap]
    truechridls = convert.(String,[string(i[1,chrcol[1]]) for i in truemarkermap])
    estchridls = convert.(String,[string(i[1,chrcol[2]]) for i in estmarkermap])
    pos = findall(length.(truelg) .> 1)
    if !isempty(pos)        
        truels = [truechridls[i] for i in truelg[pos]]
        msgls = map((x,y) -> string(x, " => ", join(y,",")), estchridls[pos], truels)
        verbose && @warn join(msgls,"; ")
    end
    ls = vcat(truelg...)
    msgls = []
    for chr in 1:length(truemarkermap)
        b = ls .== chr
        if sum(b)>1
            estindices = findall([in(chr,i) for i in truelg])            
            msg = string(truechridls[chr] , " => ", join(estchridls[estindices],","))            
            push!(msgls,msg)
        end
    end    
    isempty(msgls) || (verbose && @warn join(msgls,"; "))
    pos = findall(isempty.(truelg))
    if !isempty(pos)        
        @error string("chromosomes not found: ", estchridls[pos])
    end
    truelg
end



function alignmarker!(truegeno::MagicGeno, magicest::Union{MagicGeno,MagicAncestry})
    for chr in eachindex(truegeno.markermap)                
        snpidls = truegeno.markermap[chr][!,:marker]
        estsnpidls = magicest.markermap[chr][!,:marker]
        if snpidls  != estsnpidls
            rule = Dict(snpidls .=> 1:length(snpidls))
            ii = [get(rule, i, missing) for i in estsnpidls]
            ii2 = collect(skipmissing(ii))                
            truegeno.markermap[chr] = truegeno.markermap[chr][ii2,:]
            truegeno.foundergeno[chr] = truegeno.foundergeno[chr][ii2,:]
            truegeno.offspringgeno[chr] = truegeno.offspringgeno[chr][ii2,:]
            nonmissls = .!ismissing.(ii)
            if !all(nonmissls)                    
                magicest.markermap[chr] = magicest.markermap[chr][nonmissls, :]
                magicest.foundergeno[chr] = magicest.foundergeno[chr][nonmissls, :]
                if isa(magicest,MagicGeno)
                    magicest.offspringgeno[chr] = magicest.offspringgeno[chr][nonmissls, :]
                elseif isa(magicest,MagicAncestry)                    
                    isviterbi = !isnothing(magicest.viterbipath)    
                    if isviterbi                        
                        # viterbipath[chr][:,offspring]: most probably path for offspring at chr.
                        magicest.viterbipath[chr] = magicest.viterbipath[chr][nonmissls, :]
                    else
                        for condprob in [magicest.diploprob,magicest.genoprob,magicest.haploprob]
                            if !isnothing(condprob)
                                chrprob = condprob[chr]
                                for o in eachindex(chrprob)
                                    chrprob[o] = chrprob[o][nonmissls,:]
                                end   
                            end
                        end     
                    end           
                else
                    @error "wrong type of magicest"
                end
            end            
        end       
    end 
    truegeno, magicest
end

function alignfounder!(truegeno::MagicGeno, magicest::Union{MagicGeno,MagicAncestry})    
    for chr in eachindex(truegeno.markermap)                
        perm = permfounder(truegeno.foundergeno[chr],magicest.foundergeno[chr])                        
        if perm != 1:length(perm)
            truegeno.foundergeno[chr] = truegeno.foundergeno[chr][:,perm]                                    
            if unique(truegeno.markermap[chr][!,:offspringformat]) == ["discretefgl"]                
                isfounderinbred = !(typeof(first(truegeno.foundergeno[chr])) <: AbstractVector)
                if isfounderinbred
                    dict= Dict(perm .=> 1:length(perm))
                    nfgl = length(perm)                    
                else
                    rules = perm .=> 1:length(perm)
                    dict = Dict(reduce(vcat, [[2*i-1,2*i] .=> [2*j-1,2*j] for (i,j) in rules]))
                    nfgl = 2*length(perm)
                end
                diploindex = MagicBase.prior_diploindex(nfgl)
                dict2 = Dict([j => [get(dict, i, nothing) for i in j] for j in diploindex])
                chrtrue = truegeno.offspringgeno[chr]
                for i in eachindex(chrtrue)
                    chrtrue[i] = get(dict2,chrtrue[i],nothing)
                end
                (nothing in chrtrue) && @error string("wrong chrtrue in alignfounder")
            end
        end
        # println("ch=",chr, ", perm=",perm, ", d=", sum(truegeno.foundergeno[chr] .!== magicest.foundergeno[chr]))
    end
    truegeno
end


function permfounder(chrftrue::AbstractMatrix,chrfest::AbstractMatrix)
    ncol = size(chrfest,2)
    perm = Vector{Int}()
    isfounderinbred = !(typeof(first(chrftrue)) <: AbstractVector)
    for i in eachcol(chrfest)
        cols = setdiff(1:ncol,perm)
        if isfounderinbred
            dls = [sum(i .!= chrftrue[:,j])  for j in cols]
        else
            dls = [begin 
                d = sum(i .!= chrftrue[:,j])  
                d2 = sum(i .!= reverse.(chrftrue[:,j]))  
                d = min(d,d2)
            end for j in cols]
        end
        push!(perm,cols[argmin(dls)])
    end
    perm
end

function alignfounder_phase!(chrftrue::AbstractMatrix,chrfest::AbstractMatrix)
    rules = []
    for f in 1:size(chrftrue,2)
        col_est = view(chrfest,:,f)
        col_true = view(chrftrue,:, f)
        col_rev_true = reverse.(col_true)
        d1 = sum(col_true .!= col_est)
        d2 = sum(col_rev_true .!= col_est)
        if d2 < d1
            chrftrue[:,f] .= col_rev_true
            append!(rules,[2*f-1,2*f] .=> [2*f, 2*f-1])
        end    
    end
    rules
end

function alignfounder_phase!(truegeno::MagicGeno, magicest::Union{MagicGeno,MagicAncestry})    
    isfounderinbred = !(typeof(truegeno.foundergeno[1][1]) <: AbstractVector)
    isfounderinbred && return truegeno
    for chr in eachindex(truegeno.markermap)                            
        chrftrue = truegeno.foundergeno[chr]
        chrfest = magicest.foundergeno[chr]
        rules = alignfounder_phase!(chrftrue, chrfest)
        if unique(truegeno.markermap[chr][!,:offspringformat]) == ["discretefgl"]          
            replace!.(truegeno.offspringgeno[chr],rules...)
        end
    end
    truegeno
end