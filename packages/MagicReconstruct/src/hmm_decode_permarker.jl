
function callogbackward_permarker!(decodetempfile, chrfhaplo::AbstractMatrix,chroffgeno::AbstractMatrix,
    popmakeup::AbstractDict,priorprocess::AbstractDict;    
    epsf::Union{Real,AbstractVector},
    epso::Union{Real,AbstractVector},    
    epso_perind::Union{Nothing,AbstractVector},     
    seqerror::Union{Real,AbstractVector},
    allelebalancemean::Union{Real,AbstractVector},
    allelebalancedisperse::Union{Real,AbstractVector},
    alleledropout::Union{Real,AbstractVector},    
    israndallele::Bool, 
    issnpGT::AbstractVector,    
    snporder::Union{Nothing,AbstractVector}=nothing)    
    nsnp,noff = size(chroffgeno)    
    isnothing(snporder) && (snporder = 1:nsnp)    
    epsfls = typeof(epsf) <: Real ? epsf*ones(nsnp) : epsf
    epsols = typeof(epso) <: Real ? epso*ones(nsnp) : epso
    seqerrorls = typeof(seqerror) <: Real ? seqerror*ones(nsnp) : seqerror
    allelebalancemeanls = typeof(allelebalancemean) <: Real ? allelebalancemean*ones(nsnp) : allelebalancemean
    allelebalancedispersels = typeof(allelebalancedisperse) <: Real ? allelebalancedisperse*ones(nsnp) : allelebalancedisperse
    alleledropoutls = typeof(alleledropout) <: Real ? alleledropout*ones(nsnp) : alleledropout
    markerincl=first(values(priorprocess)).markerincl
    tseq = findall(markerincl)    
    hmmalg = "logbackward"
    jldopen(decodetempfile,"w") do file 
        write(file, "hmmalg", hmmalg)
        write(file, "snporder", snporder)        
        write(file, "noffspring", noff)    
        write(file, "tseq", tseq)          
        dataprobls = MagicReconstruct.init_dataprobls_singlephase(popmakeup)    
        logbwprob = [zeros(length(i)) for i in dataprobls]          
        write(file, string("t",tseq[end]),logbwprob)        
        for kk in length(tseq)-1:-1:1    
            tnow = tseq[kk]        
            tback = tseq[kk+1]        
            snp = snporder[tback]
            MagicReconstruct.calsitedataprob_singlephase!(dataprobls, chrfhaplo[snp,:], 
                chroffgeno[snp,:],popmakeup;
                epsf=epsfls[snp], epso=epsols[snp], epso_perind, seqerror=seqerrorls[snp], 
                allelebalancemean=allelebalancemeanls[snp], allelebalancedisperse=allelebalancedispersels[snp],
                alleledropout=alleledropoutls[snp], israndallele, issiteGT = issnpGT[snp])            
            calnextlogbackward!(logbwprob, tback, tnow, dataprobls,popmakeup, priorprocess)        
            write(file, string("t",tseq[kk]),logbwprob)
        end
    end
    decodetempfile
end

function calinitforward(dataprobls::AbstractVector,popmakeup::AbstractDict)
    noff = length(dataprobls)    
    fwprob = [similar(i) for i in dataprobls]
    fwlogl = zeros(noff)
    for popid in keys(popmakeup)        
        initprob = popmakeup[popid]["initprob"]        
        offls = popmakeup[popid]["offspring"]        
        for off in offls
            fwprob[off] .= initprob .* dataprobls[off]
            scale = sum(fwprob[off])
            fwprob[off] ./= scale
            fwlogl[off] = log(scale)
        end
    end
    [fwprob,fwlogl]
end

# input fwprob_logl refers to tnow, and output fwprob_logl refers to tnext
# dataprobls refers to tnext
function calnextforward!(fwprob_logl::AbstractVector,tnow::Integer, 
    dataprobls::AbstractVector, popmakeup::AbstractDict,priorprocess::AbstractDict;
    ttrandis::Union{Nothing,Real}=nothing,
	ttran::Real=tnow)
    fwprob, fwlogl = fwprob_logl    
    for popid in keys(popmakeup)
        hashcode = popmakeup[popid]["hashcode"]
        if isnothing(ttrandis)
            tranprob = priorprocess[hashcode].tranprobseq[ttran]
        else
            tranprob = MagicReconstruct.get_tranprobmatrix(ttrandis,priorprocess[hashcode].tranrate)            
        end        
        offls = popmakeup[popid]["offspring"]        
        for off in offls
            dataprob = dataprobls[off]
            if isnothing(tranprob)
                fwprob[off] .= fwprob[off] .* dataprob
            else
                fwprob[off] .= (tranprob' * fwprob[off]) .* dataprob
            end
            scale = sum(fwprob[off])
            fwprob[off] ./= scale
            fwlogl[off] = fwlogl[off] + log(scale)
        end
    end
    fwprob_logl
end


# input logbwprob refers to tback, and output logbwprob refers to tnow
# dataprobls refers to tback not tnow!
# tnow < tback backwardly
function calnextlogbackward!(logbwprob::AbstractVector,tback::Integer, tnow::Integer,
    dataprobls::AbstractVector, popmakeup::AbstractDict,priorprocess::AbstractDict)    
    tnow < tback || @error string("tnow=", tnow, ", tback=",tback, ", not statisfy tnow < tback for backward")
    for popid in keys(popmakeup)
        hashcode = popmakeup[popid]["hashcode"]
        tranprob = priorprocess[hashcode].tranprobseq[tnow]        
        offls = popmakeup[popid]["offspring"]        
        for off in offls
            logpmax = maximum(logbwprob[off])
            logbwprob[off]  .= exp.(logbwprob[off] .- logpmax)
            logbwprob[off] .*= dataprobls[off] 
            if !isnothing(tranprob)
                logbwprob[off] .= tranprob * logbwprob[off]
            end
            logbwprob[off] .= log.(logbwprob[off]) .+ logpmax
        end
    end
    logbwprob
end


function callogbackward!(decodetempfile::AbstractString,chrfhaplo::AbstractMatrix,chroffgeno::AbstractMatrix,
    popmakeup::AbstractDict,priorprocess::AbstractDict;
    epsf::Union{Real,AbstractVector},    
    epso::Union{Real,AbstractVector},
    epso_perind::Union{Nothing,AbstractVector},     
    seqerror::Union{Real,AbstractVector},
    allelebalancemean::Union{Real,AbstractVector},
    allelebalancedisperse::Union{Real,AbstractVector},
    alleledropout::Union{Real,AbstractVector},
    issnpGT::AbstractVector,
    israndallele::Bool, 
    snporder::Union{Nothing,AbstractVector}=nothing)
    hmmdecode_chr(chrfhaplo,chroffgeno,
        popmakeup,priorprocess;
        epsf, epso,epso_perind, seqerror,allelebalancemean, allelebalancedisperse,alleledropout,
        hmmalg="logbackward", decodetempfile = decodetempfile,
        israndallele, issnpGT, snporder)    
    tomarker_major_order!(decodetempfile)
    decodetempfile
end


function tomarker_major_order!(decodetempfile::AbstractString)
    hmmalg, snporder, decode = read_decodefile(decodetempfile)	
    noffspring = length(decode)
    tlength, nstate = size(decode[1][2])     
    nzcolls = first.(decode)
    decodels = last.(decode)
    assparse = issparse(decode[1][2])
    jldopen(decodetempfile, "w") do file        
        write(file, "hmmalg", hmmalg)
        write(file, "noffspring", noffspring)
        write(file, "snporder", snporder)
        write(file, "tlength", tlength)
        write(file, "nstate", nstate)
        write(file, "nzcol", nzcolls)
        for t in 1:tlength
            res = assparse ? spzeros(noffspring, nstate) : zeros(noffspring, nstate)
            for off in 1:noffspring
                nzcol = nzcolls[off]    
                res[off, nzcol] .= decodels[off][t, nzcol]
            end
            write(file,string("t", t),res)
        end
    end
	decodetempfile
end

function tomarker_major_order2!(decodetempfile::AbstractString)    
    tempid = tempname(dirname(decodetempfile),cleanup=true)
    outtempfile = string(tempid, "_", basename(decodetempfile))
    try 
        jldopen(decodetempfile, "r") do file
            hmmalg = file["hmmalg"]
            if !in(hmmalg,["forwardbackward","logforward","logbackward"]) 
                error(string("tomarker_major_order! does not work for hmmalg=",hmmalg))
            end
            snporder = file["snporder"]
            noffspring = file["noffspring"]
            tlength = file["tlength"]
            nstate = file["nstate"]
            jldopen(outtempfile, "w") do outfile
                write(outfile, "hmmalg", hmmalg)
                write(outfile, "snporder", snporder)
                write(outfile, "noffspring", noffspring)
                write(outfile, "tlength", tlength)
                write(outfile, "nstate", nstate)
                for off in 1:noffspring
                    nzcol, off_decode = file[string(off)]
                    tlength2,nstate2 = size(off_decode)
                    tlength == tlength2 || @error string("ïnconsistent #markers")            
                    nstate == nstate2 || @error string("ïnconsistent #states")            
                    off_decode2 = off_decode'
                    for t in 1:tlength
                        outfile[string("t",t, "/offspring",off)] = off_decode2[:,t]
                    end
                    outfile[string("nzcol/offspring",off)] = nzcol
                end
            end
        end
        jldopen(outtempfile, "r") do file
            hmmalg = file["hmmalg"]
            snporder = file["snporder"]
            noffspring = file["noffspring"]
            tlength = file["tlength"]
            nstate = file["nstate"]
            nzcol = file["nzcol"]
            nzcolls = [nzcol[string("offspring",off)] for off in 1:noffspring] 
            jldopen(decodetempfile, "w") do outfile
                write(outfile, "hmmalg", hmmalg)
                write(outfile, "snporder", snporder)
                write(outfile, "noffspring", noffspring)
                write(outfile, "tlength", tlength)
                write(outfile, "nstate", nstate)        
                write(outfile, "nzcol", nzcolls)        
                for t in 1:tlength
                    res = file[string("t",t)]
                    res2 = permutedims(reduce(hcat,[res[string("offspring",off)] for off in 1:noffspring]))
                    write(outfile,string("t",t),res2)
                end
            end
        end
    finally
        rm(outtempfile,force=true)
    end
    decodetempfile
end

