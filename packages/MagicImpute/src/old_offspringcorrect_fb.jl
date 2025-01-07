
function offspringcorrect!(magicgeno::MagicGeno;
    model::AbstractString="jointmodel",
    detectthreshold::Real=0.9,
    epsf::Real=0.005,
    epso::Real=0.005,    
    epso_perind::Union{Nothing,AbstractVector}, 
    isfounderinbred::Bool=true,
    israndallele::Bool=true, 
    isfounderinbred::Bool=true,        
    snporder::Union{Nothing,AbstractVector}=nothing,
    isparallel::Bool=true,
    io::Union{Nothing,IO}=nothing,
    verbose::Bool=true)
    isnothing(io) && (io=IOBuffer(append=true))
    magicprior=MagicReconstruct.calmagicprior(magicgeno,model; isfounderinbred)
    nchr = length(magicgeno.markermap)
    jltempdir = mktempdir(tempdirectory; prefix="jl_MagicImpute_", cleanup=true)
    tempid = tempname(jltempdir,cleanup=false)
    refinetempfilels = [string(tempid,"_offspringcorrect_chr",chr,".jld2")
        for chr in 1:nchr]
    try
        if isparallel && nprocs()>1
            magicgenols =[submagicgeno(magicgeno,chrsubset=[chr])
                for chr=1:nchr]
            resls = pmap((x,y)->offspringcorrect_chr!(x, 1, magicprior;
                model, detectthreshold, epsf,epso,epso_perind, 
                isfounderinbred,issnpGT,snporder,
                refinetempfile=y,io, verbose), magicgenols,refinetempfilels)
            for chr=1:nchr
                iobuffer = resls[chr][2]
                write(io,String(take!(iobuffer)))
                flush(io)
                close(iobuffer)
            end
        else
            for chr =1:nchr
                offspringcorrect_chr!(magicgeno, chr, magicprior;
                    model, detectthreshold, epsf,epso,epso_perind, 
                    isfounderinbred,issnpGT, snporder,
                    refinetempfile=refinetempfilels[chr],io, verbose)
            end
        end
    finally
        rm.(refinetempfilels;force=true)
        rm(jltempdir;force=true,recursive=true)
    end
    magicgeno
end

function offspringcorrect_chr!(magicgeno::MagicGeno,chr::Integer,
    magicprior::NamedTuple;
    model::AbstractString="jointmodel",
    detectthreshold::Real=0.9,
    epsf::Real=0.005,
    epso::Real=0.005,
    epso_perind::Union{Nothing,AbstractVector}, 
    isfounderinbred::Bool=true,
    israndallele::Bool=true, 
    isfounderinbred::Bool=true,        
    snporder::Union{Nothing,AbstractVector}=nothing,
    refinetempfile::AbstractString=tempdir(),
    io::Union{Nothing,IO}=nothing,
    verbose::Bool=true)
    isnothing(io) && (io=IOBuffer(append=true))
    popmakeup,priorprocess = MagicReconstruct.calpriorprocess(magicgeno,chr,
        model,magicprior; isfounderinbred)
    chrfoundergeno = magicgeno.foundergeno[chr]
    chroffgeno = magicgeno.offspringgeno[chr]
    ischrx = lowercase(magicgeno.markermap[chr][1,:linkagegroup]) == "chrx"
    itmax = 20
    chrid = magicgeno.markermap[chr][1,:linkagegroup]
    for it=1:itmax
        errorpos= offspringcorrect_chr(chrfoundergeno,
            chroffgeno,  popmakeup, priorprocess;
            epsf,epso, epso_perind, detectthreshold,isfounderinbred,
            issnpGT,ischrx,snporder,refinetempfile)
        nerror = length(errorpos)
        errratio = round(Int,1000*nerror/*(size(chroffgeno)...))
        if nerror >  0
            # TODO: rewrite updating of error correcting
            missinggeno = any(issnpGT) ? [ones(4)/4 for i in 1:length(errorpos)] : "NN"
            chroffgeno[errorpos] .= missinggeno
        end
        msg = string("chr=", chrid, ", it=",it,
            ", #error =", nerror," (", errratio/10.0, "%)")
        printconsole(io,verbose,msg)
        nerror == 0 && break
    end
    io
end

function offspringcorrect_chr(chrfoundergeno::AbstractMatrix,chroffgeno::AbstractMatrix,
    popmakeup::AbstractDict,priorprocess::AbstractDict;
    epsf::Union{Real,AbstractVector}, 
    epso::Union{Real,AbstractVector},
    epso_perind::Union{Nothing,AbstractVector}, 
    callthreshold::Real,
    isfounderinbred::Bool,
    israndallele::Bool=true, 
    isfounderinbred::Bool=true,        
    ischrx::Bool,
    snporder::Union{Nothing,AbstractVector}=nothing,
    refinetempfile::AbstractString)
    isnothing(snporder) && (snporder = 1:size(chroffgeno,1))
    _, chrcondprob = MagicReconstruct.hmmdecode_chr(chrfoundergeno,chroffgeno,
        popmakeup,priorprocess;
        epsf,epso, epso_perind, hmmalg="forwardbackward",decodetempfile = refinetempfile,
        isfounderinbred,issnpGT,snporder)
    tindex = sortperm(snporder) # tindex[i]: t index for inputsnp index i
    chrcondprob = [i[tindex,:] for i in chrcondprob]
    ishaploidls = [i["ishaploid"] for i in values(popmakeup)]
    ishaploprob = all(ishaploidls)
    snpgenoprob = MagicImpute.calsnpgenoprob(chrfoundergeno,chrcondprob,
        ishaploprob, ischrx, epsf)
    # calculuate called offpsring genotypes
    chrcalledinput = copy(chroffgeno)
    if !all(issnpGT)
        readcallthreshold = callthreshold
        chrcalledinput[.!issnpGT,:] .= MagicBase.callfromprob.(chroffgeno[.!issnpGT,:],
            readcallthreshold; isphased=false)
    end
    errpos, newcalledgeno = offspringcorrect_chr(snpgenoprob,chrcalledinput; callthreshold)
    errpos, newcalledgeno
end

function offspringcorrect_chr(snpgenoprob::Matrix{Vector{T}} where T<:AbstractFloat,
    chrcalledinput::AbstractMatrix;  callthreshold::AbstractFloat)
    rule = Dict(["11","12","22","1N","2N","NN"] .=> [[1],[2,3],[4],[1,2,3],[2,3,4],[1,2,3,4]])
    obsii = [get(rule,i,[]) for i in chrcalledinput]
    obsprob = map((x,y)->sum(x[y]),snpgenoprob, obsii)
    errpos = findall(obsprob .< 1-callthreshold)
    # callthreshold2 = 0.55
    newcalledgeno = MagicBase.callfromprob.(snpgenoprob[errpos],callthreshold;isphased=false)
    errpos, newcalledgeno
end

function offspringcorrect_chr(chrestgeno::Matrix{T} where T<:AbstractFloat, chrcalledinput::AbstractMatrix)
    errors = Vector{CartesianIndex{2}}()
    pos = findall(chrestgeno .!= chrcalledinput)
    pos = pos[chrcalledinput[pos] .!= "NN"]
    pos = pos[chrestgeno[pos] .!= "NN"]
    posxx = pos[chrestgeno[pos] .== "11"]
    if !isempty(posxx)
        push!(errors, posxx[chrcalledinput[posxx] .== "12"]...)
        push!(errors, posxx[chrcalledinput[posxx] .== "21"]...)
        push!(errors, posxx[chrcalledinput[posxx] .== "22"]...)
        push!(errors, posxx[chrcalledinput[posxx] .== "2N"]...)
        push!(errors, posxx[chrcalledinput[posxx] .== "N2"]...)
        pos = setdiff(pos, posxx)
    end
    posxx = pos[chrestgeno[pos] .== "22"]
    if !isempty(posxx)
        push!(errors, posxx[chrcalledinput[posxx] .== "12"]...)
        push!(errors, posxx[chrcalledinput[posxx] .== "21"]...)
        push!(errors, posxx[chrcalledinput[posxx] .== "11"]...)
        push!(errors, posxx[chrcalledinput[posxx] .== "1N"]...)
        push!(errors, posxx[chrcalledinput[posxx] .== "N1"]...)
        pos = setdiff(pos, posxx)
    end
    posxx = pos[chrestgeno[pos] .== "12"]
    in("21", unique(chrestgeno[pos])) && @error("unexpected 21")
    if !isempty(posxx)
        push!(errors, posxx[chrcalledinput[posxx] .== "11"]...)
        push!(errors, posxx[chrcalledinput[posxx] .== "22"]...)
        pos = setdiff(pos, posxx)
    end
    posxx = pos[chrestgeno[pos] .== "1N"]
    in("N1", unique(chrestgeno[pos])) && @error("unexpected N1")
    if !isempty(posxx)
        push!(errors, posxx[chrcalledinput[posxx] .== "22"]...)
        pos = setdiff(pos, posxx)
    end
    posxx = pos[chrestgeno[pos] .== "2N"]
    in("N2", unique(chrestgeno[pos])) && @error("unexpected N2")
    if !isempty(posxx)
        push!(errors, posxx[chrcalledinput[posxx] .== "11"]...)
        pos = setdiff(pos, posxx)
    end
    isempty(pos) || @error string("unexpected ", unique(chrestgeno[pos]))
    newestgeno = chrestgeno[errors]
    errors, newestgeno
end
