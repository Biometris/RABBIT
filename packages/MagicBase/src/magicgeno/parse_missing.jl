
function getismissing(genovec::AbstractVector, format::AbstractString)
    missval = MagicBase.get_missingcode(format)
    [i == missval for i in genovec]
end

function getismissing(genomtx::AbstractMatrix, formatvec::AbstractVector)
    formatset = unique(formatvec)
    ismiss = BitMatrix(undef, size(genomtx)...)
    for format in formatset
        ii = findall(formatvec .== format)
        missval = MagicBase.get_missingcode(format)
        ismiss[ii,:] .= [i == missval for i in genomtx[ii,:]]
    end
    ismiss
end

function count_missing(genomtx::AbstractMatrix, formatvec::AbstractVector;
    permarker::Bool=true)
    ismiss = getismissing(genomtx,formatvec)
    if permarker
        sum(ismiss,dims=2)[:,1]
    else
        sum(ismiss,dims=1)[1,:]
    end
end

# used for info_missing
function count_missing(magicgeno::MagicGeno; targetfounder::Bool=true)
    nchr = length(magicgeno.markermap)
    nmissls = zeros(nchr)
    ntotalls = zeros(nchr)
    for chr in 1:nchr
        if targetfounder
            genomtx = magicgeno.foundergeno[chr]
            formatvec = magicgeno.markermap[chr][!,:founderformat]
        else
            genomtx = magicgeno.offspringgeno[chr]
            formatvec = magicgeno.markermap[chr][!,:offspringformat]
        end
        nmissls[chr] = sum(getismissing(genomtx, formatvec))
        ntotalls[chr] = length(genomtx)
    end
    chridls = [isempty(i) ? "missing" : i[1,:linkagegroup] for i in magicgeno.markermap]
    fracls = nmissls ./ ntotalls
    DataFrame(chromosome=chridls,nmiss = nmissls, ntotal=ntotalls, missfrac=fracls)
end

function get_missingcode(formatls::AbstractVector)
    dict = Dict([i=>get_missingcode.(i) for i in unique(formatls)])
    [dict[i] for i in formatls]
end

function get_missingcode(format::AbstractString)
    if format=="GT_haplo"
        "N"
    elseif format=="GT_unphased"
        # half-missings such as "1N" and "2N" are regarded as non-missing
        "NN"
    elseif format == "GT_phased"
        ["N","N"]
    elseif format == "AD"
        [0,0]
    elseif format == "GL"
        [0.0,0.0,0.0]
    elseif format == "GP"
        [0.25,0.5,0.25]
    else
        @error string("TODO for format =",format)
    end
end

function info_missing(magicgeno::MagicGeno;
    io::Union{IO,Nothing}=nothing, verbose::Bool=true)
    res = count_missing(magicgeno;targetfounder=true)
    fmiss = sum(res[!,:nmiss])/sum(res[!,:ntotal])
    msg = string("founder frac_missing = ", round(fmiss,digits=5))
    res = count_missing(magicgeno;targetfounder=false)
    fmiss = sum(res[!,:nmiss])/sum(res[!,:ntotal])
    msg *= string(", offspring frac_missing = ", round(fmiss,digits=5))
    MagicBase.printconsole(io,verbose,msg)
end

function merge_missprogeny!(magicgeno::MagicGeno,missgeno::MagicGeno)    
    # dict for chr
    chrls = [i[1,:linkagegroup] for i in magicgeno.markermap]
    misschrls = [isempty(i) ? nothing : i[1,:linkagegroup] for i in missgeno.markermap]
    deleteat!(misschrls,isnothing.(misschrls))
    if !issubset(string.(misschrls),string.(chrls))
        msg = string("misschrls =",misschrls, ", is not a sbuset of chrls=",chrls)
        @warn msg
    end
    misschrdict = Dict(misschrls .=> 1:length(misschrls))
    chr2misschr = [get(misschrdict,chr,nothing) for chr in chrls]
    # calc inputmap and initialize
    inputmap = groupby(missgeno.misc["inputmarkermap"],:linkagegroup)    
    is_miss_bp = any([in("", i[!,:physposbp]) || any(ismissing.(i[!,:physposbp])) for i in inputmap])
    poscol = is_miss_bp ? :poscm : :physposbp
    # @info string("merge_missprogeny! according to poscol=",poscol)
    nf = size(magicgeno.magicped.founderinfo,1)
    noff = size(magicgeno.magicped.offspringinfo,1)
    for chr in eachindex(chrls)
        isnothing(chr2misschr[chr]) && continue
        misschr = chr2misschr[chr]
        addmap = missgeno.markermap[misschr]
        isempty(addmap) && continue        
        if poscol == :physposbp  
            addmap[!,:poscm] .= missing
        end        
        nowmap = magicgeno.markermap[chr]      
        newmap = vcat(nowmap,addmap)
        # restore input marker ordering
        # sort!(newmap,poscol)         
        newdict = Dict(newmap[!,:marker] .=> 1:size(newmap,1))        
        inputindices = [get(newdict,i,missing) for i in inputmap[chr][!,:marker]]        
        deleteat!(inputindices,ismissing.(inputindices))
        newmap = newmap[inputindices, :]
        # interpolate poscm
        b = ismissing.(newmap[!,:poscm])
        if poscol == :physposbp            
            interp_mono = extrapolate(interpolate(newmap[.!b,:physposbp],Float64.(newmap[.!b,:poscm]),  LinearMonotonicInterpolation()), Flat());
            newmap[b,:poscm] .= [interp_mono(i) for i in newmap[b,:physposbp]]  
        else            
            any(b) && @warn string("Could not iterpolate poscm for missgeno!")
            issorted(newmap,poscol) || @error string(poscol, " is not sorted")
        end      
        # save
        newnsnp = size(newmap,1)
        newdict = Dict(newmap[!,:marker] .=> 1:newnsnp)        
        nowindices = [newdict[i] for i in nowmap[!,:marker]]
        addindices = [newdict[i] for i in addmap[!,:marker]]
        fgeno = Matrix(undef,newnsnp, nf)
        fgeno[nowindices,:] .= magicgeno.foundergeno[chr]
        fgeno[addindices,:] .= missgeno.foundergeno[misschr]
        offgeno = Matrix(undef,newnsnp, noff)
        offgeno[nowindices,:] .= magicgeno.offspringgeno[chr]
        offgeno[addindices,:] .= missgeno.offspringgeno[misschr]
        # up magicgeno
        magicgeno.markermap[chr] = newmap
        magicgeno.foundergeno[chr] = fgeno
        magicgeno.offspringgeno[chr] = offgeno
    end
    magicgeno
end

function split_missprogeny!(magicgeno::MagicGeno)
    misc = Dict("inputmarkermap"=>reduce(vcat,magicgeno.markermap))
    markermap = Vector()
    fgeno = Vector()
    offgeno = Vector()
    noff = size(magicgeno.magicped.offspringinfo,1)
    for chr in 1:length(magicgeno.markermap)
        chrgeno = magicgeno.offspringgeno[chr]
        formatls = magicgeno.markermap[chr][!,:offspringformat]
        missfrac = MagicBase.count_missing(chrgeno,formatls) ./ noff
        b_del = missfrac .== 1.0
        push!(markermap,magicgeno.markermap[chr][b_del,:])
        push!(fgeno,magicgeno.foundergeno[chr][b_del,:])
        push!(offgeno,magicgeno.offspringgeno[chr][b_del,:])
        if any(b_del)
            b_keep = .!b_del
            magicgeno.markermap[chr] = magicgeno.markermap[chr][b_keep,:]
            magicgeno.foundergeno[chr] = magicgeno.foundergeno[chr][b_keep,:]
            magicgeno.offspringgeno[chr] = magicgeno.offspringgeno[chr][b_keep,:]
        end
    end
    missgeno = MagicGeno(magicgeno.magicped,markermap,fgeno,offgeno,misc)
    missgeno
end

