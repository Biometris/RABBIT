
function tocontfgl(magicancestry::MagicAncestry)
    if isnothing(magicancestry.viterbipath)
        @error string("no viterbi path")
        return nothing
    end
    chrlen = [i[end,:poscm]-i[1,:poscm] for i in magicancestry.markermap]
    chrid = [i[1,:linkagegroup] for i in magicancestry.markermap]
    isfounderinbred = get_isfounderinbred(magicancestry)
    markermap = contfgl_markermap(chrlen,chrid)
    nfounder = size(magicancestry.magicped.founderinfo,1)
    fgeno = contfgl_foundergeno(chrlen,isfounderinbred,nfounder)
    offgeno = get_off_contfgl(magicancestry)
    misc = Dict{String, DataFrame}()
    contfgl  = MagicGeno(magicancestry.magicped,markermap,fgeno,offgeno,misc)
    contfgl
end

function contfgl_markermap(chrlen::AbstractVector,chrid::AbstractVector)    
    length(chrlen) == length(chrid) || @error string("mismatch lengths between chrlen and chrid")
    [DataFrame(
        "marker"=>["maternal_homolog_chr"*string(chr),"paternal_homolog_chr"*string(chr)],
        "linkagegroup"=>chrid[chr],
        "poscm"=>[0.0, chrlen[chr]],
        "physchrom"=>missing,
        "physposbp"=>missing,
        "info"=>missing,
        "founderformat" => "contfgl",
        "offspringformat" => "contfgl",
        "foundererror" => missing,
        "offspringerror" => missing,
        "seqerror" => missing,
        "allelebalancemean" => missing,
        "allelebalancedisperse" => missing,
        "alleledropout" => missing,
    ) for chr in eachindex(chrlen)]
end

function contfgl_foundergeno(chrlen::AbstractVector,isfounderinbred::Bool,nfounder::Integer)
    if isfounderinbred
        fgls= [i for i = 1:nfounder, j = 1:2]
    else
        fgls= [2*(i-1)+j for i = 1:nfounder, j = 1:2]
    end
    [[[[0.0,fgls[f,o]],[chrlen[chr],NaN]] for o in 1:2, f in 1:nfounder] 
        for chr in eachindex(chrlen)]
end


function get_isfounderinbred(magicancestry::MagicAncestry)
    nfounder = size(magicancestry.magicped.founderinfo,1)
    fgls = magicancestry.statespace["haplotype"]
    length(fgls) in [nfounder,2*nfounder] || @error string("fgls=",fgls, " mismatch nfounder=",nfounder)
    length(fgls) == nfounder
end

function get_off_contfgl(magicancestry::MagicAncestry)
    is_path_diplo, path_fgls = get_path_fgls(magicancestry)
    [begin 
        poscm = magicancestry.markermap[chr][!,:poscm]
        poscm .-= first(poscm) # set positon of 1st marker to 0
        viterbi = magicancestry.viterbipath[chr]
        contfgl = [begin 
            jpath = MagicBase.tojumppath(dpath)
            dbreaks = jpath[1,:]
            breaks = poscm[dbreaks[1:end-1]]
            push!(breaks,last(poscm))
            states = jpath[2,1:end-1]
            if is_path_diplo    
                fgls = vcat(path_fgls[states],[NaN,NaN])
                homolog1 = map((x,y)->[x,y],breaks,first.(fgls))
                homolog2 = map((x,y)->[x,y],breaks,last.(fgls))
                shrink_homolog!(homolog1)
                shrink_homolog!(homolog2)
            else
                fgls = vcat(states,[NaN])
                homolog1 = map((x,y)->[x,y],breaks,fgls)
                homolog2 = map((x,y)->[x,y],breaks,fgls)
            end
            [homolog1,homolog2]
        end for dpath in eachcol(viterbi)]
        reduce(hcat, contfgl)
    end for chr in eachindex(magicancestry.markermap)]
end

function get_path_fgls(magicancestry::MagicAncestry)
    haploindex = magicancestry.statespace["haploindex"]
    nfgl = length(haploindex)
    is_path_diplo = false
    for viterbi in magicancestry.viterbipath
        if any(unique(viterbi) .> nfgl)
            is_path_diplo = true
            break
        end    
    end
    path_fgls = is_path_diplo ? MagicBase.prior_diploindex(nfgl) : haploindex 
    is_path_diplo, path_fgls
end


function shrink_homolog!(homolog::AbstractVector)
    n = length(homolog)
    b = falses(n)
    for i in 2:n-1
        b[i] = homolog[i][2] == homolog[i-1][2]
    end
    deleteat!(homolog,b)
end
