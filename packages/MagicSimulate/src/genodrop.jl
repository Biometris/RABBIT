
function genedrop(magicfgl::MagicGeno; isfounderinbred::Bool=true, isphased::Bool=true)
    markermap = deepcopy(magicfgl.markermap)
    foundergeno = deepcopy(magicfgl.foundergeno)
    nchr = length(foundergeno)
    offspringgeno = [begin
        fgeno = foundergeno[chr]
        if isfounderinbred
            fhaplo = permutedims(fgeno)
            markermap[chr][!,:founderformat] .= "GT_haplo"
        else
            fhaplo = vcat([hcat(col...) for col = eachcol(fgeno)]...)
            markermap[chr][!,:founderformat] .= "GT_phased"
        end        
        offfgl = permutedims(magicfgl.offspringgeno[chr])
        nsnp = size(offfgl,2)
        if isphased
            trueoff = [[fhaplo[i, snp] for i=offfgl[:, snp]] for snp=1:nsnp]
            markermap[chr][!,:offspringformat] .= "GT_phased"
        else
            trueoff = [[join(fhaplo[i,snp]) for i in offfgl[:,snp]] for snp = 1:nsnp]
            markermap[chr][!,:offspringformat] .= "GT_unphased"
        end
        permutedims(reduce(hcat,trueoff))
    end for chr=1:nchr]
    MagicGeno(magicfgl.magicped,markermap,foundergeno,offspringgeno,Dict())
end
