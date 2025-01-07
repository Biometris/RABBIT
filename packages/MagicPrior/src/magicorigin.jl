
"""
    magicorigin(pedigree::Pedigree; kwargs...)

calculate expected ancestral probabilities and junction densities  in pedigree.

# Keyword Argument

`memberlist::Union{Nothing,AbstractVector}=nothing`: a list of pedigree members.
By default, memberlist contains the last pedigree member.

`isautosome::Bool=true`: autosome rather than sex chromosome,

`isconcise::Bool=false`: if false, calculate quantities: phi12, R(a)(or R^m), R(b)(or R^p),
ρ, J1112, J1121, J1122, J1211, J1213, J1222, J1232. If true, calculate
only the first 4 quantities.

`isfglexch::Bool=false`: if true, assume that  founder genotype lables (FGLs) are exchangeable.

"""
function magicorigin(pedigree::Pedigree;
    memberlist::Union{Nothing,AbstractVector}=nothing,
    isfounderinbred::Bool=true,
    isautosome::Bool=true,
    isfglexch::Bool=false,
    isconcise::Bool=false)
    # `founderfgl::Union{Nothing,AbstractMatrix}=nothing`: founderfgl[i,j] gives the
    # matenally(j=1, denoted by ^m) or patenally (j=2, denoted by ^p) derived
    # halploid genome label in founder i. if founderfgl is nothing, assign a unique
    # label to diploid genome of each founder.
    if isnothing(memberlist)
        memberlist=[pedigree.member[end]]
    else
        b = [!in(i,pedigree.member) for i in memberlist]
        if any(b)
            @error string(memberlist[b], " in memberlist are not members of pedigree")
        end
    end
    founderfgl = get_founderfgl(pedigree.nfounder;isfounderinbred)    
    if isfglexch
        res=pedidentity(pedigree,founderfgl,memberlist,isautosome,isconcise)
    else
        res=pedancestry(pedigree,founderfgl,memberlist,isautosome,isconcise)
    end    
    dememoize()
    res
end
