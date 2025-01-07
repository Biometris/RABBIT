
struct TraitQtl{T<:AbstractFloat}        
    qtl::Vector{Vector{T}}
    founder_haplo::Vector{BitMatrix} # false=allele1, true = allele2
    allele2_effect::Vector{Vector{T}}
    dominance::T
    function TraitQtl(qtl::Vector{Vector{T}}, founder_haplo::Vector{BitMatrix},
        allele2_effect::Vector{Vector{T}},dominance::T) where T<:AbstractFloat
        # TODO: check consistency
        new{T}(qtl,founder_haplo,allele2_effect,dominance)
    end
end

function TraitQtl(nf::Integer, chrlen::Vector,nqtl::Integer,dominance::Real)    
    p = chrlen ./ sum(chrlen)
    nqtl_ls = rand(Multinomial(nqtl,p))
    qtl = map((x,y)->sort(round.(rand(Uniform(0,x),y),digits=4)),chrlen,nqtl_ls)
    founder_haplo = [falses(i,nf) for i in nqtl_ls]
    for chr in eachindex(nqtl_ls)
        for i in 1:nqtl_ls[chr]
            pos = sample(1:nf,rand(DiscreteUniform(1,nf-1)) ,replace=false)
            founder_haplo[chr][i,pos] .= true
        end
    end
    allele2_effect = [round.(rand(Normal(),i),digits=4) for i in nqtl_ls]
    dominance2 = convert(eltype(first(qtl)),dominance)
    TraitQtl(qtl,founder_haplo,allele2_effect,dominance2)
end

function trait_genetic_val(zygote::AbstractVector, traitqtl::TraitQtl)
    g = [begin 
        homolog = [gridhomologfgl(i,traitqtl.qtl[chr]) for i in zygote[chr]]
        fgl = vcat.(homolog...)
        fhaplo = traitqtl.founder_haplo[chr]
        dose = [sum(fhaplo[i,fgl[i]]) for i in eachindex(fgl)]
        eff = traitqtl.allele2_effect[chr]
        a = dot(eff, dose)
        d = traitqtl.dominance*sum(abs.(eff[dose .== 1]))
        a+d 
    end for chr in eachindex(zygote)]
    sum(g)
end

