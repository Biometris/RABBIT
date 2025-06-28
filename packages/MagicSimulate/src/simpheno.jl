

function simpheno(contfgl::MagicGeno,magicfgl::MagicGeno;
    pheno_nqtl::Integer=10,    
    pheno_h2::Real= 0.5)
    fgldf = magicfgl.misc["fgl"]
    fglls = Vector.(eachrow(fgldf[!,2:3]))
    finfo = contfgl.magicped.founderinfo
    finfo[!,1] == fgldf[!,1] || @error "inconsistent founders"
    nfgl = length(unique(reduce(vcat,unique.(fglls))))
    chrlen = [i[end,:poscm] for i in contfgl.markermap]

    # TODO: add dom for sim_trait
    pheno_dominance = 0.0
    sim_trait = TraitQtl(nfgl, chrlen, pheno_nqtl,pheno_dominance)
    nchr = length(contfgl.markermap)
    # fgl_qtl,geno_qtl,dose_qtl
    fgl_qtl = [begin
        homolog = [gridhomologfgl(i,sim_trait.qtl[chr]) for i in contfgl.offspringgeno[chr]]
        fgl = [vcat.(c...) for c in eachcol(homolog)];
        fgl2 = reduce(hcat,fgl)
        isempty(fgl2) ? nothing : fgl2
    end for chr in 1:nchr]
    geno_qtl = [if isnothing(fgl_qtl[chr])
            nothing
        else
            fhaplo = sim_trait.founder_haplo[chr]
            dose =[[fhaplo[i,j] .+ 1 for j in fgl_qtl[chr][i,:]]  for i in 1:size(fhaplo,1)]
            reduce(hcat,dose)'
        end for chr in 1:nchr]
    dose_qtl = [isnothing(i) ? i : sum.(i) .- 2 for i in geno_qtl]
    # gadd_qtl
    gadd_qtl = [if isnothing(dose_qtl[chr])
            nothing
        else
            reduce(hcat,sim_trait.allele2_effect[chr] .* eachrow(dose_qtl[chr]))'
        end for chr in 1:nchr]
    gadd_off = sum(sum.(gadd_qtl[.!isnothing.(gadd_qtl)],dims=1))[1,:]
    vadd = pheno_h2/(1.0-pheno_h2) # set var_noise = 1
    genetic_scale = vadd/StatsBase.std(gadd_off)
    gadd_qtl = [isnothing(i) ? i : round.(i*genetic_scale,digits=5) for i in gadd_qtl];
    # sim phenotype
    offspring = contfgl.magicped.offspringinfo[!,:individual]
    members = contfgl.magicped.offspringinfo[!,:member]
    pheno_gadd = sum(sum.(gadd_qtl[.!isnothing.(gadd_qtl)],dims=1))[1,:]
    pheno_noise = rand(Normal(),length(offspring))
    pheno_noise ./= StatsBase.std(pheno_noise)
    pheno_sim = pheno_gadd + pheno_noise
    pheno_df = DataFrame(:individual=>offspring,:population => members, :phenotype =>pheno_sim)
    # map_qtl
    chridls = [i[1,:linkagegroup] for i in contfgl.markermap]
    founderls = contfgl.magicped.founderinfo[!,:individual]
    map_qtl = get_map_qtl(sim_trait,chridls,founderls,fglls)
    map_qtl[!,:allele2_effect] .*= genetic_scale
    var_qtl = reduce(vcat,var.(gadd_qtl[.!isnothing.(gadd_qtl)],dims=2))[:,1]
    insertcols!(map_qtl, 5,:var_qtl=>var_qtl)
    insertneighbor!(map_qtl,magicfgl)
    # save
    qtlid = map_qtl[!,:marker]
    fgl_qtl2,geno_qtl2,gadd_qtl2 = [begin
        data = hcat(qtlid,join.(reduce(vcat,res[.!isnothing.(res)]),"|"))
        DataFrame(data,vcat(["marker"],offspring))
    end  for res in [fgl_qtl,geno_qtl,gadd_qtl]]
    (fgl_qtl=fgl_qtl2, geno_qtl=geno_qtl2,gadd_qtl=gadd_qtl2,
        map_qtl=map_qtl,pheno_df=pheno_df)
end

function get_map_qtl(sim_trait::TraitQtl,chridls::AbstractVector,
    founderls::AbstractVector,fglls::AbstractVector)
    linkagegroup = reduce(vcat,[repeat([chridls[i]],length(sim_trait.qtl[i])) for i in eachindex(chridls)])
    marker = [string("qtl",i) for i in 1:length(linkagegroup)]
    position = reduce(vcat, sim_trait.qtl)
    allele2_effect = reduce(vcat, sim_trait.allele2_effect)
    df = DataFrame(marker=marker, linkagegroup=linkagegroup,
        poscm=position,allele2_effect=allele2_effect)
    fhaplo = reduce(vcat, sim_trait.founder_haplo) .+ 1
    fgeno = [join.(eachrow(fhaplo[:,fgl]),"|") for fgl in fglls]
    df= hcat(df,DataFrame(fgeno,Symbol.(founderls)))
    df
end


function insertneighbor!(map_qtl::DataFrame, magicfgl::MagicGeno)
    chrls = [i[1,:linkagegroup] for i in magicfgl.markermap]
    nbr = [begin
        chrid, pos = map_qtl[i, [:linkagegroup, :poscm]]
        chr = findfirst(chrls .== chrid)
        posls = magicfgl.markermap[chr][!,:poscm]
        leftnbr = findlast(posls .<= pos)
        leftsnp = isnothing(leftnbr) ? "nothing" : magicfgl.markermap[chr][leftnbr,:marker]
        rightnbr = findfirst(posls .> pos)
        rightsnp = isnothing(rightnbr) ? "nothing" : magicfgl.markermap[chr][rightnbr,:marker]
        [leftsnp,rightsnp]
    end for i in 1:size(map_qtl,1)]
    nbr2 = reduce(hcat,nbr)
    insertcols!(map_qtl,4, :leftsnp=>nbr2[1,:])
    insertcols!(map_qtl,5, :rightsnp=>nbr2[2,:])
    map_qtl
end
