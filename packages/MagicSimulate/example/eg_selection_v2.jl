
using Revise
using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
using MagicBase, MagicSimulate
using Pedigrees
using DataFrames
using Distributions
using CSV
cd(@__DIR__)
pwd()
workdir = pwd()

select_proportion = 0.1
select_nqtl = 50

fhaplofile = "ArabMAGIC_founderhaplo.csv"
@time obsgeno, truegeno, magicfgl, contfgl= magicsimulate(fhaplofile,"pedinfo.csv"; 
    select_nqtl,select_proportion,    
    outstem=nothing
)


@time magicsimulate(fhaplofile,"ArabMAGIC_ped.csv";
    select_nqtl,select_proportion
)

rm.(filter(x->occursin("outstem",x), readdir()))

propertynames(contfgl)

chrlen= [i[end,:poscm] for i in contfgl.markermap]

pheno_nqtl = 10
nf = size(contfgl.magicped.founderinfo,1)
phenoqtl = TraitQtl(chrlen, pheno_nqtl, nf) 



zygote = [i[:,1] for i in contfgl.offspringgeno]


homolog = [gridhomologfgl(i,traitqtl.qtl[chr]) for i in chrzygote]
fgl = vcat.(homolog...)
fhaplo = traitqtl.founder_haplo[chr]
dose = [sum(fhaplo[i,fgl[i]]) for i in 1:length(fgl)]
sum(traitqtl.allele2_effect[chr] .* dose)


trait_gadd(zygote,phenoqtl)
