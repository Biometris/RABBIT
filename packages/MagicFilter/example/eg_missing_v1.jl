
using Revise
# using Pkg
# Pkg.activate(joinpath(@__DIR__, ".."))
using MagicFilter
using MagicBase

cd(@__DIR__)
pwd()

genofile = "sim_geno.vcf.gz"
pedfile = "sim_ped.csv"
@time magicgeno = filter_missing(genofile, pedfile;
    outstem = replace(genofile,".vcf.gz"=>""),
    max_offspringmiss = 0.8,
    max_foundermiss = 0.8,
)
