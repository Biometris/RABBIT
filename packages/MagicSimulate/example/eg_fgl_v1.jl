using Revise
using MagicSimulate
using MagicBase
using Distributions
cd(@__DIR__)
pwd()

using Plots

pedinfo = "multipop_designcode.csv"

# pedinfo = designcode
@time contfgl,magicped = magicsimulate(pedinfo;    
    chrlen = 100*ones(5),
    outstem=nothing
)
# rm.(filter(x->occursin("sim",x), readdir()))


fhaplofile = "sim_fhaplo.vcf.gz"
founderhaplo = MagicBase.formfhaplo(fhaplofile);
