# using Revise
using MagicBase
using MagicReconstruct 
using Test
cd(@__DIR__)

# haplotype reconstruct
dataid =  "4ril-self3_outbred"
genofile = dataid*"_magicimpute_geno.vcf.gz"
designcode = "4ril-self3"
magicancestry = magicreconstruct(genofile,designcode;
    isfounderinbred = false,
    outstem=nothing,
    verbose=false
)

truefgl = formmagicgeno(string(dataid,"_magicsimulate_truefgl.csv.gz"),
    designcode; isfounderinbred=false);
acc = magicaccuracy!(truefgl, magicancestry)
@test acc["offspring"].assignerr < 0.07
@test acc["offspring"].callerr < 0.05

#  clean up
# rm.(filter(x->occursin(oustem,x), readdir()))
