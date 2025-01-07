using MagicBase
using MagicReconstruct 
using Test
cd(@__DIR__)

# haplotype reconstruct
designcode = "4ril-self3"
dataid = designcode*"_array"
genofile = dataid*"_magicsimulate_geno.vcf.gz"
magicancestry = magicreconstruct(genofile,designcode;
    outstem=nothing,
    verbose=false
);

# test results
truefgl = formmagicgeno(string(dataid,"_magicsimulate_truefgl.csv.gz"),designcode)
acc = magicaccuracy!(truefgl, magicancestry)
# println(acc)
@test acc["offspring"].assignerr < 0.05
@test acc["offspring"].callerr < 0.05

#  clean up
# rm.(filter(x->occursin(oustem,x), readdir()))
