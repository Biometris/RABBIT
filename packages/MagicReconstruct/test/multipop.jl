using MagicBase
using MagicReconstruct
using Test
cd(@__DIR__)

# haplotype reconstruct
dataid = "multipop"
genofile = dataid*"_magicsimulate_geno.vcf.gz"
pedfile = dataid*"_magicsimulate_ped.csv"
magicancestry = magicreconstruct(genofile,pedfile;
    # hmmalg = "viterbi",
    outstem=nothing,
    verbose=false
);

# test results
truefgl = formmagicgeno(dataid*"_magicsimulate_truefgl.csv.gz", pedfile)
acc = magicaccuracy!(truefgl, magicancestry)
# println(acc)
@test acc["offspring"].assignerr < 0.05
@test acc["offspring"].callerr <= acc["offspring"].assignerr

#  clean up
# rm.(filter(x->occursin(oustem,x), readdir()))