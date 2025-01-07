using MagicBase
using MagicImpute
using Test
cd(@__DIR__)

# genotype impute
designcode = "4ril-self3"
dataid = designcode*"_mix"
genofile = dataid*"_magicsimulate_geno.vcf.gz"
magicgeno = magicimpute(genofile, designcode;
    outstem = nothing,
    isinfererror = false,
    verbose=false
)

# test results
# magicgeno = formmagicgeno(dataid*"_magicimpute_geno.vcf.gz",designcode)
truegeno = formmagicgeno(string(dataid,"_magicsimulate_truegeno.csv.gz"),designcode)
acc = magicaccuracy!(truegeno, magicgeno)

# println(acc)
@test acc["founder"].genoerr < 0.05
@test acc["offspring"].genoerr < 0.05

#  clean up
# rm.(filter(x->occursin(outstem*"_magicimpute",x), readdir()))
