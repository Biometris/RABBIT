using MagicBase
using MagicImpute
using Test
cd(@__DIR__)

# genotype impute
designcode = "4ril-self3"
dataid = designcode*"_gbs"
genofile = dataid*"_magicsimulate_geno.vcf.gz"
magicgeno = magicimpute(genofile, designcode;
    outstem = nothing,
    isinfererror = true, 
    verbose=false
)

# test results
truegeno = formmagicgeno(string(dataid,"_magicsimulate_truegeno.csv.gz"),designcode)
acc = magicaccuracy!(truegeno, magicgeno)
# println(acc)
@test acc["founder"].genoerr < 0.05
@test acc["offspring"].genoerr < 0.05
