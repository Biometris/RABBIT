# using Revise
using MagicImpute, MagicBase
using Test
cd(@__DIR__)

# genotype impute
dataid = "multipop"
genofile = dataid*"_magicsimulate_geno.vcf.gz"
pedfile = dataid*"_magicsimulate_ped.csv"
magicgeno = magicimpute(genofile, pedfile;
    # iscorrectfounder = true,
    # isdelmarker = true,
    outstem=nothing,
    isinfererror = false,
    verbose=false
);

# test results
truegeno = formmagicgeno(string(dataid,"_magicsimulate_truegeno.csv.gz"),pedfile)
acc = magicaccuracy!(truegeno, magicgeno)

# println(acc)
@test acc["founder"].genoerr < 0.05
@test acc["offspring"].genoerr < 0.05

#  clean up
# rm.(filter(x->occursin("outstem",x), readdir()))
