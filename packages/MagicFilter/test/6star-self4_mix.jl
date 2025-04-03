using MagicBase
using MagicFilter
using Test
cd(@__DIR__)

# filtering
dataid = "6star-self4_mix"
genofile = dataid*"_magicsimulate_geno.vcf.gz"
pedfile = dataid*"_magicsimulate_ped.csv"
outstem=dataid
magicgeno = magicfilter(genofile, pedfile;
    minmaf = 0.01,
    missfilter = (f,o) -> f<=0.0 || o <= 0.5,
    outstem,
    verbose=false
);

# test results
ingeno = formmagicgeno(genofile, pedfile)
truegeno = formmagicgeno(string(dataid,"_magicsimulate_truegeno.csv.gz"),pedfile)

acc0 = magicaccuracy!(truegeno, ingeno)
acc = magicaccuracy!(truegeno, magicgeno)
# println(acc0)
# println(acc)
@test acc["founder"].ngenoerr <= acc0["founder"].ngenoerr
@test acc["offspring"].ngenoerr <= acc0["offspring"].ngenoerr

#  clean up
cd(@__DIR__)
dataid = "6star-self4_mix"
outstem=dataid*"_magicfilter"
rm.(filter(x->occursin(outstem,x), readdir()))

