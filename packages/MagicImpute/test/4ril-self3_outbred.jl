using MagicBase
using MagicImpute
using Test
cd(@__DIR__)


# genotype impute
isfounderinbred = false
designcode = "4ril-self3"
dataid = designcode*"_outbred"
genofile = dataid*"_magicsimulate_geno.vcf.gz"
magicgeno = magicimpute(genofile,designcode;
    isfounderinbred,
    isdelmarker = true, 
    iscorrectfounder = false, 
    isinfererror = false,
    verbose=false,
    outstem = nothing, 
);

truegeno = formmagicgeno(string(dataid,"_magicsimulate_truegeno.csv.gz"),designcode;isfounderinbred);
acc = magicaccuracy!(truegeno, magicgeno)

@test acc["founder"].genoerr < 0.07
@test acc["offspring"].genoerr < 0.05

#  clean up
# rm.(filter(x->occursin(dataid*"_magicimpute",x), readdir()))
