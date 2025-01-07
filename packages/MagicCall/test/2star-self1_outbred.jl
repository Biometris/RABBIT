using MagicBase
using MagicCall
using Test
cd(@__DIR__)


# single marker genotype call
isfounderinbred = false
designcode = "2star-self1"
dataid = designcode*"_outbred"
genofile = dataid*"_magicsimulate_geno.vcf.gz"

magiccall(genofile,designcode;
    isfounderinbred,    
    verbose=false,
    priorlikeparameters = PriorLikeParameters(seqerror = MagicCall.Beta(1,9)),
    outstem = dataid, 
)

truefile = string(dataid,"_magicsimulate_truegeno.csv.gz")
estfile = string(dataid, "_magiccall_geno.vcf.gz")
acc = magicaccuracy(truefile,estfile,designcode;isfounderinbred)

@test acc["founder"].genoerr < 0.1
@test acc["offspring"].genoerr < 0.07


#  clean up
rm.(filter(x->occursin(dataid*"_magiccall",x), readdir()))

