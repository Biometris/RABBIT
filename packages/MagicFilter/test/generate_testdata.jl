
using MagicSimulate, MagicBase
cd(@__DIR__)

# simulate data
designcode = "6star-self4"
fhaplofile = "fhaplo.vcf.gz"
ncluster = 2
nsnpchr = 200
simfhaplo(nsnp=ncluster*nsnpchr, nparent=6,
    chrlen=100*ones(ncluster),
    outfile=fhaplofile
)

magicsimulate(fhaplofile,designcode;
    popsize = 200,
    outstem = designcode*"_mix",
    foundermiss = MagicSimulate.Beta(1,2),
    offspringmiss = MagicSimulate.Beta(1,2),
    seqfrac = 0.5,
    seqdepth = MagicSimulate.Gamma(1,2),
)

#  clean up
rm(fhaplofile)
rm.(filter(x->occursin(".log",x), readdir()))
rm.(filter(x->occursin("_ped.png",x), readdir()))
rm.(filter(x->occursin("_truecontfgl",x), readdir()))
rm.(filter(x->occursin("_truefgl",x), readdir()))
