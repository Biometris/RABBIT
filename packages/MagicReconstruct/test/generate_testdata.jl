
using MagicSimulate, MagicBase
cd(@__DIR__)

# simulate data
designcode = "4ril-self3"
fhaplofile = "fhaplo.csv"
ncluster = 1
nsnpchr = 200
simfhaplo(nsnp=ncluster*nsnpchr, nparent=4,
    chrlen=100*ones(ncluster),
    outfile=fhaplofile
)
magicsimulate(fhaplofile,designcode;
    popsize = 2,
    outstem = designcode*"_array",
)
magicsimulate(fhaplofile,designcode;
    popsize = 2,
    outstem = designcode*"_gbs",
    seqfrac = 1.0,
    seqdepth = MagicSimulate.Gamma(2,5),
)
magicsimulate(fhaplofile,designcode;
    popsize = 2,
    outstem = designcode*"_mix",
    seqfrac = 0.5,
    seqdepth = MagicSimulate.Gamma(2,5),
)

#  clean up
rm(fhaplofile)
rm.(filter(x->occursin(".log",x), readdir()))
rm.(filter(x->occursin("_ped",x) && occursin(designcode,x), readdir()))
rm.(filter(x->occursin("_truecontfgl",x), readdir()))
rm.(filter(x->occursin("_truegeno",x), readdir()))
