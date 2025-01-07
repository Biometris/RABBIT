
using MagicSimulate, MagicBase
cd(@__DIR__)
pwd()

# simulate data
designcode = "4ril-self3"
fhaplofile = "fhaplo.vcf.gz"
ncluster = 2
nsnpchr = 200
simfhaplo(nsnp=ncluster*nsnpchr, nparent=4,
    chrlen=100*ones(ncluster),
    outfile=fhaplofile
)
magicsimulate(fhaplofile,designcode;
    popsize = 50,
    seqfrac = 0.0,
    outstem = designcode*"_array",
)
magicsimulate(fhaplofile,designcode;
    popsize = 50,
    seqfrac = 0.5,
    outstem = designcode*"_mix",
)

#  clean up
rm(fhaplofile)
rm.(filter(x->occursin(".log",x), readdir()))
rm.(filter(x->occursin("_ped",x) && occursin(designcode,x), readdir()))
rm.(filter(x->occursin("_truecontfgl",x), readdir()))
rm.(filter(x->occursin("_truefgl",x), readdir()))





