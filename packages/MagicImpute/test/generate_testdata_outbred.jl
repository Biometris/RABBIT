
using MagicSimulate, MagicBase
cd(@__DIR__)

# simulate data
designcode = "4ril-self3"
fhaplofile = "fhaplo.vcf"
ncluster = 1
nsnpchr = 200
simfhaplo(nsnp=ncluster*nsnpchr, nparent=8,
    isfounderinbred = false,
    chrlen=100*ones(ncluster),
    outfile=fhaplofile
)
magicsimulate(fhaplofile,designcode;
    popsize = 50,
    isfounderinbred = false,
    outstem = designcode*"_outbred",
)

rm(fhaplofile)
rm.(filter(x->occursin(".log",x), readdir()))
rm.(filter(x->occursin("_ped",x) && occursin(designcode,x), readdir()))
rm.(filter(x->occursin("_truecontfgl",x), readdir()))
