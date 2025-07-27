
using MagicSimulate
cd(@__DIR__)

# simulate data
designcode = "2star-self1"
fhaplofile = "fhaplo.vcf"
ncluster = 1
nsnpchr = 100
simfhaplo(nsnp=ncluster*nsnpchr, nparent=8,
    isfounderinbred = false,
    chrlen=100*ones(ncluster),
    outfile=fhaplofile
)
using Distributions
magicsimulate(fhaplofile,designcode;
    popsize = 100,
    isfounderinbred = false, 
    seqfrac = 1.0,
    allelicbias = Beta(5,5),
    allelicoverdispersion = Exponential(0.05),
    allelicdropout = Beta(2,18),
    seqdepth = Gamma(2,10),    
    outstem = designcode*"_outbred",
)

rm(fhaplofile)
rm.(filter(x->occursin(".log",x), readdir()))
rm.(filter(x->occursin("_ped",x), readdir()))
rm.(filter(x->occursin("_truecontfgl",x), readdir()))
