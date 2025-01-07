using MagicSimulate, MagicBase
using CSV, DataFrames
cd(@__DIR__)


# simulate founder haplotypes
fhaplofile = "fhaplo.vcf.gz"
ncluster = 1
nsnpchr = 200
simfhaplo(nsnp=ncluster*nsnpchr, nparent=6,
    chrlen=100*ones(ncluster),
    outfile=fhaplofile
)

magicped = MagicBase.generate_magicped(;
    designcodes = ["bc3-self2","P2/P3=>DH","P4/3/P2/P3//P5/P6=>3","2ril-self6"],
    founders = ["P1||P2","NA", "NA", "P5||P6"],
    subpopsizes =10*ones(Int, 4));
pedfile = "multipop_ped.csv"
savemagicped(pedfile, magicped)
# plotmagicped(magicped)

using Distributions
magicsimulate(fhaplofile,pedfile;    
    foundermiss = Beta(2,6),
    offspringmiss = Beta(2,6),    
    outstem = "multipop",
)

#  clean up
rm(fhaplofile)
rm.(filter(x->occursin(".log",x), readdir()))
rm.(filter(x->occursin("_ped.png",x), readdir()))
rm.(filter(x->occursin("_truecontfgl",x), readdir()))
rm.(filter(x->occursin("_truefgl",x), readdir()))
