
using MagicSimulate, MagicBase
cd(@__DIR__)

# simulate data
using MagicSimulate
using Distributions
dataid="sim"
fhaplofile = "sim_fhaplo.vcf.gz"
ncluster = 5
nsnpchr = 250

@time simfhaplo(
    isfounderinbred =false,
    nsnp=ncluster*nsnpchr,
    nparent=20,
    chrlen=100*ones(ncluster),
    outfile=fhaplofile
)

designcode = "6star-self1"
designinfo = MagicBase.parsedesign(designcode)
magicped = formmagicped(designinfo,400)
pedfile = dataid*"_ped.csv"
savemagicped(pedfile,magicped)

epsf = epso = 0.01
@time magicsimulate(fhaplofile,pedfile;
    isfounderinbred = false,
    # seqfrac = 0.5,
    foundererror = Uniform(epsf,epsf+0.001),
    offspringerror = Uniform(epso,epso+0.001),
    foundermiss = Beta(1,4),
    offspringmiss = Beta(1,4),
    outstem= dataid,
)

#  clear up
rm(fhaplofile)
rm.(filter(x->occursin("fgl.",x), readdir()))
