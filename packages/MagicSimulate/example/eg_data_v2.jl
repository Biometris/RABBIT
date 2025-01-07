
using MagicSimulate, MagicBase
cd(@__DIR__)


# simulate data
using MagicSimulate
using Distributions
dataid="sim"
fhaplofile = "sim_fhaplo.vcf.gz"
ncluster = 5
nsnpchr = 500

@time simfhaplo(
    nsnp=ncluster*nsnpchr,
    nparent=6,
    chrlen=100*ones(ncluster),
    outfile=fhaplofile
)



designcode = "6star-self1"
designinfo = MagicBase.parsedesign(designcode)
magicped = formmagicped(designinfo,500)
pedfile = dataid*"_ped.csv"
savemagicped(pedfile,magicped)


epsf = epso = 0.01
@time magicsimulate(fhaplofile,pedfile;
    foundererror = Uniform(epsf,epsf+0.001),
    offspringerror = Uniform(epso,epso+0.001),
    # foundermiss = Beta(1,4),
    # offspringmiss = Beta(1,4),
    seqfrac = 0.5,
    outstem= dataid,
    error_randallele = 1.0,
    # verbose = true,
)

#  clear up
cd(@__DIR__)
rm(fhaplofile)
rm.(filter(x->occursin("sim",x), readdir()))



