using Revise
cd(@__DIR__)
pwd()

# simulate data
using MagicBase
using MagicSimulate
using Distributions
dataid="sim"

fhaplofile = dataid*"_fhaplo.vcf.gz"
ncluster = 2
nsnpchr = 250
@time simfhaplo(
    nsnp=ncluster*nsnpchr, 
    nparent=20,
    isfounderinbred = true,
    chrlen=100*ones(ncluster),
    outfile=fhaplofile
)

epsf = epso = 0.01
pedfile = "multipop_designcode.csv"
@time magicsimulate(fhaplofile,pedfile;
    # popsize=100,
    isfounderinbred = true,
    # foundererror = Uniform(epsf,epsf+0.001),
    # offspringerror = Uniform(epso,epso+0.001),
    foundererror = Beta(2, 2/epsf-2.0),
    offspringerror = Beta(2, 2/epso-2.0),
    foundermiss = Beta(2,8),
    offspringmiss = Beta(2,8),
    seqfrac = 1.0,
    allelebalance = Beta(2,4),
    # allelebalance = Uniform(10^(-6),10^(-5)),
    seqdepth = Gamma(2,10),
    outstem= dataid,
    isplot=true
)

# rm(fhaplofile)
