using Revise
cd(@__DIR__)
pwd()

# simulate data
using MagicBase
using MagicSimulate
using Distributions
dataid="sim"

isfounderinbred = true

fhaplofile = dataid*"_fhaplo.vcf.gz"
ncluster = 2
nsnpchr = 500
@time simfhaplo(;
    nsnp=ncluster*nsnpchr, 
    nparent=4,
    isfounderinbred,
    chrlen=100*ones(ncluster),
    outfile=fhaplofile
)

epsf = epso = 0.01
# pedcode = "6star-self1"
# pedcode = "2ril-self1"
pedcode = "2ril-self6"
@time magicsimulate(fhaplofile,pedcode;
    popsize=200,
    isfounderinbred,
    # foundererror = Uniform(epsf,epsf+0.001),
    # offspringerror = Uniform(epso,epso+0.001),
    foundererror = Beta(2, 2/epsf-2.0),
    offspringerror = Beta(2, 2/epso-2.0),
    foundermiss = Beta(2,18),
    offspringmiss = Beta(2,18),
    seqfrac = 0.0,
    allelicbias = Beta(5,5),
    allelicoverdispersion = Exponential(0.05),
    allelicdropout = Beta(2,18),
    seqdepth = Gamma(2,20),
    outstem= dataid,    
)


rm(fhaplofile)


