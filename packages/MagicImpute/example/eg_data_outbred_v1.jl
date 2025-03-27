using Revise
cd(@__DIR__)
pwd()

# simulate data
using MagicBase
using MagicSimulate
using Distributions
dataid="sim"

isfounderinbred = false

fhaplofile = dataid*"_fhaplo.vcf.gz"
ncluster = 1
nsnpchr = 200
nparent = 4
@time simfhaplo(;
    nsnp=ncluster*nsnpchr, 
    nparent,
    isfounderinbred,
    chrlen=100*ones(ncluster),
    outfile=fhaplofile
)

epsf = 0.01
epso = 0.01
missf = 0.1
misso = 0.1
# pedcode = string(nparent, "star-self0")
pedcode = string(nparent, "ril-self0")
@time magicsimulate(fhaplofile,pedcode;
    popsize=200,
    isfounderinbred,    
    foundererror = Beta(2, 2/epsf-2.0),
    offspringerror = Beta(2, 2/epso-2.0),
    foundermiss = Beta(2,2/missf-2),
    offspringmiss = Beta(2,2/misso-2),
    error_randallele = 1.0, 
    seqfrac = 0.0,
    seqerror = Beta(2,2/0.005-2),
    allelebalancemean = Beta(3,3),
    allelebalancedisperse = Exponential(0.2),
    # alleledropout = Beta(1,199),
    seqdepth = Gamma(2, 10),
    outstem= dataid,    
)

rm(fhaplofile)

