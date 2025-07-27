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
ncluster = 5
nsnpchr = 200
nparent = 2
@time simfhaplo(;
    nsnp=ncluster*nsnpchr, 
    nparent,
    isfounderinbred,
    chrlen=100*ones(ncluster),
    outfile=fhaplofile
)

epsf = 0.05
epso = 0.05
missf = 0.05
misso = 0.05
seqerr = 0.002
# pedcode = string(nparent, "star-self0")
pedcode = string(nparent, "ril-self0")
@time magicsimulate(fhaplofile,pedcode;
    popsize=200,
    isfounderinbred,    
    foundererror = Beta(1, 1/epsf-1),
    offspringerror = Beta(1, 1/epso-1),
    foundermiss = Beta(1,1/missf-1),
    offspringmiss = Beta(1,1/misso-1),
    error_randallele = 1,
    seqfrac = 1.0,
    baseerror = Beta(1,1/seqerr-1),
    allelicbias = Beta(3,3),
    allelicoverdispersion = Exponential(0.3),    
    seqdepth = Gamma(1, 40),
    outstem= dataid,    
)

rm(fhaplofile)

