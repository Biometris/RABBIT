cd(@__DIR__)
pwd()

# simulate data
using MagicBase
using MagicSimulate
using Distributions
dataid="sim"

isfounderinbred = true


fhaplofile = dataid*"_fhaplo.vcf"
ncluster = 1
nsnpchr = 200
@time simfhaplo(; 
    isfounderinbred, 
    nsnp=ncluster*nsnpchr, nparent=7, 
    chrlen=100*ones(ncluster),
    outfile=fhaplofile
)

magicped = MagicBase.generate_magicped(;
    designcodes = ["P1/P2=>DH","P1/P3=>2","P2/P3=>1","P4/P5//P6/P7=>1"],    
    subpopsizes = 2*ones(Int, 4));
plotmagicped(magicped)

pedfile = dataid*"_ped.csv"
savemagicped(pedfile,magicped)

epsf = 0.02
epso = 0.02
missf = 0.5
@time magicsimulate(fhaplofile,pedfile;    
    isfounderinbred, 
    foundererror = Beta(1, 1/epsf-1.0),
    offspringerror = Beta(1, 1/epso-1.0),
    foundermiss = Beta(1,1/missf-1),
    offspringmiss = Beta(1,10),
    error_randallele = 1.0,
    seqfrac = 0.0,    
    allelebalancemean = Beta(5,5),
    allelebalancedisperse = Exponential(0.05),
    alleledropout = Beta(1,199),
    seqerror = Beta(1,100),
    seqdepth = Gamma(2,10),
    outstem= dataid,    
)

rm(fhaplofile)
