cd(@__DIR__)
pwd()

# simulate data
using MagicBase
using MagicSimulate
using Distributions
dataid="sim"

fhaplofile = dataid*"_fhaplo.vcf"
ncluster = 1
nsnpchr = 200
@time simfhaplo(
    nsnp=ncluster*nsnpchr, nparent=32,
    chrlen=100*ones(ncluster),
    outfile=fhaplofile
)

# designcode = "2star-self1"
designcode = "2ril-self1"
designinfo = MagicBase.parsedesign(designcode)
magicped = formmagicped(designinfo,50)

# magicped = MagicBase.generate_magicped(;
#     designcodes = ["bc2-self1","P2/P3=>DH","P4/3/P2/P3//P5/P6=>2","2ril-self1"],
#     founders = ["P1||P2","NA", "NA", "P5||P6"],
#     subpopsizes = 100*ones(Int, 4));
# plotmagicped(magicped)

pedfile = dataid*"_ped.csv"
savemagicped(pedfile,magicped)

epsf = 0.01
epso = 0.01
fmiss = 0.999
@time magicsimulate(fhaplofile,pedfile;    
    foundererror = Beta(1, 2/epsf-1.0),
    offspringerror = Beta(1, 2/epso-1.0),
    foundermiss = Beta(2,2/fmiss-2),
    offspringmiss = Beta(1,99),    
    seqfrac = 0.0,    
    allelebalancemean = Beta(3,3),
    allelebalancedisperse = Exponential(0.2),
    alleledropout = Beta(1,199),
    seqerror = Beta(1,199),
    seqdepth = Gamma(2,20),
    outstem= dataid,    
)

rm(fhaplofile)
