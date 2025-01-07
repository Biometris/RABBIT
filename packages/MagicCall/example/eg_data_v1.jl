
using MagicSimulate, MagicBase
cd(@__DIR__)

# simulate data
using MagicSimulate
using Distributions
dataid="sim"
fhaplofile = "sim_fhaplo.vcf.gz"
ncluster = 2
nsnpchr = 300

isfounderinbred = true

@time simfhaplo(;
    isfounderinbred,
    nsnp=ncluster*nsnpchr,
    nparent=20,
    chrlen=100*ones(ncluster),
    # chrlen = [100,200,50,150,120], 
    outfile=fhaplofile
)

designcode = "6star-self4"
# designcode = "4ril-self4"
designinfo = MagicBase.parsedesign(designcode)
magicped = formmagicped(designinfo,200)
pedfile = dataid*"_ped.csv"
savemagicped(pedfile,magicped)

epsf = 0.05
epso = 0.05
@time magicsimulate(fhaplofile,pedfile;
    isfounderinbred,    
    foundererror = Beta(2, 2/epsf-2.0),
    offspringerror = Beta(2, 2/epso-2.0),    
    foundermiss = Beta(1,4),
    offspringmiss = Beta(1,4),    
    seqfrac = 1.0,
    allelebalancemean = Beta(5,5),
    allelebalancedisperse = Exponential(0.05),
    alleledropout = Beta(2,18),
    seqdepth = Gamma(2,10),
    outstem= dataid,
)

#  clear up
rm(fhaplofile)
rm.(filter(x->occursin("fgl.",x), readdir()))
