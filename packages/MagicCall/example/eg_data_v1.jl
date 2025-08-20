
using MagicSimulate, MagicBase
cd(@__DIR__)

# simulate data
using MagicSimulate
using Distributions
dataid="sim"
fhaplofile = "sim_fhaplo.vcf.gz"
ncluster = 2
nsnpchr = 200
nparent = 2
isfounderinbred = false

@time simfhaplo(;
    isfounderinbred,
    nsnp=ncluster*nsnpchr,
    nparent,
    chrlen=100*ones(ncluster),
    # chrlen = [100,200,50,150,120], 
    outfile=fhaplofile
)

designcode = isfounderinbred ? string(nparent, "star-self1") : string(nparent, "star-self0")
# designcode = string(naprent, "ril-self6")
designinfo = MagicBase.parsedesign(designcode)
magicped = formmagicped(designinfo,200)
pedfile = dataid*"_ped.csv"
savemagicped(pedfile,magicped)

epsf = 0.02
epso = 0.02
baseerror = 0.002
missfreq = 0.2
depth = 40
@time magicsimulate(fhaplofile,pedfile;    
    isfounderinbred,
    foundererror = Beta(1, 1/epsf-1.0),
    offspringerror = Beta(1, 1/epso-1.0),
    foundermiss = Beta(1,1/missfreq-1),
    offspringmiss = Beta(1,1/missfreq-1),    
    error_randallele = 1.0, 
    seqfrac = 0.0,    
    allelicbias = Beta(3,3),
    allelicoverdispersion = Exponential(0.3),
    # allelicdropout = Beta(1,1/0.02-1),
    baseerror = Beta(1,1/baseerror-1),
    seqdepth = Gamma(1,depth),
    outstem= dataid,    
    nplot_subpop = 1, 
)

#  clear up
rm(fhaplofile)
rm.(filter(x->occursin("fgl.",x), readdir()))
