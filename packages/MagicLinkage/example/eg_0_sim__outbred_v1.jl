
using MagicSimulate, MagicBase
cd(@__DIR__)

# using MagicFilter
# initprob = 0.25*ones(4)
# initprob = 0.5*ones(2)
# freq_g22 = MagicFilter.possible_freq_g22(initprob,0.01)
# MagicFilter.cal_min_popsize(freq_g22,0.99)


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
magicped = formmagicped(designinfo,200)
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
