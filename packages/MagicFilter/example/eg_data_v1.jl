
using MagicBase, MagicSimulate
cd(@__DIR__)
pwd()

# reset pedfile so that founders are consistent with simfhaplo.
pedfile = "DH_magicped.csv"
magicped = readmagicped(pedfile)
ped = magicped.designinfo
newfounderid = [string("P",i) for i in 1:ped.nfounder]
magicped.designinfo = MagicBase.setfounderid(ped,newfounderid)
pedinfo = "sim_DH_magicped.csv"
savemagicped(pedinfo,magicped)
MagicBase.plotped(magicped.designinfo)
plotmagicped(magicped)

# sim founder genotypes
using Distributions
dataid="sim"
fhaplofile = dataid*"_fhaplo.csv"
nparent = 9
ncluster = 5
nsnpchr = 200
@time simfhaplo(;
    nsnp=ncluster*nsnpchr, nparent,
    chrlen=100*ones(ncluster),
    outfile=fhaplofile
)

# simulate offspring genotypes
@time magicsimulate(fhaplofile,pedinfo;
    foundermiss = Beta(1,4),
    offspringmiss = Beta(1,4),
    seqfrac = 0.5,
    outstem= dataid,
)
