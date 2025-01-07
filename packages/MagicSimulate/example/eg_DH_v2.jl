
using Revise
# using Pkg
# Pkg.activate(joinpath(@__DIR__,".."))
using MagicBase, MagicSimulate
cd(@__DIR__)
pwd()
workdir = pwd()

# reset pedfile so that founders are consistent with simfhaplo. 
pedfile = "DH_magicped.csv"
# pedfile = "DH_magicped5_p31_pp79_o4256.csv"
# pedfile = "DH_magicped7_p26_pp29_o5194.csv"
magicped = readmagicped(pedfile)
ped = magicped.designinfo
newfounderid = [string("P",i) for i in 1:ped.nfounder]
magicped.designinfo = MagicBase.setfounderid(ped,newfounderid)
pedinfo = "outstem_DH_magicped.csv"
magicped
savemagicped(pedinfo,magicped)
MagicBase.plotped(magicped.designinfo)
plotmagicped(magicped)

# pedinfo = "4ril-self4"
# sim founder genotypes
using MagicSimulate
using Distributions
fhaplofile = "outstem_fhaplo.csv"
nparent = ped.nfounder
# nparent = 4
ncluster = 5
nsnpchr = 200
@time simfhaplo(;
    nsnp=ncluster*nsnpchr, nparent,
    chrlen=100*ones(ncluster),
    outfile=fhaplofile
)


# simulate offspring genotypes
@time magicsimulate(fhaplofile,pedinfo;  
    popsize = 200,      
    foundermiss = Beta(1,4),
    offspringmiss = Beta(1,4),    
)

rm.(filter(x->occursin("outstem",x), readdir()))
