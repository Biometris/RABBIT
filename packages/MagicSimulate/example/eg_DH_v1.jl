
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
pedinfo = "sim_DH_magicped.csv"
magicped
savemagicped(pedinfo,magicped)
plotped(magicped.designinfo)
plotmagicped(magicped)

# pedinfo = "4ril-self4"
# sim founder genotypes
using MagicSimulate
using Distributions
dataid="sim"
fhaplofile = dataid*"_fhaplo.csv"
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


# rm.(filter(x->occursin("outstem",x), readdir()))

popsize= 100
chrlen = [150,120,80,100,60]
isobligate = false
interference = 0
isfounderinbred =  true
foundererror = offspringerror = Beta(1,95)
foundermiss = offspringmiss = Beta(1,4)
seqfrac = 0.5
seqdepth = Gamma(2,10)
baseerror = 0.001
select_proportion = 0.1
select_nqtl = 50
fhaplofile = "sim_fhaplo.csv"
pedinfo = "sim_DH_magicped.csv"

founderhaplo = MagicBase.formfhaplo(fhaplofile; workdir,isfounderinbred);
if pedinfo[end-3:end]==".csv"
    magicped = readmagicped(pedinfo; workdir)
else
    magicped = formmagicped(pedinfo,popsize)
    MagicBase.setfounderid!(magicped,founderhaplo.magicped.founderinfo[!,:individual])
end


nondiv = allunique(magicped.offspringinfo[!,:member])
isplotped = (!nondiv) || (length(magicped.designinfo.member)<1000)
if isplotped
    gped = plotmagicped(magicped)
end

contfgl, magicfgl = simfgl(founderhaplo,magicped;
    isfounderinbred,isobligate,interference,
    select_nqtl,select_proportion)
0
typeof(res)
contfgl.offspringgeno[1]
contfgl.offspringgeno[1][:,2]

ishomo = magicped.offspringinfo[!,:ishomozygous]

chr = 1
offgeno = contfgl.offspringgeno[chr]




truegeno = genedrop(magicfgl; isfounderinbred,isphased=true)
obsgeno = simgeno(truegeno;foundererror,offspringerror,
    foundermiss,offspringmiss,seqfrac,seqdepth,baseerror)
# 0
obsgeno.markermap[1]
obsgeno.foundergeno[1]
obsgeno.offspringgeno[1]
savegenodata("test.vcf",obsgeno)

