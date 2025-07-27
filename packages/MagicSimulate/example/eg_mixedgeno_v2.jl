
using Revise
# using Pkg
# Pkg.activate(joinpath(@__DIR__,".."))
using MagicBase, MagicSimulate
cd(@__DIR__)
pwd()
workdir = pwd()

MagicBase.formfhaplo("sim_fgeno_haplo.csv.gz")

@time magicsimulate("sim_fgeno_haplo.csv.gz","pedinfo.csv";
    seqfrac = 0.5,
    # outext = ".vcf"
    # outstem = nothing
)

@time magicsimulate("sim_fgeno_phased.csv.gz","pedinfo.csv";
    seqfrac = 0.5,
    isfounderinbred = false
)

rm.(filter(x->occursin("outstem",x), readdir()))

using Distributions
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
fhaplofile = isfounderinbred ? "sim_fgeno_haplo.csv.gz" : "sim_fgeno_phased.csv.gz"
pedinfo = "pedinfo.csv"


founderhaplo = MagicBase.formfhaplo(fhaplofile; workdir,isfounderinbred);
founderhaplo

savegenodata("test.csv",founderhaplo)

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
founderhaplo
magicped
contfgl, magicfgl = MagicSimulate.simfgl(founderhaplo,magicped;
    isfounderinbred,isobligate,interference,
    select_nqtl,select_proportion)
truegeno = MagicSimulate.genedrop(magicfgl; isfounderinbred,isphased=true)
obsgeno = MagicSimulate.simgeno(truegeno;foundererror,offspringerror,
    foundermiss,offspringmiss,seqfrac,seqdepth,baseerror)

0

MagicSimulate.checkfounder!(founderhaplo,magicped)
markermap = deepcopy(founderhaplo.markermap)
# chrlen+1: ensure chrlen is larger than the last marker position
chrid = [i[1,:linkagegroup] for i=markermap]
chrlen = [ceil(i[end,:poscm]-i[1,:poscm])+1 for i=markermap]
markerpos = [i[!,:poscm] .- i[1,:poscm] for i=markermap]
contfgl = MagicSimulate.simcontfgl(magicped; chrid,chrlen,isfounderinbred,isobligate,interference,
    select_nqtl,select_proportion)
0



obsgeno.markermap[1]
obsgeno.foundergeno[1]
obsgeno.offspringgeno[1]
savegenodata("test.vcf",obsgeno)
