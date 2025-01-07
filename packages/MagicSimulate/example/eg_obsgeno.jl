using Revise
using MagicSimulate
using MagicBase
using Distributions
cd(@__DIR__)
pwd()

dataid = "sim"
fhaplofile = dataid*"_fhaplo.vcf.gz"
ncluster = 2
nsnpchr = 250
@time simfhaplo(
    nsnp=ncluster*nsnpchr, 
    nparent=20,
    isfounderinbred = true,
    chrlen=100*ones(ncluster),
    outfile=fhaplofile
)



fhaplofile = "sim_fhaplo.vcf"
pedinfo = "4star-self3"
chrlen = 100*ones(5)
isfounderinbred,isobligate,interference = false, false,0
missingstring = "NA"
workdir = pwd()
popsize= 100
foundermiss=Beta(1,9)
offspringmiss=Beta(1,9)
error_randallele=0.5
foundererror=Beta(1,199)
offspringerror=Beta(1,199)   
seqfrac=0.5
seqerror = Beta(1,199)
allelebalance = Beta(1,99)    
seqdepth = Gamma(2, 1)
seqdepth_disperse = 1.0

founderhaplo = MagicBase.formfhaplo(fhaplofile; 
    formatpriority = ["GT"], isphysmap = nothing, 
    missingstring = unique(vcat(missingstring, "missing")),
    workdir,isfounderinbred)
magicped = MagicSimulate.pedinfo2magicped(pedinfo; popsize,workdir)       
MagicBase.setfounderid!(magicped,founderhaplo.magicped.founderinfo[!,:individual])
# contfgl, magicfgl = MagicSimulate.simfgl(founderhaplo,magicped;
#     isfounderinbred,isobligate,interference,select_prop=1.0)
# truegeno = MagicSimulate.genedrop(magicfgl; isfounderinbred,isphased=true)
# obsgeno = MagicSimulate.simgeno(truegeno,magicfgl)

##############################################################

select_nqtl = 10
select_prop = 0.5
isfounderinbred = false
MagicSimulate.checkfounder!(founderhaplo,magicped)
markermap = deepcopy(founderhaplo.markermap)
# chrlen+1: ensure chrlen is larger than the last marker position
chrid = [i[1,:linkagegroup] for i=markermap]
chrlen = [ceil(i[end,:poscm]-i[1,:poscm])+1 for i=markermap]
markerpos = [i[!,:poscm] .- i[1,:poscm] for i=markermap]
# contfgl = MagicSimulate.simcontfgl(magicped; chrid,chrlen,isfounderinbred,isobligate,interference,
#     select_nqtl,select_prop)  

##############################################################

