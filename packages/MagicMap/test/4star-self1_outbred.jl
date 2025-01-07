
# using Revise
using MagicBase
using MagicMap
using Test
cd(@__DIR__)

# construct map
dataid = "4star-self1_outbred"
genofile = dataid*"_magicfilter_geno.vcf.gz"
pedfile = dataid*"_magicfilter_ped.csv"
outstem=dataid
ncluster = 2
outfiles = magicmap(genofile,pedfile; 
    maxncluster = 10,
    isfounderinbred = false,    
    ncluster = 2,     
    outstem,    
    verbose=false,
)

# test results
truemapfile = string(dataid,"_magicsimulate_truegeno.csv.gz")
estmapfile = outstem*"_magicmap_construct_map.csv.gz"
acc = mapaccuracy(truemapfile,estmapfile; isgroupacc=true)

# println(acc)
@test acc.groupacc > 0.9
@test sum(abs.(acc.orderacc))/length(acc.orderacc) > 0.6

#  clean up
rm.(filter(x->occursin(outstem*"_magicmap",x), readdir()))
