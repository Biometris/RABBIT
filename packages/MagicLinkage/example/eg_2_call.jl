using MagicCall
cd(@__DIR__)
pwd()


isfounderinbred = true 

dataid = "sim"
genofile=string(dataid,"_magicfilter_geno.vcf.gz")
pedfile = string(dataid,"_magicfilter_ped.csv")
outstem = dataid
@time magiccall(genofile,pedfile;    
    isfounderinbred,
    outstem 
)