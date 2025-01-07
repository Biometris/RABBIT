
using Revise
using MagicBase
using MagicLinkage
cd(@__DIR__)
pwd()



isfounderinbred = false
genofile = "sim_magiccall_geno.vcf.gz"
pedfile = "sim_magicfilter_ped.csv"
ldfile = "sim_magicld.csv.gz"


# magicgeno = formmagicgeno(genofile,pedfile; isfounderinbred);
# merge_chromosome!(magicgeno);
# MagicLinkage.magiclinkage!(magicgeno;isfounderinbred,ldfile); 


magiclinkage(genofile,pedfile; isfounderinbred, ldfile, outstem="sim")
