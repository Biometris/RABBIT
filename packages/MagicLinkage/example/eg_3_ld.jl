
using MagicLD
cd(@__DIR__)
pwd()

genofile = "sim_magiccall_geno.vcf.gz"
skipind = 6

magicld(genofile,skipind; outstem="sim")
