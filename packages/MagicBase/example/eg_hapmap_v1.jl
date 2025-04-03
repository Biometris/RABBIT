using Revise
using MagicBase
cd(@__DIR__)
pwd()


hapmapfile = "data_hmp.txt"
isdelmultiallelic = true
delim = r" +"
outstem  = "outstem"
workdir = pwd()
verbose=true

isfile(MagicBase.hapmap2vcf(hapmapfile))
