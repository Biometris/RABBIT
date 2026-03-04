
# using Distributed
# addprocs(5) # add worker processes on local machine
# @info string("nprocs=", nprocs())
# @everywhere  using MagicFilter

using Revise
using MagicBase, MagicFilter
cd(@__DIR__)
pwd()
using CSV, DataFrames

genofile = "sim_magicsimulate_geno.vcf.gz"
pedfile = "sim_magicsimulate_ped.csv"
isfounderinbred = false
magicgeno = formmagicgeno(genofile,pedfile; isfounderinbred);
@time MagicFilter.filter_marker!(magicgeno;
    isfounderinbred, 
    # maxhetero = 0.9,
    isparallel = false,
);
0

# cd(@__DIR__)
# outid = "sim"
# rm.(filter(x->occursin(outid,x), readdir()))

