
# registry rm General
# registry add https://github.com/JuliaRegistries/General
# registry up General
using Revise
using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
using MagicSimulate
using Pedigrees
using DataFrames
using Distributions
using DelimitedFiles
using CSV
cd(@__DIR__)
pwd()
comment_char = '#'
delim = ','
workdir = pwd()
outstem = "outstem"
popsize= 100
chrlen = [150,120,80,100,60]
isobligate = false
interference = 0


# TODO: write savemagicgeno for contfgl/magicfgl
# TODO: test read/save magicfgl/magicgeno
# TODO: sim snp-array data
# TODO: sim gbs data
# TODO: for outbred founders

# TODO: sim random mating schemes
# TODO: sim progeny selection
# TODO: sim genomic re-arrangement

designcode = "8ril-self4"
# designcode = "5star-self3"
@time magicped = formmagicped(designcode,popsize)
plotmagicped(magicped)
outfile = "test.csv"
savemagicped(outfile,magicped)
readmagicped(outfile)
@time contfgl = simcontfgl(magicped;
    chrlen,isobligate,interference,outstem
)

#
# genofile = "Soybean_NAM_geno.csv"
# pedfile = "Soybean_NAM_ped.csv"
# @time magicgeno = formmagicgeno(genofile,pedfile)

fhaplofile = "ArabMAGIC_founderhaplo.csv"
magicpedfile ="ArabMAGIC_magicped.csv"
@time magicfgl = simfgl(fhaplofile,magicpedfile;
    isobligate,interference,outstem
)
