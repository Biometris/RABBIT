
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

# TODO: sim gbs data
# TODO: sim mixed genotypign systems among markers/individuals
# TODO: sim random mating schemes
# TODO: sim related founders
# TODO: sim progeny selection
# TODO: sim genomic re-arrangement

# designcode = "8ril-self4"
# designcode = "5star-self3"
# magicped = formmagicped(designcode,popsize)
# plotmagicped(magicped)
# outfile = "test.csv"
# savemagicped(outfile,magicped)
# magicped2 = readmagicped(outfile)
# plotmagicped(magicped2)

fhaplofile = "ArabMAGIC_founderhaplo.csv"
magicpedfile ="ArabMAGIC_magicped.csv"
@time magicsimulate(fhaplofile,magicpedfile)

fhaplofile = "founderhaplo_inbred.csv"
# designcode = "4ril-self4"
designcode = "6star-self5"
@time magicsimulate(fhaplofile,designcode)

fhaplofile = "founderhaplo_outbred.csv"
designcode = "8ril-self4"
# designcode = "5star-self3"
@time magicsimulate(fhaplofile,designcode)
