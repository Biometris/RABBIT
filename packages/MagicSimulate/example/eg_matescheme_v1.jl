
using Revise
using MagicBase, MagicSimulate
using Distributions
cd(@__DIR__)
pwd()

using Pedigrees
using DataFrames
using StatsBase
using Random
using GraphRecipes, Plots

matescheme = MateScheme(2,["Pairing","DH"],[1,1])
matescheme = MateScheme(4,["Pairing","Selfing"],[2,6])
matescheme = MateScheme(8,["Pairing","Selfing"],[3,6])
# matescheme = MateScheme(16,["Pairing","Selfing"],[4,3])
# matescheme = MateScheme(8,["Pairing","Sibling"],[2,20])
# matescheme = MateScheme(19,["RM1_NE_20", "RM1_E", "Selfing"],[1,4,3])

magicped = rand_magicped(matescheme)
ped = magicped.designinfo
g = plotped(ped)
# savefig(g,"ped.png")
# g2 = plotmagicped(magicped)
# savefig(g2,"magicped.png")
