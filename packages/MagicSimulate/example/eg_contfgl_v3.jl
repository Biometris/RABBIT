using Revise
using MagicSimulate
using MagicBase
using Distributions
cd(@__DIR__)
pwd()

using Plots

# matescheme = MateScheme(2,["Pairing","DH"],[1,1])
matescheme = MateScheme(2,["Pairing","Selfing"],[1,1])
matescheme = MateScheme(4,["Pairing","DH"],[2,1])
matescheme = MateScheme(4,["Pairing","Selfing"],[2,4])
# matescheme = MateScheme(8,["Pairing","Selfing"],[3,1])
# matescheme = MateScheme(8,["Pairing"],[3])
# matescheme = MateScheme(8,["Pairing","Sibling"],[2,10])
# matescheme = MateScheme(19,["FullDiallel", "RM1_E", "Selfing"],[1,4,6])
# designcode = "4ril-self1"

pedinfo = matescheme
# pedinfo = designcode
@time outfiles = magicsimulate(pedinfo;
    popsize= 200,        
    # isfounderinbred = false,
    chrlen = 100*ones(5),
    # outstem=nothing
)

# rm.(filter(x->occursin("sim",x), readdir()))

