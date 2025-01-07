using Revise
using MagicIO
cd(@__DIR__)
pwd()
using Pedigrees

design = MagicIO.parsedesign("4ril-self3")
plotped(design.pedigree)

design = MagicIO.parsedesign("11star-self3")
# magicped = formmagicped(ped,100)
plotped(design.pedigree)

magicped = readmagicped("multipop_designcode.csv")
plotmagicped(magicped)
Pedigrees.ped2df(magicped.designinfo)

magicped = readmagicped("multipop_magicped.csv")
plotmagicped(magicped)
# savemagicped("multipop_magicped.csv",magicped)
using Plots
savefig("multipop.png")

pedfile = "multipop_designcode.csv"
peddict = MagicIO.readdlm2dict(pedfile)
designdf = peddict["designinfo"]

