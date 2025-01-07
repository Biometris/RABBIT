using Revise
using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
using Pedigrees
using Plots


cd(@__DIR__)
pwd()


pedfile="pedigree.csv"
ped=Pedigrees.readped(pedfile)
Pedigrees.plotped(ped)
founders = [string("P",i) for i=1:ped.nfounder]
ped2 = Pedigrees.setfounderid(ped,founders)
g = Pedigrees.plotped(ped,fontsize=10)
g2 = Pedigrees.plotped(ped2,fontsize=10)
plot(g,g2)

Pedigrees.plotped(ped)