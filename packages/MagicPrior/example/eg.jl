using Revise
# using Pkg
# Pkg.activate(joinpath(@__DIR__,".."))
using MagicPrior
cd(@__DIR__)
pwd()
using Pedigrees

names(MagicPrior)

filename="CC_ped.csv"
ped = Pedigrees.readped(filename)
nfounder = ped.nfounder
founderfgl= [i for i = 1:nfounder, j = 1:2]
# founderfgl= [2*(i-1)+j for i = 1:nfounder, j = 1:2]
indarray = ped.member[end-1:end]
indarray2 = [i for i = indarray, j = 1:2]
@time for isautosome=[true,false],isfglexch=[true,false]
    res = magicorigin(ped,founderfgl=founderfgl,memberlist=indarray,
        isautosome=isautosome,isfglexch=isfglexch,isconcise=false)
end

@time for isautosome=[true,false],isfglexch=[true,false]
    prior =magicprior(ped,founderfgl=founderfgl,memberlist=indarray,
        isautosome=isautosome,isfglexch=isfglexch,isconcise=false)
end

# MagicPrior.isexchcompatible(ped,isautosome=true)
# MagicPrior.isexchcompatible(ped,isautosome=false)

isautosome=true
isfglexch=false
prior =magicprior(ped,founderfgl=founderfgl,memberlist=indarray,
    isautosome=isautosome,isfglexch=isfglexch,isconcise=false)

prior["53"]
println(prior["53"])

a=prior["53"].contmarkov[3]

using SparseArrays
i,v = findnz(a[1])

i2,j2,v2 = findnz(a[2])


issubset(i2,i)
unique(i2)

keys(prior["53"])
