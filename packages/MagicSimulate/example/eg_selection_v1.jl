
using Revise
using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
using MagicBase, MagicSimulate
using Pedigrees
using DataFrames
using Distributions
using CSV
cd(@__DIR__)
pwd()
workdir = pwd()

popsize= 100
chrlen = [150,120,80,100,60]
isobligate = false
interference = 0
isfounderinbred =  true
select_proportion = 0.1
select_nqtl = 50

fhaplofile = "ArabMAGIC_founderhaplo.csv"
pedinfo = "pedinfo.csv"
@time magicsimulate(fhaplofile,pedinfo; select_proportion)

rm.(filter(x->occursin("outstem",x), readdir()))



pedinfo = "pedinfo.csv"
@time magicped = readmagicped(pedinfo)
plotmagicped(magicped)
@time contfgl = simcontfgl(magicped;chrlen,isobligate,interference)



designped = magicped.designinfo
indexped = Pedigrees.toindexped(designped)
nf = indexped.nfounder
seltrait = select_proportion ==1.0 ? nothing : TraitQtl(chrlen, select_nqtl, nf)    
nondiv = allunique(magicped.offspringinfo[!,:member])

pedfgl = setfounderfgl(indexped,chrlen,isfounderinbred)        
noff  = size(magicped.offspringinfo,1)
offspringfgl = Vector(undef,noff)
peddict = getpeddict(magicped)

subpop = first(keys(peddict))
offls, subped = peddict[subpop]
node_popsize = length(offls)
ped = subped;

popsize_bef = round(Int,node_popsize/select_proportion)
@assert eltype(ped.member) <: Integer
isoogamy = (2 in unique(ped.gender) && length(chrlen)>1)
pedfgl = setfounderfgl(ped,chrlen,isfounderinbred)
nf = ped.nfounder
memberindices = zeros(Int,max(ped.member...))
memberindices[ped.member] = 1:length(ped.member)    


ind = nf+1
motherindex = memberindices[ped.mother[ind]]
fatherindex = memberindices[ped.father[ind]]

mother = (motherindex,rand(1:length(pedfgl[motherindex])))
father = (fatherindex,rand(1:length(pedfgl[fatherindex])))
egg  = crossover.(pedfgl[mother[1]][mother[2]].fgl, isobligate,interference)
sperm  = crossover.(pedfgl[father[1]][father[2]].fgl, isobligate,interference)
zygote = map((x,y)->[x,y],egg,sperm)
if isoogamy
    # if isoogamy, gender 1=female or 2=male; 0=notapplicable is not allowed
    # sex chromsome is in the last
    gender = ped.gender[ind] 
    sexsperm = pedfgl[father[1]][end][gender] 
    zygote[end][2] = sexsperm
end
selpheno = trait_gadd(zygote,seltrait)
simpheno = 0.0
fglpheno = FglPheno(zygote, mother,father, selpheno, simpheno)        
fglpheno
