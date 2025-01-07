
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
cd(@__DIR__)
pwd()

designcode = "8ril-self4"
# designcode = "5star-self6"
popsize=100
chrlen = [150,120,80,100,60]
isobligate = false
interference = 0
comment_char = '#'
delim = ','
workdir = pwd()

# @time magicped = formmagicped(designcode,popsize)
# @time plotPedigree(magicped.designinfo)
pedfile="ArabMAGIC_magicped.csv"
@time magicped = readMagicPed(pedfile)
@time offspringfgl = simoffspringfgl(magicped,chrlen,isobligate,interference)
# hcat(offspringfgl[55][1][1]...)
# hcat(offspringfgl[55][1][2]...)
# memb = rand(magicped.offspringinfo[!,:member])
# @time plotPedigree(getsubpedigree(magicped.designinfo,memb))
fhaplofile = "ArabMAGIC_founderhaplo.csv"
@time founderhaplo = formfhaplo(fhaplofile)

parents = magicped.founderinfo[!,:individual]
parents2 = founderhaplo.magicped.founderinfo[!,:individual]
diff = setdiff(parents,parents2)
if !isempty(diff)
    msg = string("parents in ", pedfile, " but not in ",fhaplofile, ":",diff)
    error(msg)
end
fdict = Dict(parents2 .=> 1:length(parents2))
indices = [get(fdict,i,nothing) for i=parents]
@assert parents2[indices] == parents
founderinfo = founderhaplo.magicped.founderinfo[indices,:]
foundergeno = [i[:,indices] for i=founderhaplo.foundergeno]
