using MagicPrior
using Pedigrees
using Test
cd(@__DIR__)

@time @testset "MagicPrior" begin
    for pedfile = ["F1.csv","8waymagic.csv","CC_ped.csv"]
        ped = Pedigrees.readped(pedfile)
        nfounder = ped.nfounder
        isautols = issubset(["female","male"], unique(ped.gender)) ? [true, false] : [true]
        indarray = length(ped.member)>=ped.nfounder+2 ? ped.member[end-1:end] : ped.member[end:end]        
        for isfounderinbred=[true,false],isautosome=isautols,isfglexch=[true,false]
            @test isa(MagicPrior.magicorigin(ped; 
                memberlist=indarray, isfounderinbred,
                isautosome=isautosome,isfglexch=isfglexch,isconcise=false),Matrix)
            @test isa(MagicPrior.magicprior(ped; 
                memberlist=indarray, isfounderinbred,
                isautosome,isfglexch,isconcise=false),Dict)
        end
    end
end
