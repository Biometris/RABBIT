# using Revise
# using Pkg
# Pkg.activate(joinpath(@__DIR__,".."))
using Pedigrees
using Test
cd(@__DIR__)

@time @testset "Pedigrees" begin
    for pedfile = ["F1.csv","8waymagic.csv","CC_ped.csv"]
        @test isa(Pedigrees.readped(pedfile),Pedigree)
        ped=Pedigrees.readped(pedfile)
        @test isa(Pedigrees.plotped(ped),Any)
        @test isa(Pedigrees.saveped("test.csv",ped),String)
        founders = [string("P",i) for i=1:ped.nfounder]
        @test isa(Pedigrees.setfounderid(ped,founders),Pedigree)
        @test isa(Pedigrees.getsubped(ped,ped.member[end]),Pedigree)
        @test isa(Pedigrees.toindexped(ped),Pedigree)
        iiped = Pedigrees.toindexped(ped)
        @test isa(Pedigrees.getsubped(iiped,iiped.member[end]),Pedigree)
    end
    rm("test.csv")
end
