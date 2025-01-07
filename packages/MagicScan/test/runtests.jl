
using MagicScan
using MagicBase
using Test
cd(@__DIR__)

@testset "MagicScan" begin
    @testset "4ril-self6" begin
        @time include("4ril-self6.jl")
    end
end
