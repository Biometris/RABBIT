# using Revise
using MagicMap
using MagicBase
using Test
cd(@__DIR__)

@testset "MagicMap" begin
    @testset "4ril-self3" begin
        @time include("4ril-self3_array.jl")
        @time include("4ril-self3_mix.jl")
    end
    @testset "outbred" begin
        @time include("4star-self1_outbred.jl")        
    end
    @testset "multipop" begin
        @time include("multipop.jl")        
    end
end
