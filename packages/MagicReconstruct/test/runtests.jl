using MagicBase
using MagicReconstruct
using Test
cd(@__DIR__)

@testset "MagicReconstruct" begin
    @testset "4ril-self3" begin
        include("4ril-self3_array.jl")
        # include("4ril-self3_gbs.jl")
        include("4ril-self3_mix.jl")
    end
    @testset "outbred" begin
        include("4ril-self3_outbred.jl")        
    end
    @testset "multipop" begin
        include("multipop.jl")        
    end
end
