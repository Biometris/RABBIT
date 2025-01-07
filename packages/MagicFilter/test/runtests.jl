# using Revise
using MagicBase
using MagicFilter
using Test
cd(@__DIR__)

@testset "MagicFilter" begin
    @testset "6star-self4" begin
        include("6star-self4_mix.jl")
    end
end
