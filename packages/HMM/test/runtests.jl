# using Revise
# using Pkg
# Pkg.activate(joinpath(@__DIR__,".."))
using HMM
using Test

tests = ["test_posteriorDecode.jl","test_posteriorsample.jl",
  "test_viterbiDecode.jl", "test_logforback.jl"]    # the test file names are stored as strings...
for t in tests
  println("testing  ",t)
  include(t)             # ... so that they can be evaluated in a loop
end
