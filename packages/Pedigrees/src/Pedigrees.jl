module Pedigrees

using DataFrames, CSV
using GraphRecipes, Plots
using Random: randstring

export
    # private
    readped, saveped, setfounderid, getsubped, plotped,
    Pedigree

include("pedigreeui.jl")

end # module
