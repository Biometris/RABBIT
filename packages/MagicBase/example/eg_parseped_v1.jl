using Revise
using MagicBase
using Pedigrees
cd(@__DIR__)
pwd()

using CSV, DataFrames

magicped = parsebreedped("example_breedped.csv")
plotmagicped(magicped)



breedcode = "p1/5/p2/3/p2/p3//p2/p4/4/p2 => 0"


df = MagicBase.parsebreedcode(breedcode)
ped = Pedigree(df)
plotped(ped)


MagicBase.generate_breedped(outfile="breedped.csv")

parsebreedped("breedped.csv")


magicped = readmagicped("breedped_magicped.csv")
split_subpop("breedped_magicped.csv")

