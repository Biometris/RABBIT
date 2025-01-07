
using Revise
using MagicBase
cd(@__DIR__)
pwd()

dataid = "multipop"
generate_breedped(; outfile=dataid*".csv",subpopsize=2)
parsebreedped(dataid*".csv"; fixed_nself=10)