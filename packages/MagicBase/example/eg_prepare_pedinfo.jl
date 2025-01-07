using MagicBase, Pedigrees
designinfo = MagicBase.parsedesign("P1/P2=>DH")
plotped(designinfo.pedigree)


designinfo = MagicBase.parsedesign("P1/4/P1/3/P1/P2//P3/P4=>2")
plotped(designinfo.pedigree)


designinfo = MagicBase.parsedesign("P1/P2//P3/P4/3/P5/P6//P7/P8/4/P1/P2//P3/P4/3/P5/P6//P7/P8=>2"; popid="pop")
designinfo = MagicBase.parsedesign("P1/P2//P3/P4/3/P5/P6//P7/P8=>2")
plotped(designinfo.pedigree)





designinfo = MagicBase.parsedesign("P1/P2//P3/P4/3/P1/P2//P3/P4=>2";popid="pop")
plotped(designinfo.pedigree)


breedcode = "P1/P2//P3/P4/3/P5/P6//P7/P8/4/P1/P2//P3/P4/3/P5/P6//P7/P8=>2"
breedcode = "P1/P2//P3/P4/3/P1/P2//P3/P4=>2"
df = MagicBase.parsebreedcode(breedcode)
ped = Pedigree(unique(df))

plotped(ped)




designinfo = parsedesign("nfounder=3||ibd=1.0||mapexpansion=5.0")


using MagicBase
cd(@__DIR__)
magicped = generate_magicped(;
  designcodes=["P1/P2=>DH", "2ril-self3", "ibd=1.0||mapexpansion=5.0","P4/3/P4//P2/P3=>3"],
  founders = ["NA","P1||P3","P3||P4||P5","NA"],
  subpopsizes=100*ones(4)
)
savemagicped("example_ped_junc.csv", magicped)
plotmagicped(magicped)


MagicBase.pedfile_designcode2ped("example_ped_junc.csv"; outfile = "example_ped.csv")