using Revise
using MagicBase
cd(@__DIR__)
pwd()

magicped = MagicBase.generate_magicped(;
    designcodes = ["bc3-self2","P2/P3=>DH","P4/3/P2/P3//P5/P6=>3","2ril-self6"],
    founders = ["P1||P2","NA", "NA", "P5||P6"],
    subpopsizes =10*ones(Int, 4));
plotmagicped(magicped)

# savemagicped("ped.csv",magicped)
