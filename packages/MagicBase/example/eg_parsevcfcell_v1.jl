using Revise
using MagicBase
cd(@__DIR__)
pwd()

formatpriority = ["AD","GT"]
# formatpriority = ["AD"]
format0 = "GT:AD"
format = Symbol.(split(format0,":"))
prioritytuple = (; zip(Symbol.(formatpriority), 1:length(formatpriority))...)
formatcode = [get(prioritytuple,i,-1) for i in format] 
rowgeno = ["0/0:3,4","1/1:0,4","./1:10,0"]
ncell = length(rowgeno)
formatls = Vector(undef,ncell)
resgeno = Vector(undef,ncell)
missingset = [".","./.", ".|."]
keepvcf = false


for i in 1:ncell
    @inbounds resgeno[i], formatls[i] = MagicBase.parse_vcfcell(rowgeno[i], formatcode,
        formatpriority, missingset,keepvcf)
end
