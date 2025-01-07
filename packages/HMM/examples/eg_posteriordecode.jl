using Revise
using HMM
using LinearAlgebra
cd(@__DIR__)


# https://en.wikipedia.org/wiki/Forward%E2%80%93backward_algorithm
initprob = [0.5, 0.5]
tranprob = [0.7 0.3; 0.3 0.7]
# datamodel M; # row=# possible hidden states; #col= # possible discrete data
# M[i,j]= prob of observing data j given state i
datamodel = [0.9 0.1; 0.2 0.8]
dataseq = [1, 1, 2, 1, 1]
tranprobseq=[tranprob for i =1:size(dataseq,1)-1]
dataprobseq=[datamodel[:, i] for i = dataseq]

tranprobseq2 = Vector{Union{Nothing,Matrix{Float64}}}(undef,4)
tranprobseq2 .= Matrix{Float64}.(tranprobseq)
# tranprobseq2[3] = nothing
eltype(tranprobseq2)


res = HMM.forward(initprob,tranprobseq,dataprobseq)
keys(res)
fwscale=res.fwscale
bwprob=HMM.backward(tranprobseq,dataprobseq,fwscale)


tranprobseq = Matrix{Float64}.(tranprobseq)
loglike,posteriorprob = HMM.posteriordecode(Vector{Float32}(initprob),tranprobseq2,
    dataprobseq
)


fwprob,fwscale = HMM.forward(initprob,tranprobseq2,dataprobseq)

logbwprob=HMM.logbackward(tranprobseq2,dataprobseq)
res = HMM.caloglike(fwprob,fwscale,logbwprob)
all(res .≈ loglike)

logfwprob=HMM.logforward(initprob, tranprobseq,dataprobseq)
res = HMM.caloglike(logfwprob,logbwprob)
all(res .≈ loglike)
