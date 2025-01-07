# https://en.wikipedia.org/wiki/Forward%E2%80%93backward_algorithm
initprob = [0.5, 0.5]
tranprob = [0.7 0.3; 0.3 0.7]
# datamodel M; # row=# possible hidden states; #col= # possible discrete data
# M[i,j]= prob of observing data j given state i
datamodel = [0.9 0.1; 0.2 0.8]
dataseq = [1, 1, 2, 1, 1]
tranprobseq=[tranprob for i =1:size(dataseq,1)-1]
dataprobseq=[datamodel[:, i] for i = dataseq]

fwprob,fwscale = HMM.forward(initprob,tranprobseq,dataprobseq)
logbwprob=HMM.logbackward(tranprobseq,dataprobseq)

res = HMM.caloglike(fwprob,fwscale,logbwprob)
@test all(res .≈ loglike)

logfwprob=HMM.logforward(initprob, tranprobseq,dataprobseq)
res = HMM.caloglike(logfwprob,logbwprob)
@test all(res .≈ loglike)
