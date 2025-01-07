# https://en.wikipedia.org/wiki/Forward%E2%80%93backward_algorithm
initprob = [0.5, 0.5]
tranprob = [0.7 0.3; 0.3 0.7]
# datamodel M; # row=# possible hidden states; #col= # possible discrete data
# M[i,j]= prob of observing data j given state i
datamodel = [0.9 0.1; 0.2 0.8]
dataseq = [1, 1, 2, 1, 1]
tranprobseq=[tranprob for i =1:size(dataseq,1)-1]
dataprobseq=[datamodel[:, i] for i = dataseq]

loglike,posteriorprob = HMM.posteriordecode(initprob,tranprobseq,dataprobseq)
loglike2, samples = HMM.posteriorsample(initprob, tranprobseq, dataprobseq; samplesize=1000)
probls = first.(posteriorprob)
probls2 = [sum(i .== 1)/length(i)  for i in eachcol(samples)]
@test  all(isapprox.(probls,probls; atol = 0.01))
