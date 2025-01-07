
# Example from Stamp2004_A revealing introduction to hidden Markov models.pdf
initprob = [0.6, 0.4]
tranprob = [0.7 0.3; 0.4 0.6]
dataseq = [1, 2, 1, 3]
# datamodel M; # row=# possible hidden states; #col= # possible discrete data
# M[i,j]= prob of observing data j given state i
datamodel = [0.1 0.4  0.5; 0.7 0.2 0.1]
tranprobseq=[tranprob for i =1:size(dataseq,1)-1]
dataprobseq=[datamodel[:, i] for i = dataseq]
tol = 0.0001

logvpathprob, vpath = HMM.viterbidecode(initprob,tranprobseq,dataprobseq)
@test logvpathprob ≈ -5.8702 atol = tol
@test vpath == [2,2,2,1]

dpath = [1,6,1,1,1,16,16,6, 16,16,11,16]
spath = HMM.tostringpath(HMM.tojumppath(dpath))
@test spath == "1-1-2-6-3-1-6-16-8-6-9-16-11-11-12-16-13"
@test HMM.todiscretepath(HMM.tojumppath(spath)) == dpath
