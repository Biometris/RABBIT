using Revise
using HMM
cd(@__DIR__)


initprob = [0.6, 0.4]
tranprob = [0.7 0.3; 0.4 0.6]
dataseq = [1, 2, 1, 3]
datamodel = [0.1 0.4  0.5; 0.7 0.2 0.1]
tranprobseq=[tranprob for i =1:size(dataseq,1)-1]
dataprobseq=[datamodel[:, i] for i = dataseq]
@time logvpathprob, vpath = HMM.viterbidecode(initprob,tranprobseq,dataprobseq)

vpath


dpath = [1,6,1,1,1,16,16,6, 16,16,11,16]
spath = HMM.tostringpath(HMM.tojumppath(dpath))
spath == "1-1-2-6-3-1-6-16-8-6-9-16-11-11-12-16-13"
HMM.todiscretepath(HMM.tojumppath(spath)) == dpath
