using Revise
using HMM
cd(@__DIR__)
pwd()

# https://www.cs.sjsu.edu/~stamp/RUA/HMM.pdf
# Mark Stamp. A revealing introduction to hidden Markov model. 2001. 

# state space of hidden temperature: {high, low}
initprob = [0.6, 0.4]
tranprob = [0.7 0.3; 0.4 0.6]
# state space of observed size of tree growth ring: {small, medium,large} or {1,2,3}
dataseq = [1,2,1,3]   
emissprob = [0.1 0.4 0.5; 0.7 0.2 0.1]  # emissprob[i,j]=Pr(ringsize=j|temperature=i)


#transform parameters
tmax = length(dataseq)
tranprobseq = repeat([tranprob],tmax-1)
dataprobseq=[emissprob[:, j] for j in dataseq]


fwprob,fwscale = HMM.forward(initprob,tranprobseq,dataprobseq)
loglike = HMM.calloglike(fwscale)
bwprob=HMM.backward(tranprobseq,dataprobseq,fwscale)

loglike,postprob = HMM.posteriordecode(initprob,tranprobseq,dataprobseq)

using DataFrames
DataFrame(Matrix(reduce(hcat,fwprob)'),:auto) 
DataFrame(Matrix(reduce(hcat,bwprob)'),:auto) 
DataFrame(Matrix(reduce(hcat,postprob)'),:auto) 

reduce(*,fwscale) ≈ exp(loglike)
round.(fwscale,digits=6)


logprob, path,his_logprob,his_index  = HMM.viterbidecode(initprob,tranprobseq,dataprobseq)

exp(logprob)
path
exp.(reduce(hcat,his_logprob))