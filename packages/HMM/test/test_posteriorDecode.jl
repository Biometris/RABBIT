# https://en.wikipedia.org/wiki/Forward%E2%80%93backward_algorithm
initprob = [0.5, 0.5]
tranprob = [0.7 0.3; 0.3 0.7]
# datamodel M; # row=# possible hidden states; #col= # possible discrete data
# M[i,j]= prob of observing data j given state i
datamodel = [0.9 0.1; 0.2 0.8]
dataseq = [1, 1, 2, 1, 1]
tranprobseq=[tranprob for i =1:size(dataseq,1)-1]
dataprobseq=[datamodel[:, i] for i = dataseq]
tol = 0.0001


fwprob, fwscale = HMM.forward(initprob,tranprobseq,dataprobseq)

begin
    @test fwprob ≈ [[0.8182,0.1818], [0.8834,0.1166], [0.1907, 0.8093],
    [0.7308, 0.2692], [0.8673, 0.1327]] atol =tol
end
@test fwscale ≈ [0.55, 0.6391, 0.3427, 0.4634, 0.6146] atol=tol

loglike,posteriorprob = HMM.posteriordecode(initprob,tranprobseq,dataprobseq)
begin
    @test posteriorprob ≈ [[0.8673, 0.1327], [0.8204, 0.1796], [0.3075, 0.6925],
    [0.8204, 0.1796], [0.8673, 0.1327]] atol = tol
end
