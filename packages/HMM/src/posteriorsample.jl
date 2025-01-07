# denote by x_t the hidden state and y_t the obseved data at time t=1...tmax
# fwprob[t] =  Pr(x[t]|y[1:t])
function backwardsample(tranprobseq::AbstractVector{<: Union{Nothing,AbstractMatrix}}, 
    fwprob::AbstractVector{<: AbstractVector};
    samplesize::Integer=1000)    
    tmax = length(fwprob)
    bwsamples= zeros(Int, samplesize, tmax)
    hidden = rand(Categorical(fwprob[end]),samplesize)
    nexthidden = similar(hidden)
    bwsamples[:,end] .= hidden
    for t = tmax-1:-1:1
        sls = unique(hidden)
        if isnothing(tranprobseq[t])
            for s in sls                
                prob = zeros(length(fwprob[t]))             
                prob[s] = 1.0                
                b = hidden .== s
                nexthidden[b] .= rand(Categorical(prob),sum(b))
            end
        else
            for s in sls
                prob = tranprobseq[t][:,s] .* fwprob[t]            
                prob ./= sum(prob)
                b = hidden .== s
                nexthidden[b] .= rand(Categorical(prob),sum(b))
            end
        end 
        hidden .= nexthidden
        bwsamples[:,t] .= hidden
    end
    bwsamples
end

"""
    posteriorsample(initprob,tranprobseq, dataprobseq; samplesize)

randomly sample hidden states for the given HMM

# Arguments

`initprob::AbstractVector`: the probablity vector at initial time t=1

`tranprobseq::AbstractVector{<: Union{Nothing,AbstractMatrix}`: tranprobseq[t] is the transition probaiblity
matrix from t to t+1. The value of nothing denotes an idenity matrix. If length of tranprobseq equals
that of dataprobseq, the last element of tranprobseq is irrelevant. 

`dataprobseq::AbstractVector{<: AbstractVector}`: dataprobseq[t] is the emission probability
vector at time t. 

# Keyword arguments

`samplesize::Integer=1000`: number of random samples

"""
function posteriorsample(initprob::AbstractVector,
    tranprobseq::AbstractVector{<: Union{Nothing,AbstractMatrix}}, 
    dataprobseq::AbstractVector{<: AbstractVector};
    samplesize::Integer=1000)
    fwprob,fwscale = forward(initprob,tranprobseq,dataprobseq)
    loglike = calloglike(fwscale)
    samples = backwardsample(tranprobseq,fwprob; samplesize)
    (loglike=loglike,samples=samples)
end

