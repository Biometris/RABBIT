
# denote by x_t the hidden state and y_t the obseved data at time t=1...tmax
# fwprob[t] =  Pr(x[t]|y[1:t])
# fwscale[t] = Pr(y[t]|y[1:t-1])
function forward(initprob::AbstractVector,    
    tranprobseq::AbstractVector{<: Union{Nothing,AbstractMatrix}}, 
    dataprobseq::AbstractVector{<: AbstractVector})
    tmax = length(dataprobseq)    
    T = eltype(initprob)
    fwprob=Vector{Vector{T}}(undef,tmax)
    fwscale=Vector{T}(undef,tmax)
    prob=initprob .* dataprobseq[1]
    scale=sum(prob)
    prob ./= scale
    fwscale[1] = scale
    fwprob[1] = prob
    for t = 2:tmax
        if isnothing(tranprobseq[t-1])
            prob = fwprob[t-1] .* dataprobseq[t]
        else
            prob = (tranprobseq[t-1]' * fwprob[t-1]) .* dataprobseq[t]
        end
        scale=sum(prob)
        prob ./= scale
        fwscale[t] = scale
        fwprob[t] = prob
    end
    (fwprob=fwprob,fwscale=fwscale)
end

# denote by x_t the hidden state and y_t the obseved data at time t=1...tmax
# bwprob[t] =  Pr(y[t+1:T]|x[t]) / Pr(y[t+1:T]|y[1:t])
function backward(tranprobseq::AbstractVector{<: Union{Nothing,AbstractMatrix}}, 
    dataprobseq::AbstractVector{<: AbstractVector},
    fwscale::AbstractVector)    
    tmax = length(dataprobseq)
    T = eltype(fwscale)
    bwprob=Vector{Vector{T}}(undef,tmax)
    bwprob[end] = ones(T, length(first(dataprobseq)))
    for t = tmax-1:-1:1
        if isnothing(tranprobseq[t])        
            prob = dataprobseq[t+1] .* bwprob[t+1]
        else
            prob = tranprobseq[t] * (dataprobseq[t+1] .* bwprob[t+1])
        end
        prob ./= fwscale[t+1]
        bwprob[t] = prob
    end
    bwprob
end

"""
    posteriordecode(initprob,tranprobseq, dataprobseq)

calculate marginal posterior probablities using the HMM posterior decoding
algorightm.

# Arguments

`initprob::AbstractVector`: the probablity vector at initial time t=1

`tranprobseq::AbstractVector{<: Union{Nothing,AbstractMatrix}`: tranprobseq[t] is the transition probaiblity
matrix from t to t+1. The value of nothing denotes an idenity matrix. If length of tranprobseq equals
that of dataprobseq, the last element of tranprobseq is irrelevant. 

`dataprobseq::AbstractVector{<: AbstractVector}`: dataprobseq[t] is the emission probability
vector at time t. 

`digits::Union{Nothing,Integer}=nothing`: if nothing, not round the probabilities, and otherwise round the probabilities 
to the specified number of digits after the decimal place. 

"""
function posteriordecode(initprob::AbstractVector,
    tranprobseq::AbstractVector{<: Union{Nothing,AbstractMatrix}}, 
    dataprobseq::AbstractVector{<: AbstractVector};
    digits::Union{Nothing,Integer}=nothing)
    fwprob,fwscale = forward(initprob,tranprobseq,dataprobseq)
    loglike = calloglike(fwscale)
    bwprob=backward(tranprobseq,dataprobseq,fwscale)
    if isnothing(digits)
        posteriorprob=[fwprob[i] .* bwprob[i] for i in 1:length(dataprobseq)]
    else
        posteriorprob=[round.(fwprob[i] .* bwprob[i]; digits) for i in 1:length(dataprobseq)]
    end
    (loglike=loglike,posteriorprob=posteriorprob)
end

function calloglike(initprob::AbstractVector,    
    tranprobseq::AbstractVector{<: Union{Nothing,AbstractMatrix}}, 
    dataprobseq::AbstractVector{<: AbstractVector})
    fwscale = last(forward(initprob,tranprobseq,dataprobseq))
    loglike = calloglike(fwscale)
    loglike
end