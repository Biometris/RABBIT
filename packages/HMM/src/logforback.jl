# denote by x_t the hidden state and y_t the obseved data at time t=1...T*
# logfwprob[t] =  log[Pr(x[t],y[1:t])]
function logforward(initprob::AbstractVector,
    tranprobseq::AbstractVector{<: Union{Nothing,AbstractMatrix}}, 
    dataprobseq::AbstractVector{<: AbstractVector})    
    tmax = length(dataprobseq)
    T = eltype(initprob)
    logfwprob = Vector{Vector{T}}(undef,tmax)
    logfwprob[1] = log.(initprob .* dataprobseq[1])
    for t=2:tmax
        logpmax  = max(logfwprob[t-1]...)
        prob  = exp.(logfwprob[t-1] .- logpmax)
        if isnothing(tranprobseq[t-1])
            prob2 = prob .* dataprobseq[t]
        else
            prob2 = (tranprobseq[t-1]' * prob) .* dataprobseq[t]
        end
        logfwprob[t] = log.(prob2) .+ logpmax
    end
    logfwprob
end

# denote by x_t the hidden state and y_t the obseved data at time t=1...tmax
# logbwprob[tmax] = backinit (=0 by default)
# logbwprob[t] =  log[Pr(y[t+1:T]|x_t)]
function logbackward(tranprobseq::AbstractVector{<: Union{Nothing,AbstractMatrix}}, 
    dataprobseq::AbstractVector{<: AbstractVector};
    backinit::Union{Nothing,AbstractVector}=nothing)
    tmax = length(dataprobseq)    
    if isnothing(backinit) 
        logbwprob = Vector{Vector{Float64}}(undef,tmax)        
    else
        T = eltype(backinit)
        logbwprob = Vector{Vector{T}}(undef,tmax)        
    end
    logbwprob[end] = isnothing(backinit) ? zeros(length(last(dataprobseq))) : backinit
    for t=tmax-1:-1:1
        logpmax  = max(logbwprob[t+1]...)
        prob  = exp.(logbwprob[t+1] .- logpmax)
        if isnothing(tranprobseq[t])
            prob2 = dataprobseq[t+1] .* prob
        else
            prob2 = tranprobseq[t] * (dataprobseq[t+1] .* prob)
        end
        logbwprob[t] = log.(prob2) .+ logpmax
    end
    logbwprob
end

calloglike(fwscale::AbstractVector) = sum(log.(fwscale))

function caloglike(logfwprob::AbstractVector{<: AbstractVector},
    logbwprob::AbstractVector{<: AbstractVector})
    tmax = length(logbwprob)
    [begin
        fw=logfwprob[t]
        fwmax = max(fw...)
        fw .-= fwmax
        bw=logbwprob[t]
        bwmax = max(bw...)
        bw .-= bwmax
        log(sum(exp.(fw) .* exp.(bw))) + bwmax + fwmax
    end for t=1:tmax]
end

function caloglike(fwprob::AbstractVector{<: AbstractVector},
    fwscale::AbstractVector,
    logbwprob::AbstractVector{<: AbstractVector})
    alike = accumulate(+,log.(fwscale))
    tmax = length(logbwprob)
    clike = [begin
        bw=logbwprob[t]
        pmax = max(bw...)
        log(sum(exp.((bw .- pmax)) .* fwprob[t])) + pmax
    end for t=1:tmax]
    alike .+ clike
end
