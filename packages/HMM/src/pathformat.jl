function splitindex(f::Function,A::AbstractVector)
    size(A,1)==1 && return [1:1]
    res=Vector{typeof(1:1)}()
    i0=1
    for i=2:size(A,1)
        if !f(A[i-1],A[i])
            push!(res,i0:i-1)
            i0=i
        end
    end
    push!(res,i0:size(A,1))
    res
end

function splitindex(A::AbstractVector)
    f(x,y)= x == y
    splitindex(f,A)
end

function tojumppath(dpath::AbstractVector)
    ii=splitindex(dpath)
    ii2=first.(ii)
    res=[ii2 dpath[ii2]]'
    hcat(res,[length(dpath)+1 missing]')
end

function todiscretepath(jpath::AbstractMatrix)
    nrow,ncol = size(jpath)
    nrow == 2 || @error("jump path must be a matrix with 2 rows")
    vcat([repeat([jpath[2,i]],jpath[1,i+1]-jpath[1,i]) for i=1:ncol-1]...)
end

function tostringpath(jpath::AbstractMatrix;delim::AbstractChar='-')
    nrow,ncol = size(jpath)
    nrow == 2 || @error("jump path must be a matrix with 2 rows")
    join(reshape(jpath,:)[1:end-1],delim)
end

function tojumppath(spath::AbstractString;delim::AbstractChar='-')
    ss=split(spath,delim)
    n=length(ss)
    n>=3 || @error("at least two delimiters are required")
    isodd(n) || @error("number of delimiters must be even")
    res=tryparse.(Int,ss)
    reshape(vcat(res,missing),2,:)
end
