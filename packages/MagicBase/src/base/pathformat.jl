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
