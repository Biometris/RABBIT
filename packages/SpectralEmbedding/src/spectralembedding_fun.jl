
function spectralembedding(similarity::AbstractMatrix; 
    nev::Integer,
    laplaciantype::AbstractString= "randomwalk")
    laplacian = callaplacian(similarity;laplaciantype)
    if size(laplacian,1) < 2*nev
        res = eigen(Matrix(laplacian))
        vals = real.(res.values)
        vecs = real.(res.vectors)
    else
        # default tol = eps(real(eltype(A)))/2
        # eps(Float32)/2: 5.96*10^(-8);  
        # ,tol=6*10^(-8)        
        res = eigs(laplacian; nev,which=:SR,maxiter=10000)
        vals = real.(res[1])
        vecs = real.(res[2])
    end
    (eigenvals = vals, eigenvecs = vecs)
end

function callaplacian(similarity::AbstractMatrix;
    laplaciantype::AbstractString= "randomwalk")
    if laplaciantype == "unnormalized"
        ls = sum(similarity,dims=2)[:,1]
        dd = issparse(similarity) ? spdiagm(ls) : diagm(ls)
        Symmetric(dd - similarity)
    elseif laplaciantype == "randomwalk"
        laplacian_rw(similarity)
    elseif laplaciantype == "symmetric"
        ls = 1 ./ sqrt.(sum(similarity,dims=2)[:,1])
        idd = issparse(similarity) ? spdiagm(ls) : diagm(ls)
        sym = idd*similarity*idd
        Symmetric(I - sym)
    else
        error(string("unknown graph Laplacian: ",laplaciantype))
    end
end

function laplacian_rw(similarity::Matrix{T}) where {T<:Real}
    rw = similarity ./ sum(similarity,dims=2)
    I - rw
end

function laplacian_rw(similarity::SparseMatrixCSC)
    iils,jjls, vvls = findnz(similarity)
    df = DataFrame(i=iils,j=jjls,v=vvls)
    gd = groupby(df,:i)
    df = combine(gd, :j, :v => (x-> x ./ sum(x)) => :vnorm)
    rw = sparse(df[!,:i],df[!,:j],df[!,:vnorm])
    I-rw
end