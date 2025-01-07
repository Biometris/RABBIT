
function knn_spectralordering(similarity::AbstractMatrix;
    laplaciantype::AbstractString= "randomwalk",
    knnmin::Integer = round(Int,sqrt(size(similarity,1))),
    mincomponentsize::Integer=1,
    alwayskeep::Union{Nothing,Real}=nothing)
    knn = first(findknn(similarity; knnmin,mincomponentsize,
        alwayskeep,ncomponent=1,io=nothing, verbose=false))
    similarity2=toknnsimilarity(similarity,knn;alwayskeep)
    similarity2, connectednodes, _ = takecomponents(similarity2; ncomponent=1)
    oo = spectralordering(similarity2; laplaciantype)
    connectednodes[oo]
end

function spectralordering(similarity::AbstractMatrix;
    laplaciantype::AbstractString= "randomwalk")
    ccsizels = length.(get_connected_components(similarity))
    nev = 1+length(ccsizels)
    laplacian = callaplacian(similarity;laplaciantype)
    if size(laplacian,1) ==1
        return [1]
    elseif size(laplacian,1) < 200
        res = eigen(Matrix(laplacian))
        fiedler_val = real(res.values[nev])
        fiedler_vec = real.(res.vectors[:,nev])        
    else        
        # https://discourse.julialang.org/t/finding-eigenvalues-with-the-smallest-magnitude-using-eigs/79130
        # to replace Arpack.jl by KrylovKit.jl
        res = eigs(laplacian; nev,which=:SR,maxiter=10000)
        min_eval = real(res[1][nev-1])
        fiedler_val = real(res[1][nev])
        # @info string("eigenvalues = ",[min_eval,fiedler_val])
        if abs(min_eval )> 1e-5
            msg = string("eigenvalues = ",[min_eval,fiedler_val],", the smallest is not close to zero")                            
            @warn msg
        end        
        if fiedler_val < 10^(-10)    
            nev += 1                    
            res = eigs(laplacian; nev,which=:SR,maxiter=10000)          
            fiedler_vec = real.(res[2][:,nev])
        else
            # res = eigs(laplacian; nev=2,which=:SR)
            fiedler_vec = real.(res[2][:,nev])
        end  
        if nev > 2
            msg = string("similarity size=",size(similarity),
                ",connected component sizes=",ccsizels)
            msg *= string(", eigenvalues = ", real.(res[1]), 
            ", ordering based the ", nev, "-th eigenvector")
            @warn msg 
        end
    end
    sortperm(fiedler_vec)
end