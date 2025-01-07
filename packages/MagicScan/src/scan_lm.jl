
function permute_threshold(yresp::AbstractVector,
    xpred::AbstractMatrix,
    haploprob::AbstractVector,
    pop2off::AbstractDict;
    nperm::Integer=200)
    res = zeros(nperm,3)    
    y = similar(yresp)
    x = similar(xpred)        
    for i in 1:nperm
        permute_resp_pred!(y, x, yresp, xpred, pop2off)
        loglikels,lodls,log10pls = scan_lm(y,x,haploprob);
        res[i,1] = max(reduce(vcat,loglikels)...)
        res[i,2] = max(reduce(vcat,lodls)...)
        res[i,3] = max(reduce(vcat,log10pls)...)
    end
    threshoulddf=DataFrame(res, [:loglike,:lod,:neg_log10p])
    peremuteid=[string("permute",i) for i in 1:nperm]
    insertcols!(threshoulddf,1,:peremute=>peremuteid)
    threshoulddf
end

function permute_resp_pred!(y::AbstractVector,x::AbstractMatrix,
    yresp::AbstractVector,xpred::AbstractMatrix,pop2off::AbstractDict)    
    for (_,offs) in pop2off
        offs2 = sample(offs, length(offs),replace=false)
        y[offs] .= yresp[offs2]
        x[offs,:] .= xpred[offs2,:]
    end
    y,x
end

function scan_lm(yresp::AbstractVector,
    xpred::AbstractMatrix,
    haploprob::AbstractVector)
    nonmiss = .!ismissing.(yresp)
    y = yresp[nonmiss]
    x = xpred[nonmiss,:]
    beta0, _,loglike0= fit_lm(y,x)
    nchr = length(haploprob)    
    loglikels = Vector(undef,nchr)
    lodls = Vector(undef,nchr)
    log10pls = Vector(undef,nchr)
    nthread = Threads.nthreads()
    for chr in 1:nchr
        chrprob = haploprob[chr]
        nsnp = length(chrprob)
        loglikels_chr = zeros(nsnp)
        lodls_chr = zeros(nsnp)
        log10pls_chr = zeros(nsnp)
        basesize = max(1, div(nsnp, 5*nthread))
        # ThreadsX.foreach(eachindex(chrprob,loglikels_chr,lodls_chr,log10pls_chr); basesize) do snp            
        for snp in 1:nsnp
            gpred = hcat(x,chrprob[snp][nonmiss,2:end])
            beta, _,loglike= fit_lm(y,gpred)
            df = length(beta) - length(beta0)            
            log_ratio = -2*(loglike0 - loglike)
            loglikels_chr[snp]  = loglike
            lodls_chr[snp]  = log_ratio * log10(ℯ)/2
            # println("snp=",snp,", df=",df,", beta=",beta, ",beta0=",beta0,",loglike=",loglike,
            #     ",loglike0=",loglike0,",df=",df)
            # The error occurs for -logccdf(Chisq(103),1800)
            # The error does not occurs for -logccdf(Chisq(104),1800)
            # ERROR: ArgumentError: Unsupported order |ν| > 50 off the positive real axis if df > 102 and log_ratio 1999
            if log_ratio > 1700
                try 
                    temp = -logccdf(Chisq(df),log_ratio) 
                catch err
                    @warn string("warning: ",err, ", df=",df, ", log_ratio=",log_ratio) maxlog=20
                    @info string("recaculate -logccdf(Chisq(df),log_ratio) with log_ratio = 1700") maxlog=20
                    temp = -logccdf(Chisq(df),1700.0) 
                end
            else
                temp = -logccdf(Chisq(df),log_ratio) 
            end 
            log10pls_chr[snp] = temp * log10(ℯ)    
        end
        loglikels[chr] = loglikels_chr
        lodls[chr] = lodls_chr
        log10pls[chr] = log10pls_chr
    end 
    loglikels,lodls, log10pls
end




# https://en.wikipedia.org/wiki/Bayesian_linear_regression
# y = X*β+ϵ, ϵ ~ N[0, σ^2]; β_i  ~ N[0,1/precision]
function fit_lm(y::AbstractVector, X::AbstractMatrix)    
    precision = 0.01
    n,k = size(X)
    beta = (X'X .+ precision*I(k))\ X'y
    ss = sum((y - X*beta).^2)
    sig = ss / n # variance σ^2 of noise
    loglike = -0.5*n*log(sig) - 0.5*ss/sig
    beta, sig, loglike
end

# magicancestry.haploprob[chr][offspring][marker, haplostate]
# return haploprob[chr][marker][offspring, haplostate]
# function get_haploprob(magicancestry::MagicAncestry,offls::AbstractVector)
#     haploprob = [begin             
#         nsnp = size(first(prob),1)
#         [reduce(hcat,[prob[i][j,:] for i in offls])' for j in 1:nsnp]        
#     end for prob in magicancestry.haploprob];
#     haploprob
# end


function get_haploprob(magicancestry::MagicAncestry,offls::AbstractVector)
    [begin 
        chrprob = magicancestry.haploprob[chr][offls];
        noffspring = length(chrprob)
        nsnp, nstate = size(chrprob[1])
        issparseprob = issparse(chrprob[1])
        nzcolls = [findall(sum(chrprob[off],dims=1)[1,:] .> 0 ) for off in 1:noffspring]
        [begin
            res = issparseprob ? spzeros(noffspring, nstate) : zeros(noffspring, nstate)
            for off in 1:noffspring
                nzcol = nzcolls[off]    
                res[off, nzcol] .= chrprob[off][t, nzcol]
            end
            res
        end for t in 1:nsnp]        
    end for chr in eachindex(magicancestry.haploprob)]
end

