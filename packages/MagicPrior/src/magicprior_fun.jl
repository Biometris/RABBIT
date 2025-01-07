
"""
    magicprior(pedigree::Pedigree,founderfgl::AbstractMatrix; kwargs...)

calculate prior distribution of recombination breakpoints in pedigree. 

For each pedigree member in `memberlist`, return continutuous Markov process
parameter values: initial probability vector and rate matrix under three models
"depmodel", "indepmodel", and "jointmodel", a relationship between the two
ancestral processes along each of two homologous chromosomes.

# Position arguments 

`pedigree::Pedigree`: pedigree struct. 

`isfounderinbred::Bool`: if true, founders are inbred, and otherwise outbred.

# Keyword Arguments

`memberlist::Union{Nothing,AbstractVector}=nothing`: a list of pedigree members.

`isautosome::Bool=true`: autosome rather than sex chromosome,

`isfglexch::Bool=false`: not assume that  founder genotype lables (FGLs) are exchangeable.

`isconcise::Bool=false`: if isconcise is true, the parameter values for
"jointmodel" are not calculated.

"""
function magicprior(pedigree::Pedigree;
    memberlist::Union{Nothing,AbstractVector}=nothing,
    isfounderinbred::Bool=true,
    isautosome::Bool=true,
    isfglexch::Bool=false,
    isconcise::Bool=false)
    if isnothing(memberlist)
        memberlist=[pedigree.member[end]]
    else
        b = [!in(i,pedigree.member) for i in memberlist]
        if any(b)
            @error string(memberlist[b], " in memberlist are not members of pedigree")
        end
    end
    if isfglexch
        founderfgl = get_founderfgl(pedigree.nfounder;isfounderinbred)
        res = identityprior(pedigree,founderfgl; memberlist,isautosome,isconcise)        
    else
        # Too slow: res = ancestryprior(pedigree,founderfgl; memberlist,isautosome,isconcise)         
        # subped is much smaller than pedigree; then tranfer from statespace of subped into original statespace of pedigree
        allfounders = pedigree.member[1:pedigree.nfounder]
        resls = [begin             
            member2 = convert(String,member)
            subped = Pedigrees.getsubped(pedigree,member2)
            res = magicsubprior(subped,allfounders; member = member2, isfounderinbred, isautosome,isconcise)             
            res 
        end for member in memberlist]
        res = reduce(merge,resls)
    end
    dememoize()
    res
end

function magicsubprior(subped::Pedigree, allfounders::AbstractVector;
    member::T where T<:Union{Integer, String},
    isfounderinbred::Bool=true,
    isautosome::Bool=true,
    isfglexch::Bool=false,
    isconcise::Bool)    
    nfounder = length(allfounders)
    founderfgl = get_founderfgl(nfounder;isfounderinbred)
    # mapping subfoundeers to allfounders    
    founderdict = Dict(allfounders .=> 1:length(allfounders))
    subfounders = subped.member[1:subped.nfounder]    
    subfounder_indices = [founderdict[i] for i in subfounders]
    subfounderfgl = founderfgl[subfounder_indices,:]
    if isfglexch
        prior = identityprior(subped,subfounderfgl; memberlist = [member],isautosome,isconcise)    
    else
        prior = ancestryprior(subped,subfounderfgl; memberlist = [member],isautosome,isconcise)    
    end
    newprior = expand_prior_statespace(prior,founderfgl,subfounder_indices)
    dememoize()
    newprior
end

function get_founderfgl(nfounder::Integer; isfounderinbred::Bool=true)
    if isfounderinbred
        founderfgl= [i for i = 1:nfounder, j = 1:2]
    else
        founderfgl= [2*(i-1)+j for i = 1:nfounder, j = 1:2]
    end
end


function expand_prior_statespace(prior::AbstractDict,founderfgl::AbstractMatrix, subfounder_indices::AbstractVector)
    memid, modells= only(prior)
    modelnames= propertynames(modells)
    nfgl = length(get_fglset(founderfgl))
    subfgl_indices = get_subfgl_indices(founderfgl,subfounder_indices)    
    substate_indices = get_substate_indices(nfgl, subfgl_indices)
    newmodells = [if in(modelname,[:depmodel, :indepmodel, :jointmodel])
        nstate = modelname == :depmodel ? nfgl : nfgl^2
        indices = modelname == :depmodel ? subfgl_indices : substate_indices
        modelval = [expand_sparse_vecormat(i,nstate, indices) for i in modells[modelname]]
        (modelname,modelval)
    elseif in(modelname,[:contmarkov])
        # markovls=[[initm,ratem],[initp,ratep],[initmp,ratemp]] if concise = false
        # markovls=[[initm,ratem],[initp,ratep]] if concise = true
        markovls = modells[modelname]
        modelval = [begin 
            nstate = i<=2 ? nfgl : nfgl^2
            indices = i<=2 ? subfgl_indices : substate_indices
            [expand_sparse_vecormat(a,nstate, indices) for a in markovls[i]]
        end for i in eachindex(markovls)]
        (modelname,modelval)
    else
        @error string("unexpected modelname=",modelname)
    end for modelname in modelnames]
    newprior = Dict(memid => (;newmodells...))
    newprior
end

function get_subfgl_indices(founderfgl::AbstractMatrix, subfounder_indices::AbstractVector)
    fglset= get_fglset(founderfgl)
    fgldict = Dict(fglset .=> 1:length(fglset))
    subfounderfgl = founderfgl[subfounder_indices,:]
    subfglset = get_fglset(subfounderfgl)
    subfgl_indices = [fgldict[i] for i in subfglset]
    subfgl_indices
end

function get_substate_indices(nfgl::Integer, subfgl_indices::AbstractVector) 
    states = [[i,j] for i in 1:nfgl for j in 1:nfgl]
    statedict = Dict(states .=> 1:length(states))
    substates = [[i,j] for i in subfgl_indices for j in subfgl_indices] 
    substate_indices = [statedict[i] for i in substates]
    substate_indices
end

function expand_sparse_vecormat(a::AbstractVecOrMat, nstate::Integer,indices::AbstractVector)
    if isa(a,AbstractVector)
        let vec = spzeros(eltype(a),nstate)
            vec[indices] .= a
            vec
        end
    else
        let mtx = spzeros(eltype(a),nstate,nstate)
            mtx[indices,indices] .= a
            mtx
        end
    end
end

function ancestryprior(ped::Pedigree,founderfgl::AbstractMatrix;
    memberlist::AbstractVector{T} where T<:Union{Integer, String},
    isautosome::Bool=true,
    isconcise::Bool=false)
    ancestry=pedancestry(ped,founderfgl,memberlist,isautosome,isconcise)
    Dict(begin
        # col 2:7= "fi^m", "fi^p","phiij^mp", "Rij^m", "Rij^p"
        initm=ancestry[i+1,2]
        initp=ancestry[i+1,3]
        # row-major order = A'[:], column-major order = A[:], where A is a matrix.
        initmp=(ancestry[i+1,4])'[:]
        densitym=ancestry[i+1,5]
        densityp=ancestry[i+1,6]
        ratem=getratematrix(initm,densitym)
        ratep=getratematrix(initp,densityp)
        if isconcise
            markov=[[initm,ratem],[initp,ratep]]
        else
            # col 8:14 = "Jiiij","Jiiji", "Jiijj","Jijii","Jijik","Jijjj","Jijkj"
            junctions=ancestry[i+1,8:14]
            densitymp = todensitymatrix(junctions)
            ratemp=getratematrix(initmp,densitymp)
            markov=[[initm,ratem],[initp,ratep],[initmp,ratemp]]
        end
        key=ancestry[i+1,1]
        sex=first(ped.gender[ped.member .== key])
        ismalex = (!isautosome) && lowercase(sex) == "male"
        key=>getmodelmarkov(markov,ismalex)
    end  for i in eachindex(memberlist))
end

function getratematrix(init::AbstractVector,density::AbstractMatrix)
    rate=Matrix(density)
    for i in eachindex(init)
        if init[i]>0
            rate[i,:] = rate[i,:] ./ init[i]
        else
            if !(sum(rate[i,:]) ≈ 0)
                error("inconsistency between ancestry coefficents and junction densities")
            end
        end
        rate[i,i] -= sum(rate[i,:])
    end
    if !(isapprox(sum(rate'*init), 0, atol=10^(-6.)))
        @warn string("non-stationary markov process", sum(rate'*init))
    end
    sparse(rate)
end

function todensitymatrix(junctions::AbstractVector)
    Jiiij, Jiiji, Jiijj, Jijii, Jijik, Jijjj, Jijkj = junctions
    n=size(Jiijj,1)
    # joint hidden state is a row-major order (ij)
    # for nfgl=4, states = 11, 12, 13, 14, 21, 22, 23, 24, 31, 32, 33, 34, 41, 42, 43, 44
    ii=[(i-1)n+i for i=1:n]
    # xi is a column-marjor order by julia convection
    xi=reshape(Vector(1:n^2),(n,n))
    ix=xi'
    density = zeros(n^2,n^2)
    for z=1:n
        density[ix[z,:],ix[z,:]] = Jijik[z,:,:]
        density[xi[z,:],xi[z,:]] = Jijkj[:,z,:]
        density[ii[z],ix[z,:]] = Jiiij[z,:]
        density[ii[z],xi[z,:]] = Jiiji[z,:]
        density[ix[z,:],ii[z]] = Jijii[z,:]
        # dedug from Jijjj[z,:] =>Jijjj[:,z], 2020/8/31*)
        density[xi[z,:],ii[z]] = Jijjj[:,z]
    end
    density[ii, ii] = Jiijj
    density
end

function todensitymatrix(n::Integer, j1122::Real,
    j1211::Real,j1213::Real,j1222::Real,j1232::Real)
    j1112=j1211
    j1121=j1222
    nondiag = ones(n,n)-I
    jiijj,jiiij,jiiji,jijii,jijjj=[j/(n*(n-1))*nondiag for j=[j1122,j1112,j1121,j1211,j1222]]
    jijik=jijkj=zeros(n,n,n)
    avg1213= n<=2 ? 0 : j1213/(n*(n-2)*(n-1))
    avg1232= n<=2 ? 0 : j1232/(n*(n-2)*(n-1))
    for i=1:n
        a=ones(n,n)-I
        a[i,:] = a[:,i] = zeros(n)
        jijik[i,:,:]=a*avg1213
        jijkj[i,:,:]=a*avg1232
    end
    todensitymatrix([jiiij, jiiji, jiijj, jijii, jijik, jijjj, jijkj])
end

function kronsum(a::AbstractMatrix,b::AbstractMatrix)
    kron(a,diagm(0=>ones(size(b,1)))) + kron(diagm(0=>ones(size(a,1))),b)
end

function getmodelmarkov(markov::Vector,ismalex::Bool)
    m, p = markov[1:2]
    if ismalex
         #male X
        dep = m
        n = length(m[1])
        ii = [(i-1)*n+i for i=1:n]
        init = spzeros(n^2)
        init[ii] = m[1]
        rate = spzeros(n^2,n^2)
        rate[ii,ii]=m[2]
        indep = [init,rate]
        joint = [init,rate]
    else
        # autosomes or female XX
        dep=(m+p) ./ 2
        indep=[[i*j for i=m[1] for j=p[1]],kronsum(m[2],p[2])]
        length(markov)>=3 && (joint = markov[3])
    end
    if length(markov)>=3
        (depmodel=dep,indepmodel=indep,jointmodel=joint,contmarkov=markov)
    else
        (depmodel=dep,indepmodel=indep,contmarkov=markov)
    end
end

function getinitrate(fgl::AbstractVector, mapR::Real)
    nfgl=length(fgl)
    n=sum(fgl)
    density= n==1 ? zeros(1,1) : (ones(n,n)-I)* (mapR/(n*(n-1)))
    init=ones(n)/n
    rate=getratematrix(init,density)
    # full range
    init2=spzeros(nfgl)
    rate2=spzeros(nfgl,nfgl)
    init2[fgl] = init
    rate2[fgl,fgl] =rate
    init2,rate2
end


# identityprior
function identityprior(ped::Pedigree,founderfgl::AbstractMatrix;
    memberlist::AbstractVector,
    isautosome::Bool,
    isconcise::Bool)
    # res[1,1:12]=["a", "phi12","R^m","R^p","rho","J1112",
    # "J1121","J1122","J1211","J1213","J1222","J1232"]
    res=pedidentity(ped,founderfgl,memberlist,isautosome,isconcise)    
    nfgl = size(pedstage.fglset,1)
    states=reshape(collect(1:nfgl^2),(nfgl,nfgl))'
    # memberfgl=[[trues(nfgl) for j=1:2] for i=1:length(pedstage.pairnolist)]
    memberfi=map(x->fi(x[1]), pedstage.pairnolist)
    memberfgl=[[.!iszero.(memberfi[i][j]) for j=1:2] for i=1:length(memberfi)]
    begin
        Dict(
        begin
            Rm,Rp=res[i+1,3:4]
            initm,ratem=getinitrate(memberfgl[i][1],Rm)
            initp,ratep=getinitrate(memberfgl[i][2],Rp)
            if isconcise
                markov=[[initm,ratem],[initp,ratep]]
            else
                iifgl=union(findall.(memberfgl[i])...)
                n=length(iifgl)
                nonibd=res[i+1,2]
                init=ones(n^2)*(nonibd/(n*(n-1)))
                for i=1:n
                    init[(i-1)n+i] =(1-nonibd)/n
                end
                # col 8:12 ="J1122","J1211","J1213","J1222","J1232"
                density = todensitymatrix(n, res[i+1,8:12]...)
                rate=getratematrix(init,density)
                mpfgl=states[iifgl,iifgl]'[:]
                initmp=spzeros(nfgl^2)
                ratemp=spzeros(nfgl^2,nfgl^2)
                initmp[mpfgl]=init
                ratemp[mpfgl,mpfgl]=rate
                markov=[[initm,ratem],[initp,ratep],[initmp,ratemp]]
            end
            key=res[i+1,1]
            sex=first(ped.gender[ped.member .== key])
            ismalex = (!isautosome) && lowercase(sex) == "male"
            key=>getmodelmarkov(markov,ismalex)
        end  for i in eachindex(memberlist))
    end
end

function identitymarkov(n::Integer, ismalex::Bool, phi12::Real, j1122::Real,
    j1211::Real,j1213::Real,j1222::Real,j1232::Real,isconcise::Bool)
    j1112=j1211
    j1121=j1222
    initm=initp=ones(n)/n
    initmp=ones(n^2)*(phi12/(n*(n-1)))
    for i=1:n
        initmp[(i-1)n+i] =(1-phi12)/n
    end
    avgRm = (j1122+j1121+j1222+j1232) / (n*(n-1))
    avgRp = (j1122+j1112+j1211+j1213) / (n*(n-1))
    densitym= (ones(n,n)-I)*avgRm
    densityp= (ones(n,n)-I)*avgRp
    ratem=getratematrix(initm,densitym)
    ratep=getratematrix(initp,densityp)
    if isconcise
        markov=[[initm,ratem],[initp,ratep]]
    else
        densitymp=todensitymatrix(n,j1122,j1211,j1213,j1222,j1232)
        ratemp=getratematrix(initmp,densitymp)
        markov=[[initm,ratem],[initp,ratep],[initmp,ratemp]]
    end
    markov
end

function identityprior(n::Integer, ismalex::Bool, phi12::Real, j1122::Real,
    j1211::Real,j1213::Real,j1222::Real,j1232::Real,isconcise::Bool)
    markov = identitymarkov(n,ismalex,phi12,j1122,j1211,j1213,j1222,j1232,isconcise)
    key = ismalex ? "malex" : "nonmalex"
    Dict([key=>getmodelmarkov(markov,ismalex)])
end

function identityprior(fglindicators::BitVector, ismalex::Bool, phi12::Real, j1122::Real,
    j1211::Real,j1213::Real,j1222::Real,j1232::Real,isconcise::Bool)
    n = sum(fglindicators)
    markov = identitymarkov(n,ismalex,phi12,j1122,j1211,j1213,j1222,j1232,isconcise)
    initm,ratem, initp,ratep = reduce(vcat,markov[1:2])
    n2 = length(fglindicators)
    initm2 = zeros(n2)
    initm2[fglindicators] .= initm
    initp2 = zeros(n2)
    initp2[fglindicators] .= initp
    ratem2 = zeros(n2,n2)
    ratem2[fglindicators,fglindicators] .= ratem
    ratep2 = zeros(n2,n2)
    ratep2[fglindicators,fglindicators] .= ratep
    if isconcise
        markov=[[initm2,ratem2],[initp2,ratep2]]
    else
        b = [i && j for i= fglindicators for j = fglindicators]
        initmp,ratemp = markov[3]        
        initmp2 = zeros(n2^2)
        initmp2[b] .= initmp
        ratemp2 = zeros(n2^2,n2^2)
        ratemp2[b,b] .= ratemp
        markov=[[initm2,ratem2],[initp2,ratep2],[initmp2,ratemp2]]
    end
    key = ismalex ? "malex" : "nonmalex"
    Dict([key=>getmodelmarkov(markov,ismalex)])
end

