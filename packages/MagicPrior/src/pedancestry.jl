
@memoize function fi(a::Integer)
    if a <= pedstage.nfounder        
        return [sparse(typeof(1.0).(pedstage.fglset .== i)) for i = pedstage.founderfgl[a]]
    end
    mo=pedstage.motherno
    fa=pedstage.fatherno
    m = mean(fi(mo[a]))
    p = pedstage.isautosome ? mean(fi(fa[a])) : fi(fa[a])[pedstage.genderno[a]]
    [m,p]
end

@memoize function phiij(a::Integer,b::Integer)
    if all([a,b] .<= pedstage.nfounder)
        # alab,blab=pedstage.founderfgl[[a,b]]
        return [sparse([k*l for k=fi(a)[i], l=fi(b)[j]]) for i=1:2 for j=1:2]
    end
    mo=pedstage.motherno
    fa=pedstage.fatherno
    if a==b
        mm,pp=[sparse(Diagonal(fi(a)[i])) for i=1:2]
        pm=mean(phiij(a,mo[b])[[3,4]])
        mp=permutedims(pm)
        return [mm, mp, pm, pp]
    end
    if a>b
        mm=mean(phiij(mo[a],b)[[1,3]])
        mp=mean(phiij(mo[a],b)[[2,4]])
        if pedstage.isautosome
            pm=mean(phiij(fa[a],b)[[1,3]])
            pp=mean(phiij(fa[a],b)[[2,4]])
        else
            if pedstage.genderno[a]==1
                pm=phiij(fa[a],b)[1]
                pp=phiij(fa[a],b)[2]
            else
                pm=phiij(fa[a],b)[3]
                pp=phiij(fa[a],b)[4]
            end
        end
        return [mm, mp, pm, pp]
    else
        # a<b
        map(x->permutedims(x), phiij(b, a)[[1, 3, 2, 4]])
    end
end

# permutedims(a,[i,j,k]) is same as Transpose[a,{i,j,k}] of MMA
# except for [i,j,k] = [2,3,1], [3,1,2]
# that is,  Transpose[a,{2,3,1}] = permutedims(a,[3,1,2])
# that is,  Transpose[a,{3,1,2}] = permutedims(a,[2,3,1])
# a=reshape(Vector(1:24), (2,3,4))
# b=permutedims(a,[2,3,1])
# Permutations[Range[3]]
# a = ConstantArray[0, {2, 3, 4}]
# Do[a[[All, All, i]] =
#   Transpose[Partition[Range[(i - 1)*6 + 1, 6* i], 2]], {i, 4}]
# Table[a[[All, All, i]] // MatrixForm, {i, Length[a[[1, 1]]]}]
# b = Transpose[a, {3, 1, 2}];
# Dimensions[b]
# Table[b[[All, All, i]] // MatrixForm, {i, Length[b[[1, 1]]]}]
@memoize function phiijk(a::Integer,b::Integer,c::Integer)
    if all([a,b,c] .<= pedstage.nfounder)
        # alab,blab,clab=pedstage.founderfgl[[a,b,c]]
        # sparse() only for 1D and 2D, but not 3D array yet;
        # to check package e.g. SimpleSparseArrays
        return  [[k*l*m for k=fi(a)[i1], l=fi(b)[i2], m=fi(c)[i3]] for i1=1:2 for i2=1:2 for i3=1:2]
    end
    mo=pedstage.motherno
    fa=pedstage.fatherno
    n=length(pedstage.fglset)
    if a==b==c
        mmm,ppp=[[(i==j==k)*fi(a)[ch][i] for i=1:n, j=1:n, k=1:n] for ch=1:2]
        mmp,ppm=[[(i==j)*phiij(a,a)[ch][i,k] for i=1:n, j=1:n, k=1:n] for ch=2:3]
        mpm=permutedims(mmp,[1,3,2])
        pmp=permutedims(ppm,[1,3,2])
        mpp=permutedims(ppm,[3,1,2])
        pmm=permutedims(mmp,[3,1,2])
        return [mmm, mmp, mpm, mpp, pmm, pmp, ppm, ppp]
    end
    if a==b>c
        mmm,mmp,ppm,ppp=[[(i==j)*phiij(a,c)[ch][i,k] for i=1:n, j=1:n, k=1:n] for ch=1:4]
        # pmm=5,ppm=7
        pmm=mean(phiijk(a,mo[a],c)[[5,7]])
        mpm=permutedims(pmm,[2,1,3])
        # pmp=6,ppp=8
        pmp=mean(phiijk(a,mo[a],c)[[6, 8]]);
        mpp=permutedims(pmp, [2,1,3]);
        return [mmm, mmp, mpm, mpp, pmm, pmp, ppm, ppp]
    end
    if a>b==c
        mmm,mpp,pmm,ppp=[[(j==k)*(phiij(a,b)[ch][i,j]) for i=1:n,j=1:n, k=1:n] for ch=1:4]
        mmp,mpm=[mean(phiijk(mo[a],b,c)[[0,4] .+ i]) for i=2:3]
        if pedstage.isautosome
            pmp,ppm=[mean(phiijk(fa[a],b,c)[[0,4] .+ i]) for i=2:3]
        else
            if pedstage.genderno[a]==1
                pmp,ppm=phiijk(fa[a],b,c)[[2,3]]
            else
                pmp,ppm=phiijk(fa[a],b,c)[[6,7]]
            end
        end
        return [mmm, mmp, mpm, mpp, pmm, pmp, ppm, ppp]
    end
    if a>b>c
        mmm,mmp,mpm,mpp=[mean(phiijk(mo[a],b,c)[[0,4] .+ i]) for i=1:4]
        if pedstage.isautosome
            pmm,pmp,ppm,ppp=[mean(phiijk(fa[a],b,c)[[0,4] .+ i]) for i=1:4]
        else
            if pedstage.genderno[a]==1
                pmm,pmp,ppm,ppp=phiijk(fa[a],b,c)[1:4]
            else
                pmm,pmp,ppm,ppp=phiijk(fa[a],b,c)[5:8]
            end
        end
        return [mmm, mmp, mpm, mpp, pmm, pmp, ppm, ppp]
    end
    if a >=c >= b
        # [mmm,mmp,mpm,mpp,pmm,pmp,ppm,ppp]
        # [mmm,mpm,mmp,mpp,pmm,ppm,pmp,ppp]
        return map(x->permutedims(x,[1,3,2]),phiijk(a, c, b)[[1, 3, 2, 4, 5, 7, 6, 8]])
    end
    if b >= a >= c
        # [mmm,mmp,mpm,mpp,pmm,pmp,ppm,ppp]
        # [mmm,mmp,pmm,pmp,mpm,mpp,ppm,ppp]
        return map(x->permutedims(x,[2,1,3]),phiijk(b, a, c)[[1, 2, 5, 6, 3, 4, 7, 8]])
    end
    if b >= c >= a
        # [mmm,mmp,mpm,mpp,pmm,pmp,ppm,ppp]
        # [mmm,mpm,pmm,ppm,mmp,mpp,pmp,ppp]
        return map(x->permutedims(x,[3,1,2]),phiijk(b, c, a)[[1, 3, 5, 7, 2, 4, 6, 8]])
    end
    if c >= a >= b
        # [mmm,mmp,mpm,mpp,pmm,pmp,ppm,ppp]
        # [mmm,pmm,mmp,pmp,mpm,ppm,mpp,ppp]
        return map(x->permutedims(x,[2,3,1]),phiijk(c, a, b)[[1, 5, 2, 6, 3, 7, 4, 8]])
    end
    if c >= b >= a
        # [mmm,mmp,mpm,mpp,pmm,pmp,ppm,ppp]
        # [mmm,pmm,mpm,ppm,mmp,pmp,mpp,ppp]
        return map(x->permutedims(x,[3,2,1]),phiijk(c, b, a)[[1, 5, 3, 7, 2, 6, 4, 8]])
    end
    error("missing scenario phiijk for [a,b,c] =", [a, b, c]);
end

@memoize function Rij(a::Integer)
    n=length(pedstage.fglset)
    if a<=pedstage.nfounder
        return [sparse(zeros(n,n)) for ch=1:2]
    end
    mo=pedstage.motherno
    fa=pedstage.fatherno
    nondiag = 1 .- Diagonal(ones(n))
    Rm=mean(Rij(mo[a]))+mean(phiij(mo[a], mo[a])[[2, 3]]) .* nondiag
    if pedstage.isautosome
        Rp=mean(Rij(fa[a]))+mean(phiij(fa[a], fa[a])[[2, 3]]) .* nondiag
    else
        Rp = Rij(fa[a])[pedstage.genderno[a]]
    end
    [sparse(Rm),sparse(Rp)]
end


@memoize function junciijj(a::Integer,b::Integer)
    n=length(pedstage.fglset)
    if all([a,b] .<= pedstage.nfounder)
        return [sparse(zeros(n,n)) for ch=1:4]
    end
    mo=pedstage.motherno
    fa=pedstage.fatherno
    if a==b
        mm, pp = Rij(a)
        mp = pm = mean(junciijj(a, mo[a])[[3, 4]])
        return [mm, mp, pm, pp]
    elseif a>b
        mm, mp = [mean(junciijj(mo[a],b)[[0,2] .+ i]) for i=1:2]
        if pedstage.isautosome
            pm, pp = [mean(junciijj(fa[a],b)[[0,2] .+ i]) for i=1:2]
        else
            if pedstage.genderno[a]==1
                pm,pp=junciijj(fa[a],b)[1:2]
            else
                pm,pp=junciijj(fa[a],b)[3:4]
            end
        end
        return [mm, mp, pm, pp]
    else
        # a<b
        junciijj(b, a)[[1, 3, 2, 4]]
    end
end

@memoize function juncijjj(a::Integer,b::Integer)
    n=length(pedstage.fglset)
    if all([a,b] .<= pedstage.nfounder)
        return [sparse(zeros(n,n)) for ch=1:4]
    end
    mo=pedstage.motherno
    fa=pedstage.fatherno
    nondiag = 1 .- Diagonal(ones(n))
    if a>=b
        ls=[sparse(hcat([diag(p[i,:,:]) for i=1:n]...)') for p = phiijk(mo[a],mo[a],b)]
        ls2=[nondiag .* mean(ls[[2,4] .+ i]) for i=1:2]
        mm,mp=[mean(juncijjj(mo[a],b)[[0,2] .+ i]) .+ ls2[i] for i=1:2]
        if pedstage.isautosome
            ls=[sparse(hcat([diag(p[i,:,:]) for i=1:n]...)') for p = phiijk(fa[a],fa[a],b)]
            ls2=[nondiag .* mean(ls[[2,4] .+ i]) for i=1:2]
            pm,pp=[mean(juncijjj(fa[a],b)[[0,2] .+ i]) .+ ls2[i] for i=1:2]
        else
            if pedstage.genderno[a]==1
                pm,pp=juncijjj(fa[a],b)[1:2]
            else
                pm,pp=juncijjj(fa[a],b)[3:4]
            end
        end
        if a==b
            mm *= 0
            pp *= 0
        end
        return [mm, mp, pm, pp]
    else
        # a<b
        map(x->permutedims(x),juncijii(b, a)[[1, 3, 2, 4]])
    end
end

@memoize function juncijii(a::Integer,b::Integer)    
    n=length(pedstage.fglset)
    if all([a,b] .<= pedstage.nfounder)
        return [sparse(zeros(n,n)) for ch=1:4]
    end
    mo=pedstage.motherno
    fa=pedstage.fatherno
    if a>=b
        mm, mp = [mean(juncijii(mo[a],b)[[0,2] .+ i]) for i=1:2]
        if pedstage.isautosome
            pm, pp = [mean(juncijii(fa[a],b)[[0,2] .+ i]) for i=1:2]
        else
            if pedstage.genderno[a]==1
                pm,pp=juncijii(fa[a],b)[1:2]
            else
                pm,pp=juncijii(fa[a],b)[3:4]
            end
        end
        if a==b
            mm *= 0
            pp *= 0
        end
        return [mm, mp, pm, pp]
    else
        # a<b
        map(x->permutedims(x),juncijjj(b, a)[[1, 3, 2, 4]])
    end
end

@memoize function juncijkj(a::Integer,b::Integer)
    n=length(pedstage.fglset)
    if all([a,b] .<= pedstage.nfounder)
        return [zeros(n,n,n) for ch=1:4]
    end
    mo=pedstage.motherno
    fa=pedstage.fatherno
    nondiag = [length(union([i,j,k]))==3 ? 1.0 : 0.0 for i =1:n, j=1:n, k=1:n]
    if a>=b
        ls=[[p[i,k,j] for i =1:n, j=1:n, k=1:n] for p = phiijk(mo[a],mo[a],b)]
        ls2=[nondiag .* mean(ls[[2,4] .+ i]) for i =1:2]
        mm,mp=[mean(juncijkj(mo[a],b)[[0,2] .+ i]) .+ ls2[i] for i=1:2]
        if pedstage.isautosome
            ls=[[p[i,k,j] for i =1:n, j=1:n, k=1:n] for p = phiijk(fa[a],fa[a],b)]
            ls2=[nondiag .* mean(ls[[2,4] .+ i]) for i =1:2]
            pm,pp=[mean(juncijkj(fa[a],b)[[0,2] .+ i]) .+ ls2[i] for i=1:2]
        else
            if pedstage.genderno[a]==1
                pm,pp=juncijkj(fa[a],b)[1:2]
            else
                pm,pp=juncijkj(fa[a],b)[3:4]
            end
        end
        if a==b
            mm *= 0
            pp *= 0
        end
        return [mm, mp, pm, pp]
    else
        # a<b
        map(x->permutedims(x,[2,1,3]),juncijik(b, a)[[1, 3, 2, 4]])
    end
end

@memoize function juncijik(a::Integer,b::Integer)
    n=length(pedstage.fglset)
    if all([a,b] .<= pedstage.nfounder)
        return [zeros(n,n,n) for ch=1:4]
    end
    mo=pedstage.motherno
    fa=pedstage.fatherno
    if a>=b
        mm,mp=[mean(juncijik(mo[a],b)[[0,2] .+ i]) for i=1:2]
        if pedstage.isautosome
            pm,pp=[mean(juncijik(fa[a],b)[[0,2] .+ i]) for i=1:2]
        else
            if pedstage.genderno[a]==1
                pm,pp=juncijik(fa[a],b)[1:2]
            else
                pm,pp=juncijik(fa[a],b)[3:4]
            end
        end
        if a==b
            mm *= 0
            pp *= 0
        end
        return [mm, mp, pm, pp]
    else
        # a<b
        map(x->permutedims(x,[2,1,3]),juncijkj(b, a)[[1, 3, 2, 4]])
    end
end

function ancestryConsistencyQ(resancestry::AbstractMatrix)
    # resancestry[1,1:14] = ["(a)","fi^m","fi^p", "phiij", "Rij^m","Rij^p","rhoij",
    #       "Jiiij","Jiiji", "Jiijj","Jijii","Jijik","Jijjj","Jijkj"]
    restest=Vector{Bool}()
    for i=1:size(resancestry,1)-1
        Ra,Rb,rho,Jiiij,Jiiji, Jiijj,Jijii,Jijik,Jijjj,Jijkj= resancestry[i+1,5:14]
        Jikjk=[permutedims(p,[1,3,2]) for p = Jijkj]
        Jisjs=[sum(p,dims=3)[:,:,1] for p=Jikjk]
        a=Jiijj .+Jiiji .+ Jijjj .+ Jisjs
        testa=isapprox(a, [Ra[1],Ra[1],Ra[2],Ra[2]], atol=10^(-6.))
        Jkikj=[permutedims(p,[2,3,1]) for p = Jijik]
        Jsisj=[sum(p,dims=3)[:,:,1] for p=Jkikj]
        Jjijj=map(permutedims,Jiiij)
        b=Jiijj .+ Jiiij .+ Jjijj .+ Jsisj
        testb=isapprox(b, [Rb[1],Rb[2],Rb[1],Rb[2]], atol=10^(-6.))
        push!(restest,testa && testb)
    end
    all(restest)
end

function pedancestry(ped::Pedigree,founderfgl::AbstractMatrix,
    indarray::AbstractMatrix, 
    isautosome::Bool,
    isconcise::Bool)
    setstage(ped,isautosome,founderfgl,indarray)
    pairnols=pedstage.pairnolist
    if isconcise
        res=Matrix{Any}(undef,size(pairnols,1)+1,7)
        res[1,:] = ["(a,b)", "fi(a)","fi(b)","phi12(a,b)","R(a)","R(b)","rho(a,b)"]
    else
        res=Matrix{Any}(undef,size(pairnols,1)+1,14)        
        res[1,:] = ["(a,b)","fi(a)","fi(b)", "phiij(a,b)", "Rij(a)","Rij(b)","rhoij(a,b)",
            "Jiiij(a,b)","Jiiji(a,b)", "Jiijj(a,b)","Jijii(a,b)","Jijik(a,b)","Jijjj(a,b)","Jijkj(a,b)"]    
    end
    res[2:end,1] = pedstage.pairlist
    res[2:end,2]=map(x->fi(x[1]), pairnols)
    res[2:end,3]=map(x->fi(x[2]), pairnols)
    res[2:end,4]=map(x->phiij(x...), pairnols)
    res[2:end,5]=Ra=map(x->Rij(x[1]), pairnols)
    res[2:end,6]=Rb=map(x->Rij(x[2]), pairnols)
    Jiijj=map(x->junciijj(x...), pairnols)
    begin
        res[2:end,7] =[[Ra[p][ch1] .+ Rb[p][ch2] .- Jiijj[p][2(ch1-1)+ch2]
        for ch1=1:2 for ch2=1:2] for p=1:size(pairnols,1)]
    end
    if !isconcise
        res[2:end,10]=Jiijj
        res[2:end,11]=Jijii=map(x->juncijii(x...), pairnols)
        res[2:end,12]=Jijik=map(x->juncijik(x...), pairnols)
        res[2:end,13]=Jijjj=map(x->juncijjj(x...), pairnols)
        res[2:end,14]=Jijkj=map(x->juncijkj(x...), pairnols)
        res[2:end,8]=Jiiij=Jijii
        res[2:end,9]=Jiiji=[map(permutedims,p) for p=Jijjj]
        if !ancestryConsistencyQ(res)
            @warn "inconsistency among junction densities"
        end
    end
    dememoize()
    res
end

function pedancestry(ped::Pedigree,founderfgl::AbstractMatrix,
    indarray::AbstractVector, 
    isautosome::Bool,
    isconcise::Bool)
    indarray2 = [i for i = indarray, j = 1:2]
    res = pedancestry(ped,founderfgl,indarray2,isautosome,isconcise)
    # res[1,:] = ["(a,b)","fi(a)","fi(b)", "phiij(a,b)", "Rij(a)","Rij(b)","rhoij(a,b)",
    # "Jiiij(a,b)","Jiiji(a,b)", "Jiijj(a,b)","Jijii(a,b)","Jijik(a,b)","Jijjj(a,b)","Jijkj(a,b)"]
    res[1,:]=[replace(i,"(a,b)"=>"^mp") for i = res[1,:]]
    res[1,1]="a"
    res[2:end,1] = getindex.(res[2:end,1],1)
    for j=[2,5]
        res[1,j]=replace(res[1,j],"(a)"=>"^m")
        res[2:end,j] = getindex.(res[2:end,j],1)
    end
    for j=[3,6]
        res[1,j]=replace(res[1,j],"(b)"=>"^p")
        res[2:end,j] = getindex.(res[2:end,j],2)
    end
    res[2:end,4] = getindex.(res[2:end,4],2)
    for j=7:size(res,2)
        res[2:end,j] = getindex.(res[2:end,j],2)
    end
    res
end

function isexchcompatible(residentity::AbstractMatrix,resancestry::AbstractMatrix)
    # residentity[1,1:12]=["(a)", "phi12","R^m","R^p","rho","J1112",
    # "J1121","J1122","J1211","J1213","J1222","J1232"]
    # resancestry[1,1:14] = ["(a)","fi^m","fi^p", "phiij", "Rij^m","Rij^p","rhoij",
    # "Jiiij","Jiiji", "Jiijj","Jijii","Jijik","Jijjj","Jijkj"]
    res = Matrix{Any}(undef,size(resancestry)...)
    if (ndims(residentity[2,2]),ndims(resancestry[2,4]))==(0,2)
        for i=1:size(res,1)-1
            res[i+1,4]=1-sum(diag(resancestry[i+1,4]))
            res[i+1,5:end] = sum.(resancestry[i+1,5:end])
        end
    elseif (ndims(residentity[2,2]),ndims(resancestry[2,4]))==(1,1)
        for i=1:size(res,1)-1
            res[i+1,4]=1 .- sum.(diag.(resancestry[i+1,4]))
            res[i+1,5:end] = map(x->sum.(x), resancestry[i+1,5:end])
        end
    else
        @error "residentity and resancestry do not match"
    end
    jmax=min(size(residentity,2)+2,size(res,2))
    isapprox(res[2:end,4:jmax],residentity[2:end,2:jmax-2],atol=10^(-6.))
end

"""
    isexchcompatible(pedigree::Pedigree;isautosome::Bool=true)

check whether expected ancestral probabilities and junction densities assuming FGLs
being non-exchangeable are compatible with those assuming FGLs being exchangeable.

"""
function isexchcompatible(pedigree::Pedigree;isautosome::Bool=true)
    nfounder=pedigree.nfounder
    # founderfgl= [2*(i-1)+j for i = 1:nfounder, j = 1:2]
    founderfgl= [i for i = 1:nfounder, j = 1:2]
    memberlist=pedigree.member[nfounder+1:end]
    isconcise=false
    residentity=pedidentity(pedigree,founderfgl,memberlist,isautosome,isconcise)
    resancestry=pedancestry(pedigree,founderfgl,memberlist,isautosome,isconcise)    
    dememoize()
    isexchcompatible(residentity,resancestry)
end
