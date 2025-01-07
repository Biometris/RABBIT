
function setstage(ped::Pedigree,isautosome::Bool, founderfgl::AbstractMatrix,
    indarray::AbstractMatrix)
    pedno = Pedigrees.toindexped(ped)    
    fglset = get_fglset(founderfgl)
    founderfgl2=Vector.(eachrow(founderfgl))
    dict=Dict(ped.member .=> pedno.member)    
    pairlist = Vector.(eachrow(indarray))
    pairnolist=map(x->[get(dict,i,0) for i in x],pairlist)
    global pedstage=(nfounder=ped.nfounder,memberno=pedno.member,
        genderno=pedno.gender,motherno=pedno.mother,fatherno=pedno.father,
        founderfgl=founderfgl2,fglset=fglset,
        isautosome=isautosome,pairlist=pairlist,pairnolist=pairnolist)
end

function get_fglset(founderfgl::AbstractMatrix)
    unique(founderfgl')
end

# phi12=two-gene non-ibd probability; phi12[a] =[mm,mp,pm,pp]
@memoize function phi12(a::Integer,b::Integer)
    if all([a,b] .<= pedstage.nfounder)        
        alab,blab=pedstage.founderfgl[[a,b]]
        return typeof(1.0)[i!=j for i = alab for j = blab]
    end
    if a>=b
        mo=pedstage.motherno
        fa=pedstage.fatherno
        mm = (a!=b)*mean(phi12(mo[a],b)[[1,3]])
        mp = mean(phi12(mo[a],b)[[2,4]])
        if pedstage.isautosome
            pm = mean(phi12(fa[a],b)[[1,3]])
            pp = (a!=b)*mean(phi12(fa[a],b)[[2,4]])
        else
            if pedstage.genderno[a]==1
                # 1=female, X chromosome
                pm = phi12(fa[a],b)[1]
                pp = (a!=b)*phi12(fa[a],b)[2]
            else
                pm = phi12(fa[a],b)[3]
                pp = (a!=b)*phi12(fa[a],b)[4]
            end
        end
        return [mm,mp,pm,pp]
    end
    phi12(b,a)[[1,3,2,4]]
end

# phi123=three-gene non-ibd probability
# phi123(a,b,c)=[mmm,mmp,mpm,mpp,pmm,pmp,ppm,ppp]
@memoize function phi123(a::Integer,b::Integer,c::Integer)
    if all([a,b,c] .<= pedstage.nfounder)        
        alab,blab,clab=pedstage.founderfgl[[a,b,c]]
        return typeof(1.0)[i != j != k != i for i = alab for j = blab for k = clab]
    end
    if a==b==c
        return zeros(8)
    end
    if a>=b>=c
        mo=pedstage.motherno
        fa=pedstage.fatherno
        mmm,mmp,mpm,mpp=[mean(phi123(mo[a],b,c)[[0,4] .+ i]) for i =1:4]
        if pedstage.isautosome
            pmm,pmp,ppm,ppp=[mean(phi123(fa[a],b,c)[[0,4] .+ i]) for i =1:4]
        else
            if pedstage.genderno[a]==1
                pmm,pmp,ppm,ppp=phi123(fa[a],b,c)[1:4]
            else
                pmm,pmp,ppm,ppp=phi123(fa[a],b,c)[5:end]
            end
        end
        mmm *= a!=b!=c
        mmp *= a!=b
        mpm *= a!=c
        mpp *= b!=c
        pmm *= b!=c
        pmp *= a!=c
        ppm *= a!=b
        ppp *= a!=b!=c
        return [mmm, mmp, mpm, mpp, pmm, pmp, ppm, ppp]
    end
    if a >=c >= b
        # [mmm,mmp,mpm,mpp,pmm,pmp,ppm,ppp]
        # [mmm,mpm,mmp,mpp,pmm,ppm,pmp,ppp]
        return phi123(a, c, b)[[1, 3, 2, 4, 5, 7, 6, 8]]
    end
    if b >= a >= c
        # [mmm,mmp,mpm,mpp,pmm,pmp,ppm,ppp]
        # [mmm,mmp,pmm,pmp,mpm,mpp,ppm,ppp]
        return phi123(b, a, c)[[1, 2, 5, 6, 3, 4, 7, 8]]
    end
    if b >= c >= a
        # [mmm,mmp,mpm,mpp,pmm,pmp,ppm,ppp]
        # [mmm,mpm,pmm,ppm,mmp,mpp,pmp,ppp]
        return phi123(b, c, a)[[1, 3, 5, 7, 2, 4, 6, 8]]
    end
    if c >= a >= b
        # [mmm,mmp,mpm,mpp,pmm,pmp,ppm,ppp]
        # [mmm,pmm,mmp,pmp,mpm,ppm,mpp,ppp]
        return phi123(c, a, b)[[1, 5, 2, 6, 3, 7, 4, 8]]
    end
    if c >= b >= a
        # [mmm,mmp,mpm,mpp,pmm,pmp,ppm,ppp]
        # [mmm,pmm,mpm,ppm,mmp,pmp,mpp,ppp]
        return phi123(c, b, a)[[1, 5, 3, 7, 2, 6, 4, 8]]
    end
    error("missing scenario phi123 for [a,b,c] =", [a, b, c]);
end

# R12(a)= [Rm,Rp] map expansion for maternally, patenrally derived chromosomes
@memoize function R12(a::Integer)
    if a<=pedstage.nfounder
        return zeros(2)
    end
    mo=pedstage.motherno
    fa=pedstage.fatherno
    Rm = mean(R12(mo[a])) + phi12(mo[a], mo[a])[2]
    if pedstage.isautosome
        Rp = mean(R12(fa[a])) + phi12(fa[a], fa[a])[2]
    else
        Rp = R12(fa[a])[pedstage.genderno[a]]
    end
    [Rm,Rp]
end

# junc1122=two-locus junction density; junc1122(a,b)=[mm,mp,pm,pp]
@memoize function junc1122(a::Integer,b::Integer)
    if all([a,b] .<= pedstage.nfounder)
        return zeros(4)
    end
    mo=pedstage.motherno
    fa=pedstage.fatherno
    if a==b
        mm, pp = R12(a)
        mp = pm = mean(junc1122(a, mo[a])[[3, 4]])
        return [mm,mp, pm,pp]
    end
    if a > b
        mm, mp = [mean(junc1122(mo[a], b)[[0, 2] .+ i]) for i = 1:2]
        if pedstage.isautosome
            pm, pp = [mean(junc1122(fa[a], b)[[0, 2] .+ i]) for i = 1:2]
        else
            if pedstage.genderno[a]==1
                pm, pp = junc1122(fa[a], b)[1:2]
            else
                pm, pp = junc1122(fa[a], b)[3:4]
            end
        end
        return [mm,mp, pm,pp]
    end
    # a < b
    junc1122(b, a)[[1, 3, 2, 4]]
end


@memoize function junc1232(a::Integer,b::Integer)
    if all([a,b] .<= pedstage.nfounder)
        return zeros(4)
    end
    mo=pedstage.motherno
    fa=pedstage.fatherno
    if a >= b
        mm = (a!=b)*(mean(junc1232(mo[a], b)[[1,3]]) + phi123(mo[a], mo[a], b)[3])
        mp = mean(junc1232(mo[a], b)[[2,4]]) + phi123(mo[a], mo[a], b)[4]
        if pedstage.isautosome
            pm = mean(junc1232(fa[a], b)[[1,3]]) + phi123(fa[a], fa[a], b)[3]
            pp = (a!=b)*(mean(junc1232(fa[a], b)[[2,4]]) + phi123(fa[a], fa[a], b)[4])
        else
            if pedstage.genderno[a]==1
                pm, pp = [1,  a!=b] .* junc1232(fa[a], b)[1:2]
            else
                pm, pp = [1,  a!=b] .* junc1232(fa[a], b)[3:4]
            end
        end
        return [mm,mp, pm,pp]
    else
        # a < b
        junc1213(b, a)[[1, 3, 2, 4]]
    end
end

@memoize function junc1213(a::Integer,b::Integer)
    if all([a,b] .<= pedstage.nfounder)
        return zeros(4)
    end
    mo=pedstage.motherno
    fa=pedstage.fatherno
    if a >= b
        mm = (a!=b)*mean(junc1213(mo[a], b)[[1,3]])
        mp = mean(junc1213(mo[a], b)[[2,4]])
        if pedstage.isautosome
            pm = mean(junc1213(fa[a], b)[[1,3]])
            pp = (a!=b)*mean(junc1213(fa[a], b)[[2,4]])
        else
            if pedstage.genderno[a]==1
                pm, pp = [1,  a!=b] .* junc1213(fa[a], b)[1:2]
            else
                pm, pp = [1,  a!=b] .* junc1213(fa[a], b)[3:4]
            end
        end
        return [mm,mp, pm,pp]
    else
        # a < b
        junc1232(b, a)[[1, 3, 2, 4]]
    end
end

# tocheck consistency between R and junc density
function pedidentity(ped::Pedigree,founderfgl::AbstractMatrix,
    indarray::AbstractMatrix, 
    isautosome::Bool,
    isconcise::Bool)
    setstage(ped,isautosome,founderfgl,indarray)
    pairnols=pedstage.pairnolist
    if isconcise
        res=Matrix{Any}(undef,size(pairnols,1)+1,5)
        res[1,:] = ["(a,b)", "phi12(a,b)","R(a)","R(b)","rho(a,b)"]
    else
        res=Matrix{Any}(undef,size(pairnols,1)+1,12)
        begin
            res[1,:] = ["(a,b)", "phi12(a,b)","R(a)","R(b)","rho(a,b)","J1112(a,b)",
            "J1121(a,b)","J1122(a,b)","J1211(a,b)","J1213(a,b)","J1222(a,b)","J1232(a,b)"]
        end
    end
    res[2:end,1] = pedstage.pairlist
    res[2:end,2]=map(x->phi12(x...), pairnols)
    res[2:end,3] =Ra=map(x->R12(x[1]), pairnols)
    res[2:end,4] =Rb=map(x->R12(x[2]), pairnols)
    j1122=map(x->junc1122(x...), pairnols)
    rho = map((x1,x2,x3)->[x1[1]+x2[1],x1[1]+x2[2],x1[2]+x2[1],x1[2]+x2[2]] .- x3,Ra,Rb,j1122)
    res[2:end,5] = rho
    if !isconcise
        res[2:end,8] =j1122
        res[2:end,10] =j1213=map(x->junc1213(x...), pairnols)
        res[2:end,12] =j1232=map(x->junc1232(x...), pairnols)
        ls=map(x->[x[1],x[1],x[2],x[2]],Ra)
        res[2:end,7] =  j1121 = res[2:end,11] = j1222 = (ls .- j1122 .- j1232) ./ 2
        ls=map(x->[x[1],x[2],x[1],x[2]],Rb)
        res[2:end,6]  = j1112 = res[2:end,9] = j1211 = (ls .- j1122 .- j1213) ./ 2
    end
    dememoize()
    res
end

function pedidentity(ped::Pedigree,founderfgl::AbstractMatrix,
    indarray::AbstractVector, 
    isautosome::Bool,
    isconcise::Bool)
    indarray2 = [i for i = indarray, j = 1:2]
    res = pedidentity(ped,founderfgl,indarray2,isautosome,isconcise)
    # res[1,:]=["(a,b)", "phi12(a,b)","R(a)","R(b)","rho(a,b)","J1112(a,b)",
    # "J1121(a,b)","J1122(a,b)","J1211(a,b)","J1213(a,b)","J1222(a,b)","J1232(a,b)"]
    res[1,:]=[replace(i,"(a,b)"=>"^mp") for i = res[1,:]]
    res[1,1]="a"
    res[1,3]=replace(res[1,3],"(a)"=>"^m")
    res[1,4]=replace(res[1,4],"(b)"=>"^p")
    # change  ["a", "phi12^mp","R^m","R^p","rho^mp","J1112^mp",
    # "J1121^mp","J1122^mp","J1211^mp","J1213^mp","J1222^mp","J1232^mp"]
    res[2:end,1] = getindex.(res[2:end,1],1)
    res[2:end,2] = getindex.(res[2:end,2],2)
    res[2:end,3] = getindex.(res[2:end,3],1)
    res[2:end,4] = getindex.(res[2:end,4],2)
    for i=5:size(res,2)
        res[2:end,i] = getindex.(res[2:end,i],2)
    end
    res
end
