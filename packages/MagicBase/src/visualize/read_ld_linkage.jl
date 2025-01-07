

function read_ld(ldfile::AbstractString;
    delim::AbstractChar=',',
    commentstring::AbstractString = "##", 
    workdir::AbstractString=pwd())
    ldfile2 = getabsfile(workdir,ldfile)            
    res= MagicBase.readmultitable(ldfile2;
        delim,
        missingstring=["missing"], # keep NA 
        commentstring,        
    )
    individuals = res["offspringinfo"][!,:individual]
    markerinfo = res["markerinfo"]
    markers = string.(markerinfo[!,:marker_represent])
    markerbinls = convert.(String, string.(markerinfo[!,:marker_bin]))
    physchromls = convert.(String, string.(markerinfo[!,:physchrom_bin]))
    physposbpls = convert.(String, string.(markerinfo[!,:physposbp_bin]))
    if any(markerbinls .== "NA") 
        dupebindict =  nothing 
    else
        dupebindict = Dict(markers .=> tuple.(markerbinls, physchromls,physposbpls))        
    end
    df = res["pairwiseld"]
    nsnp = length(markers)    
    # only snp1<snp2 were calculated and saved in the ldfile    
    # LD r2,squred allelic correlation
    iils, jjls = df[!,1], df[!,2]
    iils, jjls = vcat(iils, jjls, 1:nsnp), vcat(jjls, iils, 1:nsnp)
    vvls = df[!,3]
    vvls = vcat(vvls,vvls,ones(nsnp))
    ldr2  = sparse(iils,jjls, vvls, nsnp,nsnp)
    dropzeros!(ldr2)
    # LD lod sores for independence measurement 
    vvls = df[!,4]
    vvls = vcat(vvls,vvls,ones(nsnp))
    lod = DataFrame(i=iils, j=jjls, v=vvls)
    # calculate diagonal as the row max
    gd = groupby(lod,[:i])
    v1max = [maximum(g[!,:v]) for g in gd]
    gd = groupby(lod,[:j])
    v2max = [maximum(g[!,:v]) for g in gd]
    vdiag = map(max,v1max,v2max)
    vlimit = 10^(ndigits(ceil(Int,maximum(vdiag)))+1)
    vdiag[vdiag .≈ 1.0] .= vlimit
    vvls[2*size(df,1)+1:end] .= vdiag
    ldlod  = sparse(iils,jjls, vvls, nsnp,nsnp)
    (individuals=individuals,markers = markers, dupebindict = dupebindict, ldr2 = ldr2, ldlod=ldlod)
end

function read_linkage(linkagefile::AbstractString;
    delim::AbstractChar=',',
    commentstring::AbstractString = "##", 
    workdir::AbstractString=pwd())
    linkagefile2 = getabsfile(workdir,linkagefile)        
    res= MagicBase.readmultitable(linkagefile2;
        delim,
        missingstring=["missing"], # keep NA 
        commentstring
    )
    individuals = res["offspringinfo"][!,:individual]
    markerinfo = res["markerinfo"]    
    markers = convert.(String, string.(markerinfo[!,:marker]) )
    physchromls = convert.(String, string.(markerinfo[!,:physchrom]))
    physposbpls = convert.(String, string.(markerinfo[!,:physposbp]))
    physmapdict = Dict(markers .=> tuple.(physchromls,physposbpls))
    nmissingls = markerinfo[!,:nmissing]
    df = res["pairwiselinkage"]
    # only snp1<snp2 were calculated and saved in the linkagefile    
    nsnp = length(markers)    
    # recombination fraction
    iils, jjls = df[!,1], df[!,2]
    iils, jjls = vcat(iils, jjls, 1:nsnp), vcat(jjls, iils, 1:nsnp)
    vvls = 1.0 .- df[!,3]
    b = vvls .< 0.0
    if any(b) 
        @warn "unexpected negagive recom_fraction"
        vvls[b] .= 0.0
    end
    vvls = vcat(vvls,vvls,ones(nsnp))
    recomnonfrac  = sparse(iils,jjls, vvls, nsnp,nsnp)
    dropzeros!(recomnonfrac)
    # recombination LOD score
    vvls = df[!,4]
    vvls = vcat(vvls,vvls,ones(nsnp))
    lod = DataFrame(i=iils, j=jjls, v=vvls)
    # calculate diagonal as the row max
    gd = groupby(lod,[:i])
    v1max = [maximum(g[!,:v]) for g in gd]
    gd = groupby(lod,[:j])
    v2max = [maximum(g[!,:v]) for g in gd]
    vdiag = map(max,v1max,v2max)
    vlimit = 10^(ndigits(ceil(Int,maximum(vdiag)))+1)
    vdiag[vdiag .≈ 1.0] .= vlimit
    vvls[2*size(df,1)+1:end] .= vdiag
    recomlod  = sparse(iils,jjls, vvls, nsnp,nsnp)
    dropzeros!(recomlod)
    (individuals=individuals,markers = markers, physmapdict=physmapdict, 
        nmissingls = nmissingls, recomnonfrac=recomnonfrac, recomlod=recomlod)
end
