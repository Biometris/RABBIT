

"""
    rabbitped_mma2jl(mmapedfile; kwargs...)

convert input pedfile of Mathematica-version RABBIT into pedfile for Julia-version RABBIT.

# keyword arguments

`ishomozygous::Bool = false`: if true, offspring is homozygous. 

`isfglexch::Bool = false`: if true, offspring is produced with rand parent ordering.

`workdir::AbstractString=pwd()`: working directory for input and output files.

`outfile::AbstractString = "outstem_ped.csv"`: output filename.

"""
function rabbitped_mma2jl(mmapedfile::AbstractString;    
    ishomozygous::Bool = false, 
    isfglexch::Bool = false, 
    outfile::AbstractString = "outstem_ped.csv",
    workdir::AbstractString=pwd())
    pedfile = getabsfile(workdir,mmapedfile)
    isfile(pedfile) || @error string(pedfile, " does not exist")
    pedinfo = MagicBase.readmultitable(pedfile;delim=',',commentstring="##",isordereddict=true);
    infols = collect(values(pedinfo))
    length(infols ) == 2 || @error string("wrong number ",length(infols), ",  of talbes in pedinfo")
    size(infols[1],2) >=5 || @error string("wrong number ",size(infols[1],2), ", of columns in designinfo of  pedinfo")
    designdf = infols[1][!,[2,4,5]]
    if names(designdf) != ["MemberID","MotherID","FatherID"] 
        @warn string("designinfo colnames  ",names(designdf), " are not the expected ", ["MemberID","MotherID","FatherID"])
    end
    rename!(designdf,["member","mother","father"])
    size(infols[2],2) >=3 || @error string("wrong number ",size(infols[2],2), ", of columns in offspring of  pedinfo")
    offspringinfo = infols[2][!,1:3]
    if names(offspringinfo) != ["ProgenyID","MemberID","Funnelcode"]
        @warn string("offspringinfo colnames ",names(offspringinfo), " are not the expected colnames ", ["ProgenyID","MemberID","Funnelcode"])    
    end
    fcodes = unique(offspringinfo[:,3])
    if length(fcodes) > 1
        @warn string("account for different funnel codes by isfglexch=true")
        if !isfglexch
            @info string("reset isfglexch=true")
            isfglexch = true
        end
    end
    rename!(offspringinfo,["individual","member","ishomozygous"])
    offspringinfo[!,3] .= ishomozygous     
    insertcols!(offspringinfo,4,:isfglexch=>isfglexch)    
    designinfo = df2designinfo(designdf)
    founderinfo = get_founderinfo(designinfo)
    magicped = MagicPed(designinfo, founderinfo, offspringinfo)
    savemagicped(outfile, magicped;workdir)    
end

"""
    rabbitgeno_mma2jl(mmagenofile; kwargs...)

convert input genofile of Mathematica-version RABBIT into genofile for Julia-version RABBIT. 

# keyword arguments

`workdir::AbstractString=pwd()`: working directory for input and output files.

`outfile::AbstractString = "outstem_geno.vcf.gz"`: output filename.

"""
function rabbitgeno_mma2jl(mmagenofile::AbstractString;
    isfounderphased::Bool=false,
    outfile::AbstractString = "outstem_geno.vcf.gz",
    workdir::AbstractString=pwd())
    ext = last(splitext(mmagenofile))
    ext == ".csv" || @error string("file ext = ",ext, ", not .csv")
    inputfile = getabsfile(workdir, mmagenofile)
    genodata = readdlm(inputfile,',');
    nf = genodata[1,2]
    @info string("#founders=",nf)
    genodf =DataFrame(permutedims(string.(genodata[5:end,2:end])), string.(genodata[5:end,1]))
    mapdf = DataFrame(permutedims(genodata[2:4,2:end]),["marker","linkagegroup","poscm"])
    mapdf[!,1] .= string.(mapdf[!,1])
    mapdf[!,2] .= [i == "NA" ? missing : string(i) for i in mapdf[!,2]]
    mapdf[!,3] .= [i == "NA" ? missing : Float32(i) for i in mapdf[!,3]]
    gset = unique(Matrix(genodf[!,1:nf]))
    if issubset(gset,["1","2","N"])
        fformat = "GT_haplo"    
    elseif issubset(gset,["NN", "N1", "1N", "N2", "2N", "11", "12", "21", "22"])
        if isfounderphased
            fformat = "GT_phased"
            for j in 1:nf
                genodf[!,j] .= [join(split(i,""),"|") for i in genodf[!,j]]
            end
        else
            fformat = "GT_unphased"        
        end        
    elseif all([sum(i .== "|")==1 for i in split.(gset,"")])
        fformat = "AD"
    else
        @error string("unknown founder genotypes: ",gset)
    end
    gset = unique(Matrix(genodf[!,nf+1:end]))
    if issubset(gset,["NN", "N1", "1N", "N2", "2N", "11", "12", "21", "22"])
        offformat = "GT_unphased"
    elseif all([sum(i .== "|")==1 for i in split.(gset,"")])
        offformat = "AD"
    elseif all([in(sum(i .== "|"),[3,4]) for i in split.(gset,"")])
        offformat = "GP"
    else
        @error string("unknown offspring genotypes: ",gset)
    end    
    middf = DataFrame(physchrom = [missing for _ in 1:size(mapdf,1)], 
        physposbp=missing, info=missing,
        founderformat=fformat,offspringformat=offformat,
        foundererror=missing,offspringerror=missing,baseerror=missing,
        allelicbias=missing,allelicoverdispersion=missing,allelicdropout=missing)
    resdf = hcat(mapdf,middf,genodf)
    isfounderinbred = fformat == "GT_haplo"
    indls = names(genodf)
    noff = length(indls) - nf
    magicped= formmagicped(JuncDist(nfounder=nf),noff)
    MagicBase.setfounderid!(magicped,indls[1:nf])
    MagicBase.setoffspringid!(magicped,indls[nf+1:end])
    magicgeno = MagicBase.formmagicgeno!(resdf, magicped;isfounderinbred, isphysmap=false)
    outfile2 = getabsfile(workdir, outfile)
    savegenodata(outfile2,magicgeno; delim='\t',commentstring = "##",workdir)
end



function magicgeno_2rabbitmma(genofile::AbstractString,pedinfo::AbstractString;
    isfounderinbred::Bool=true,
    formatpriority::AbstractVector=["AD","GT"],    
    isphysmap::Bool=false,
    recomrate::Real=1.0,
    isdelmultiallelic::Bool=true,
    commentstring::AbstractString="##",
    missingstring=["NA","missing"],    
    outstem::AbstractString = "outstem",
    workdir::AbstractString=pwd())
    magicgeno = formmagicgeno(genofile, pedinfo; isphysmap,recomrate,formatpriority,isfounderinbred, isdelmultiallelic, commentstring, missingstring, workdir); 
    magicgeno_2rabbitmma(magicgeno; outstem, workdir)
end


function magicgeno_2rabbitmma(magicgeno::MagicGeno; 
    outstem::AbstractString = "outstem",
    workdir::AbstractString=pwd())    
    fls = magicgeno.magicped.founderinfo[!,:individual]
    offfls = magicgeno.magicped.offspringinfo[!,:individual]
    nfounder = length(fls)
    samplels = vcat(fls,offfls)
    outfile = joinpath(workdir, string(outstem,"_mma_geno.csv"))
    open(outfile,"w") do io    
        delim2 = ','
        write(io, string("nfounder,",nfounder),"\n")
        # save markermap
        markermap = reduce(vcat, magicgeno.markermap)
        markermap2 = markermap[!,1:3]
        mapmtx = hcat(reshape(names(markermap2),:,1),permutedims(Matrix(markermap2)))
        writedlm(io,mapmtx,delim2)    
        # save genodata    
        colls = _col_1stsample:(_col_1stsample + length(samplels) -1)
        resmtx = reduce(hcat,[begin 
            mtx = togenomtx(magicgeno, chr; isvcf=false, missingstring = "NA", target="all")            
            mtx2 = permutedims(mtx[:,colls])
            for i in eachindex(mtx2)
                g = mtx2[i]
                if occursin("&",g)
                    mtx2[i] = replace(g, "&"=>"|")
                elseif occursin("|",g)
                    ls = split(g, "|")
                    mtx2[i] = join(replace(x-> x== "0" ? "1" : "2", ls))
                elseif occursin("/",g)                    
                    ls = split(g, "/")
                    mtx2[i] = join(replace(x-> x== "0" ? "1" : "2", ls))
                elseif in(g, ["11","12","21","22","N1","1N","N2","2N","NN","1","2","N"])
                    0 # do nothing
                else
                    @error string("unknown genotype = ", g) maxlog=10                    
                end
            end
            mtx2
        end for chr in eachindex(magicgeno.markermap)])    
        resmtx = hcat(reshape(samplels, :,1), resmtx)
        writedlm(io,resmtx,delim2)
        flush(io)    
    end    
    ped =  magicgeno.magicped.designinfo
    if !isa(ped,Pedigree) 
        @error "designinfo is not Pedigree type"
        return outfile, nothing
    end
    df = Pedigrees.ped2df(ped)
    df = df[!,[:generation, :member, :gender,:mother,:father]]
    df[!,:gender] .= [if i == "notapplicable" 
        0
    elseif i == "female" 
        1
    elseif i == "male"
        2
    else
        @error string("failed to parse gender=", i)
    end for i in df[!, :gender]]
    rename!(df, ["gender" => "Female=1/Male=2/Hermaphrodite=0"])
    outpedfile = joinpath(workdir, string(outstem,"_mma_ped.csv"))
    open(outpedfile,"w") do io    
        delim2 = ','
        write(io, string("Pedigree-Information, pedinfo\n"))
        CSV.write(io, df; delim=delim2, writeheader=true,append=true)
        write(io, "Pedigree-Information, offspringinfo\n")
        offinfo  = magicgeno.magicped.offspringinfo[!,1:2]
        code = join(1:nfounder,"-")
        insertcols!(offinfo,3,:FunnelCode => code) 
        CSV.write(io, offinfo; delim=delim2, writeheader=true,append=true)
    end
    outfile,outpedfile
end
