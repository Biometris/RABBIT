

function fix_nonsegregate(genofile::AbstractString,pedfile::AbstractString; 
    isfounderinbred::Bool=true,
    workdir::AbstractString=pwd(),    
    outstem::Union{Nothing,AbstractString}="outstem",
    verbose::Bool=true)
    magicgeno = formmagicgeno(genofile, pedfile; isfounderinbred,formatpriority=["GT"]);    
    fix_nonsegregate!(magicgeno;isfounderinbred,verbose)
    outfile = getabsfile(workdir, string(outstem,"_fixnonsegregate.vcf.gz"))
    savegenodata(outfile,magicgeno)
end

function fix_nonsegregate!(magicgeno::MagicGeno; isfounderinbred::Bool=true,verbose::Bool=true)
    model = "jointmodel"
    MagicBase.info_magicgeno(magicgeno;verbose)
    magicprior=MagicReconstruct.calmagicprior(magicgeno,model; isfounderinbred)
    missingcode = isfounderinbred ? "N" : "NN"
    fillcode = isfounderinbred ? "1" : "11"
    nchr = length(magicgeno.markermap)    
    @info "impute/fill genotypes for founders whose progeny's genotypes are all missing"
    @info "#foundergeno_impute: number of founder gneotypes are imputed based on cofounders"
    @info string("#foundergeno_fill: number of founder gneotypes are arbitrarily filled to ", fillcode)    
    for chr in 1:nchr
        popmakeup = MagicReconstruct.calpopmakeup(magicgeno,chr,model, magicprior; isfounderinbred)    	
        fprogenls = MagicBase.get_fprogenyls(popmakeup) # list of offspring for each founder
        popfounderls = unique([i["founder"] for i=values(popmakeup)])
        cofounderls = [sort(setdiff(unique(reduce(vcat,popfounderls[in.(i,popfounderls)])),[i])) for i in 1:length(fprogenls)]
        nmarker = size(magicgeno.markermap[chr],1)
        offgeno = magicgeno.offspringgeno[chr]
        fgeno = magicgeno.foundergeno[chr]
        nimpute = nfill = 0
        for snp in 1:nmarker
            ls = [unique(offgeno[snp,offs]) for offs in fprogenls]
            fls_allmiss = findall([i == [missingcode] for i in ls])
            while true
                oldcount = nimpute
                for f in fls_allmiss
                    alleles = setdiff(unique(split(join(fgeno[snp,cofounderls[f]]),"")),["N"])
                    if length(alleles) == 1 && fgeno[snp,f] == missingcode                
                        fgeno[snp,f] = isfounderinbred ? only(alleles) : only(alleles)^2
                        nimpute += 1                
                    end
                end
                oldcount == nimpute && break
            end
            for f in fls_allmiss
                # alleles = setdiff(unique(split(join(fgeno[snp,cofounderls[f]]),"")),["N"])
                if fgeno[snp,f] == missingcode                
                    fgeno[snp,f] = fillcode
                    nfill += 1
                end
            end
        end
        chrid = magicgeno.markermap[chr][1,:linkagegroup]        
        msg = string("chr=",chrid, ", #foundergeno_impute=",nimpute, ", #foundergeno_fill=",nfill)
        verbose && @info msg
    end    
    magicgeno
end

