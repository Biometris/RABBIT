
function shift_outliers(calledgenofile, outgenofile, delgenofile; 
    tukeyfence::Real=2, 
    softthreshlikeparam::SoftThreshLikeParam= SoftThreshLikeParam(), 
    threshlikeparam::ThreshLikeParam= ThreshLikeParam(),     
    commentstring::AbstractString="##",    
    workdir::AbstractString=pwd())
    filels = [calledgenofile, outgenofile, delgenofile]
    calledgenofile2, outgenofile2, delgenofile2 = getabsfile.(workdir, filels)
    if !isfile(calledgenofile2)
        @error string("calledgenofile not found: ", calledgenofile)
    end
    if !isfile(delgenofile2)
        @error string("delgenofile not found: ",delgenofile)
    end
    outliers = get_outlier_markers(calledgenofile2; tukeyfence, softthreshlikeparam, threshlikeparam, commentstring)    
    if isempty(outliers)
        mv(calledgenofile2, outgenofile2; force=true)
    else
        nheader = MagicBase.findlastcomment(calledgenofile2; commentstring) + 1    
        in_open, out_open, del_open = [begin 
            inext = last(MagicBase.split_allext(f))
            inext in [".vcf",".vcf.gz"] || @error string("genofile ext must be .vcf or .vcf.gz")
            inext in [".vcf.gz"] ? GZip.open : open        
        end for f in filels]    
        in_open(calledgenofile2,"r") do inio                    
            out_open(outgenofile2, "w") do outio                
                del_open(delgenofile2,"a") do delio
                    shift_outlier_line(inio, outio, delio; nheader, outliers)            
                end
            end
        end    
        rm(calledgenofile2)
    end
    outliers
end


function shift_outlier_line(inio::IO,outio::IO,delio::IO; nheader, outliers::AbstractVector)
    # nheader: line no. for the header row
    for _ in 1:nheader
        line = readline(inio,keep=true)        
        write(outio, line)        
    end
    while !eof(inio)        
        rowstring = readline(inio,keep=false)
        # vcf colnames = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL","FILTER", "INFO", "FORMAT"]
        if isempty(outliers)
            isoutlier = false
        else
            markerid = strip(split(rowstring, "\t")[3])
            isoutlier = markerid in outliers
        end
        if isoutlier
            write(delio,rowstring,"\n")
        else
            write(outio,rowstring,"\n")
        end
    end
end    

function get_outlier_markers(mapfile::AbstractString;
    tukeyfence::Real=2, 
    softthreshlikeparam::SoftThreshLikeParam= SoftThreshLikeParam(), 
    threshlikeparam::ThreshLikeParam= ThreshLikeParam(), 
    missingstring::AbstractString="NA",
    commentstring::AbstractString="##",
    workdir::AbstractString = pwd())
    markermap = readmarkermap(mapfile; del_ungrouped = false,        
        commentstring,missingstring, workdir)          
    outliermarkers = get_outlier_markers(markermap; tukeyfence, softthreshlikeparam, threshlikeparam)
    outliermarkers
end


function get_outlier_markers(markermap::MagicBase.AbstractDataFrame; 
    tukeyfence::Real=2, 
    softthreshlikeparam::SoftThreshLikeParam= SoftThreshLikeParam(), 
    threshlikeparam::ThreshLikeParam= ThreshLikeParam(), 
    )    
    outls = Int[]
    colls = [:foundererror,:offspringerror,:baseerror,:allelicbias, :allelicbias_2,:allelicoverdispersion,:allelicdropout]
    for col = colls
        if col == :allelicbias_2
            errname = :allelicbias
            errls = 1.0 .- maximum.(markermap[!,errname]) # remove low-outlier            
        else
            errname = col
            errls = maximum.(markermap[!,errname]) # remove up-outlier
        end
        allequal(errls) && continue 
        softthresh = getproperty(softthreshlikeparam,errname)        
        hardthresh = getproperty(threshlikeparam,errname)
        if errname in [:allelicbias,:allelicbias_2]
            softthresh = max(softthresh,1.0 - softthresh)
            hardthresh = max(hardthresh,1.0 - hardthresh)
        end
        delindices = del_error_indices2(errls; 
            tukeyfence, errorscale = errname == :allelicoverdispersion ? log : logit,
            softthresh, hardthresh, 
        )				
        append!(outls,delindices)
    end
    sort!(union!(outls))
    markermap[outls,:marker]
end


function del_error_indices2(errorls::AbstractVector;	
	tukeyfence::Real=2,
	softthresh=0.05,
	hardthresh=0.25,
	errorscale::Function)	
	ls = errorscale.(errorls)
    q1,q3 = quantile(ls,[0.25,0.75])
	upbound = q3+tukeyfence*(q3-q1)	    	
	hardthresh2 = errorscale == logit ? min(1.0,hardthresh) : hardthresh
	upbound = min(max(upbound,errorscale(softthresh)),errorscale(hardthresh2))
    findall(ls .> upbound)
end
