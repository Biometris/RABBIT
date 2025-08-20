

function plotmarkererror(mapfile::AbstractString;
    tukeyfence::Real=2,
    missingstring=["NA","missing"],
    commentstring::AbstractString="##",
    workdir::AbstractString = pwd())
    markermap = readmarkermap(mapfile; del_ungrouped = false,        
        commentstring,missingstring, workdir)  
    plotmarkererror(markermap;tukeyfence)
end

function plotmarkererror(markermap::AbstractDataFrame;
    tukeyfence::Real=2)
    colls = [:foundererror,:offspringerror,:baseerror,:allelicbias,:allelicoverdispersion,:allelicdropout]
    fencels = tukeyfence*ones(length(colls))
    b = [begin         
        isanymiss = any(ismissing.(markermap[!,i]))         
        if isanymiss 
            true
        else
            ls = unique(round.(unique(markermap[!,i]),digits=4))
            # length(ls) <=1 && println("col=",i, ",unqiue values=",ls)
            length(ls) <=1 
        end
    end for i in colls]
    deleteat!(colls,b)
    deleteat!(fencels,b)
    isempty(colls) && return nothing    
    gls = [begin     
        col = colls[c]        
        whisker_range = fencels[c]        
        # row = [!(ismissing(i) || isnan(i)) for i in markermap[!,col]]
        row = !
        yls = Vector(markermap[row, col])
        yls[yls .== 0.0] .= 1e-6
        if in(col, [:foundererror,:offspringerror,:baseerror,:allelicoverdispersion,:allelicdropout]) 
            @. yls = log10(yls)            
            ylabel = string("log10 ", col)
        else
            # @. yls = log10(abs(yls))
            # ylabel = string("log10 |", col,"|")            
            ylabel = string(col)
        end                     
        xls = string.(markermap[row,:linkagegroup])                
        g =violin(xls, yls;
            fillalpha=0.75, linewidth=0.5,                
            legend=false,            
            ylabel
        )              
        boxplot!(xls, yls;
            fillalpha=0.75, linewidth=2,
            whisker_range, 
            legend=false, 
            title = string(col, ": tukeyfence=",whisker_range)
        )   
        g
    end for c in eachindex(colls)]
    nchr = length(unique(markermap[!,:linkagegroup]))
    plot(gls...;
        layout=(length(gls),1),
        size = (max(1000,nchr*60),300*length(gls)),
        left_margin=20Plots.mm,         
        bottom_margin=15Plots.mm,
    )
end


function plot_peroffspringerror(perofffile::AbstractString;
    tukeyfence::Real=2,    
    workdir::AbstractString = pwd())
    peroffspringerr = CSV.read(getabsfile(workdir, perofffile),DataFrame)    
    plot_peroffspringerror(peroffspringerr; tukeyfence)
end


function plot_peroffspringerror(peroffspringerr::AbstractDataFrame;
    tukeyfence::Real=2)
    whisker_range = tukeyfence
    colnames = names(peroffspringerr)
    colls = findall(occursin.(r"^peroffspringerr",colnames))
    yls = vec(Matrix(peroffspringerr[!,colls]))
    yls[yls .< 1e-6] .= 1e-6
    yls .= log10.(yls)
    labs = replace.(colnames[colls],"peroffspringerror_"=>"")
    xls = repeat(labs,inner=size(peroffspringerr,1))
    fig =violin(xls, yls;
        fillalpha=0.75, linewidth=0.5,                
        legend=false,            
        ylabel = "log10 peroffspringerr"
    )              
    boxplot!(xls, yls;
        fillalpha=0.75, linewidth=2,            
        legend=false, 
        title = string("tukeyfence=",whisker_range)
    )   
    fig
end




function ploterrorcompare(truegeno, magicgeno)
    colls = [:foundererror, :offspringerror,:baseerror,:allelicbias,:allelicoverdispersion,:allelicdropout]
    gls = [begin 
        trueerrls = reduce(vcat, [[ismissing(i) ? i : mean(i) for i in df[!,col]] for df in truegeno.markermap])
        esterrls = reduce(vcat, [[ismissing(i) ? i : mean(i) for i in df[!,col]] for df in magicgeno.markermap])       
        length(trueerrls) == length(esterrls) || error("inconsistent #markers")
        b = .!(ismissing.(trueerrls) .|| ismissing.(esterrls))      
        if any(b) && length(unique(esterrls[b])) > 1
            r = cor(trueerrls[b],esterrls[b])
            mtrue = round(mean(trueerrls[b]),digits=3)
            mest = round(mean(esterrls[b]),digits=3)
            @info string("col=", col, ",mtrue=",mtrue,",mest=",mest)
            g = scatter(trueerrls[b],esterrls[b];
                xlab = string("True ", col),
                ylab = string("Est ", col),
                label=string("r=", round(r,digits=3)),
            )
            g
        else
            nothing
        end
    end for col in colls]
    deleteat!(gls, isnothing.(gls))
    isempty(gls) && return nothing
    nrow = ceil(Int,length(gls)/2)
    ncol = length(gls) > 1 ? 2 : 1
    fig = plot(gls...,
        layout=(nrow,ncol),
        size=(400*ncol,260*nrow),
        leftmargin = 10Plots.mm,
        bottommargin = 10Plots.mm,
    )
    fig
end

