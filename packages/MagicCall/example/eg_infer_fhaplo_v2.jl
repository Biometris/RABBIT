
using Revise
using MagicBase
using MagicCall
cd(@__DIR__)
pwd()

isfounderinbred = false
dataid = "sim"
genofile=string(dataid,"_magicsimulate_geno.vcf.gz")
pedfile = string(dataid, "_magicsimulate_ped.csv")

pedinfo = pedfile
outstem = ""
outext = ".vcf"
workdir = pwd()
commentstring = "##"

###########################################################
using GZip
outstem = ""
outstem *= "_magiccall"
logio = open(outstem*".log","w")
inio = GZip.open(genofile,"r") 
out_open = outext in [".vcf.gz",".csv.gz"] ? GZip.open : open
outfile = outstem*"_geno"*outext
outfile2 = getabsfile(workdir,outfile)            
delmarkerfile = outstem*"_delmarker.vcf.gz"
delmarkerfile2 = getabsfile(workdir,delmarkerfile)            
outio = out_open(outfile2, "w") 

lastcomment = MagicBase.findlastcomment(genofile; commentstring)
nheader = lastcomment + 1
for _ in 1:nheader-1
    line = readline(inio,keep=true)        
    write(outio, line)
end


using MagicReconstruct
model = "jointmodel"

# parse title row
titlerow = split(readline(inio,keep=false),"\t")
magicped, fcols, offcols,newtitlerow = MagicBase.parse_titlerow(titlerow, pedinfo;
    keepvcf=true, commentstring,workdir,logio)        
write(outio,newtitlerow,"\n")
# write(delio,newtitlerow,"\n")
isautosome = true    
prioraa = MagicReconstruct.calprior(magicped,model; isfounderinbred,isautosome)    
popmakeup = MagicReconstruct.calpopmakeup(magicped, model,prioraa; isfounderinbred,isautosome)
nstate, nfgl = MagicReconstruct.hmm_nstate_nfgl(popmakeup)    

###########################################################


MagicBase.exportall(MagicCall)
likeparameters = MagicBase.LikeParam()
formatpriority = ["AD","GT"]
missingset = [".","./.", ".|."]    
isfounderinbred,isdelmonomorphic = false, false

rowstring = readline(inio,keep=false)
rowgeno = split(rowstring,"\t")
ismultiallele = length(split(rowgeno[5],",")) > 1 # col5=alternative                     
ismultiallele && isdelmultiallelic && return ("multiallelic", rowstring)

exportall(MagicCall)

res = transform_rowgeno!(rowgeno, fcols,offcols,likeparameters,formatpriority,
    missingset,isfounderinbred,isdelmonomorphic)    
res

inputformat, fgeno, offgeno, offspringformat = res
fixedfounders = get_fixedfounders(magicped,offgeno,offspringformat)     
fhaplo, esterrors = infer_singlesite(fgeno,fixedfounders, offgeno,offspringformat; 
    model, popmakeup, israndallele,isfounderinbred, isinfererror, byfounder, 
    likeparameters, priorlikeparameters
)    




fhaploset,fhaploweight = priorfhaplo_singlesite(fgeno,fformat,fixedfounders; 
    genoerror=0.005,misserror=0.1, allowmissing=true,allowcorrect=true
)      