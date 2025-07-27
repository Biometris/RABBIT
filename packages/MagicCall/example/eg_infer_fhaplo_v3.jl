
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
outstem = "outstem"
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

exportall(MagicCall)
likeparameters = MagicBase.LikeParam()
formatpriority = ["AD","GT"]
missingset = [".","./.", ".|."]    
isdelmonomorphic = false

rowstring = readline(inio,keep=false)
rowgeno = split(rowstring,"\t")
ismultiallele = length(split(rowgeno[5],",")) > 1 # col5=alternative                     
ismultiallele && isdelmultiallelic && return ("multiallelic", rowstring)
baseerror = MagicBase.get_baseerror(likeparameters)
res = extract_rowgeno!(rowgeno, fcols,offcols,formatpriority,
    missingset,isfounderinbred,isdelmonomorphic,baseerror)    
res
res == "monomorphic" && return ("monomorphic",join(rowgeno,"\t"))
inputformat, fgeno, founderformat, offgeno, offspringformat = res
isoffmiss = rowgeno[offcols] .== "NA"

fgeno
offgeno

popidls = keys(popmakeup)
logl = callogl_singlesite(offgeno; popmakeup, popidls, epso, baseerror,
                allelicbias, allelicoverdispersion, 
                allelicdropout,offspringformat,israndallele
)

if offspringformat == "AD"                
    like = MagicReconstruct.diplolikeGBS(view(offgeno,offls),epso,baseerror,allelicbias,allelicoverdispersion,allelicdropout; israndallele)    
else
    like = MagicReconstruct.diplolike(view(offgeno,offls),epso;israndallele)
end


###########################################################

using Distributions

# imputation/correction and error estimations    
likeparameters=LikeParam(offspringerror=0.005, peroffspringerror=0.0)
priorlikeparameters = PriorLikeParam(baseerror=Beta(1,99))
israndallele = true
isinfererror = true
byfounder = 0
threshfounder = 0.9
@time fhaplo_GT, fhaplo_GP, esterrors = infer_fhaploerror_singlesite(fgeno, founderformat, offgeno,offspringformat; 
    model, popmakeup, israndallele,isinfererror, byfounder, isfounderinbred, 
    samplesize=100, likeparameters, priorlikeparameters,threshfounder,
)    

fhaplo_GT
fhaplo_GP
esterrors

###########################################################

# geno_postprob = fhaplo_postprob
# fhaplo_GT, fhaplo_GP = postprob2vcfgeno(fhaplo_postprob; callthreshold=threshfounder,digits=2)

###########################################################

liketargetls, epsf, epso, pereoffspringerror, baseerror, allelicbias,allelicoverdispersion,allelicdropout = MagicBase.extract_likeparameters(likeparameters)		        
@time fhaploset,fhaploweight = priorfhaplo_singlesite(fgeno,founderformat; foundererror=epsf, 
    baseerror, allelicbias,allelicoverdispersion,allelicdropout, israndallele
)
fhaplo = map((x,y)->x[rand(Categorical(y))], fhaploset, fhaploweight)       
nowlogp = callogl_singlesite(fhaplo, offgeno; popmakeup, epsf,epso, baseerror,
    allelicbias, allelicoverdispersion, allelicdropout,offspringformat,israndallele
) 

