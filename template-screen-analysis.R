#============================================#
#---------> USER DEFINED VARIABLES <---------#
#============================================#
screen.name='pswat'
# Full plate definition file (with all plates)
plate.def.path = "/media/elmicro/meta/PromoterSWAT/screens/R1.1/bin/platedef_promoterSWAT-clean2.csv"
# No need for "/" at the end
screen.path = "/media/elmicro/meta/PromoterSWAT/screens/R1.1" 
output.path = paste0(screen.path,"/results")
channels.order = c("GFP","RFP","BFP")
#============================================#

general   = attach(what = NULL,name = 'mylib::general')
multiwell = attach(what = NULL,name = 'mylib::multiwell')
toolbox   = attach(what = NULL,name = 'mylib::toolbox')
yeastprot = attach(what = NULL,name = 'mylib::yeastprot')
sys.source("/data/elevy/70_R_Data/bin/RToolBox/general.R",general)
sys.source("/data/elevy/70_R_Data/bin/RToolBox/Multiwell.R",multiwell)
sys.source("/data/elevy/70_R_Data/bin/RToolBox/MicroscopeToolBox.R",toolbox)
sys.source("/data/elevy/70_R_Data/bin/RToolBox/yeastProteome.R",yeastprot) 

# Remove all local packages attached and local environments
#mypkg <- grep(search(),pattern = 'mylib',v=T)
#lapply(mypkg, detach, character.only = TRUE, unload = TRUE)
#myenv = gsub(pattern="mylib::",x=mypkg,repl='')
#rm(list = myenv)



library(gtools)
my.timestamp = format(Sys.time(), "%F_%Hh%M")
# Check plate directories by looking for name containing "plate_X" (X = number of plate)
plates.path = mixedsort( grep(x=list.dirs(path=screen.path,full.names = T,recursive = F),pattern = 'plate_',ignore.case = T,value = T) )
# Get image directories of each plate for path ending with "images" without case-sensitivity)
plates.img.dirs = mixedsort(grep(list.dirs(plates.path,recursive = F,full.names = T), pattern='Images$', ignore.case = T, value=T))

if( !file.exists(plate.def.path) ){ stop("Plate definition CSV file not found !") }
#if( !all( check.plates.path ) ){ stop(sprintf("Image directory for plate(s) [%s] not found !",toString(which(!check.plates.path)))) }

design.screen = microscope.get.design(
  F=plates.img.dirs,
  D=seq(plates.path),
  PDEF=plate.def.path,
  FORMAT=384,
  OUT="_output",
  CHANELS=channels.order,
  MEASURE.COL = "_int_b5",
  KEY=c(),
  DIR.res = output.path
)

data.screen.raw = microscope.load.data(design.screen) 
uscope.count.cells(data.screen.raw)
lsos() ## check how much memory the object takes
design.screen   = uscope.process.estimate.background(data.screen.raw,design.screen)
data.screen.raw = uscope.process.reorder(data.screen.raw, design=design.screen)
data.screen.raw = uscope.process.remove.first.pic(data.screen.raw)
data.screen     = uscope.process.remove.background(data.screen.raw, design.screen)

uscope.count.cells(data.screen)

data.screen    = uscope.process.remove.small(data.screen, MIN.size=800,MAX.size=2000)
data.screen2   = uscope.process.BF(data.screen)
data.screen2   = uscope.process.remove.BF.outliers(data.screen2, cutoff=0.8)
data.screen3   = uscope.process.add.ncells(data = data.screen2)
uscope.count.cells(data.screen3)
init.result.folders(design.screen)

### Let's save the image data here.
save(list=c("data.screen.raw","design.screen","data.screen3"),file=paste0(design.screen$DIR.res,"/screen_",screen.name,"-imagedata_",my.timestamp,".Rdata"))
## load(paste0(design.screen$DIR.res,"screen_",screen.name,"-raw-design-data.Rdata"))

#
## First we make a couple of diagnostic plots to make sure things are OK 
##
## BACKGROUND INTENSITY ~ PLATE 
pdf(file= paste(design.screen$DIR.res,"plate_diag_GFP_int_b9_new.pdf",sep=""),width=20, height=15)
diagnostic.intensity(data=data.screen3, design=design.screen, col.of.interest="GFP_int_b9")
dev.off()
## CELL HIGHEST INTENSITY ~ PLATE 
pdf(file= paste(design.screen$DIR.res,"plate_diag_GFP_int_b1_new.pdf",sep=""),width=20, height=15)
diagnostic.intensity(data=data.screen3, design=design.screen, col.of.interest="GFP_int_b1$")
dev.off()
## NUMBER OF CELLS ~ PLATE 
pdf(file= paste(design.screen$DIR.res,"plate_diag_ncells_new.pdf",sep=""),width=20, height=15)
diagnostic.intensity(data=data.screen3, design=design.screen, col.of.interest="ncells$", fun2use=length)
dev.off()

## 
## Output of cell id for each picture (rquired input for Montage)
##
uscope.get.cellid = function(dat,des,p){
  if( !(p %in% names(dat)) ){ cat(sprintf('plate %s do not exist! %s\n',p)); return(NA); }
  tmp  = lapply(dat[[p]],'[',c('cell','pic'))
  for(W in names(tmp)){
    well2snum = get.sNUM(des,p,W)
    tmp[[W]]$snum = well2snum[as.numeric(tmp[[W]]$pic)]
  }
  res = do.call(rbind,tmp)
  res$well = des$WELLS[[p]][ res$snum ]
  res$plate = p
  return(res)
}

## JAVASCRIPT OBJECT FORMAT (key-value associative array = input of imageJ Montage script)
library(rjson)
for( p in 1:length(plates.path) ){
  cellid = uscope.get.cellid(data.screen3,design.screen,p=p)
  CELLS =  split(cellid$cell,f=cellid$snum)
  fname = sprintf("%s/plate_%s-SELECTED-CELLID-SNUM.JSON",plates.path[p],p)
  print(fname)
  FILE=file(fname, open = "w")  
  CELLS.JSON = toJSON(CELLS)
  writeLines(text = gsub(CELLS.JSON,pat='\\],\\"',repl='\\],\n\\"'), con = FILE)
  close(FILE)
}

options(stringsAsFactors = FALSE)

## 
## Now gets info on ORF abundance and localization (input required for Webserver)
##
source("/data/elevy/70_R_Data/bin/RToolBox/yeastProteome.R") 
#library(org.Sc.sgd.db)
#library(Biostrings)
REF = get.REF()
SC  = get.proteome.table()
SC.loc = get.loc.maya(REF)

SC.loc$cyto = 0 + (SC.loc$Control.Localization == "cytosol")

SC.loc$punc = 0 + (SC.loc$Control.Localization == "punctate" |
                     SC.loc$DTT.Localization == "punctate" |
                     SC.loc$H2O2.Localization == "punctate" |
                     SC.loc$Starvation.Localization == "punctate" 
)

# ##
# ## General localization info
# ## --> 1 Manual GO annotation = score of 100
# ## --> 1 high-throughput anno = score of 10
# ## --> 1 prediction annotation= score of 1
# ## So a score of 212 means two manual annotations, 1 HT, and two prediction for a particular loc
# ##
# SC.GO = load.sc.GO()
# ## Localization info in the form of a matrix
# SC.GO.CC = format.GO(SC.GO, REF=REF, type="CC", min=20)
# ## Merges all the localizations that are "cytoplasm" and removes non-informative ones.
# ribo.index = grep("ribos|complex",colnames(SC.GO.CC))
# SC.GO.CC = SC.GO.CC[,-ribo.index]
# SC.GO.CC.filt = SC.GO.CC[, -which( colnames(SC.GO.CC) %in% c("cytoplasm", "intracellular","cytosol", "plasma membrane enriched fraction", "soluble fraction", "cytoplasmic vesicle","polysome"))]
# SC.GO.CC.filt = cbind(SC.GO.CC.filt, cyto=rowSums(SC.GO.CC[, which( colnames(SC.GO.CC) %in% c("cytoplasm", "intracellular","cytosol", "plasma membrane enriched fraction", "soluble fraction", "cytoplasmic vesicle","polysome"))]))
# SC.GO.candidates = rowSums(SC.GO.CC.filt>10)==1 & SC.GO.CC.filt[,c("cyto")] > 10

# Summarize useful info in single dataframe
scinfo = merge(SC,SC.loc)

## 
## Output of image analysis (required for webserver)
##

# MAKE A UNIQUE ID FOR EACH WELL IN FULL PLATE DEFINITION AND PLATE DEFINITION SPLITTED BY PLATE
# To check if GFP, RFP and BFP are duplicated columns :
# duplicated( t(design.screen$PDEF[,c('GFP','RFP','BFP')]) )
design.screen$PDEF$ID = apply(design.screen$PDEF[,c('plate','well','GFP')],1,function(x){ paste(x,collapse='-') })
design.screen$PDEF.split = split(design.screen$PDEF,design.screen$PDEF$plate)
empty.plate = data.frame(well=get.384.by.row(), ord=1:384, stringsAsFactors=FALSE)
## Now makes sure that lines in each plate match a well, even if empty or if no information is available.
for( each.plate in names(design.screen$PDEF.split)){
  design.screen$PDEF.split[[each.plate]] = merge(empty.plate, design.screen$PDEF.split[[each.plate]], by="well", all.x=TRUE)
}

# Mean, Noise, Median for bin 5 of each fluorescent channel
noise = function(x){ var(x) / mean(x)^2 }
NC  = uscope.process.get.one.stat.allp(design = design.screen, data = data.screen3, col.of.interest=c("GFP_int_b0$"), info=c(channels.order,"ID","well",'plate'),fun2use = length)
colnames(NC) = c('ncells',paste0('id.',channels.order),'id.ID','id.well','id.plate')
NC$SCREEN='PromoterSWAT'
NC$PLATE=paste0('plate_',NC$plate)
NC$CHANNELS=paste0(c("BF",channels.order),collapse='-')

for( chan in channels.order ){
  #cat(sprintf("Channel %s\n",chan))
  AVG = uscope.process.get.one.stat.allp(design = design.screen, data = data.screen3, col.of.interest=sprintf('%s_int_b5$',chan), info=c(chan,"ID","well",'plate'), fun2use = mean)
  colnames(AVG)[1]=paste0('avg',chan)
  VAR = uscope.process.get.one.stat.allp(design = design.screen, data = data.screen3, col.of.interest=sprintf('%s_int_b5$',chan), info=c(chan,"ID","well",'plate'), fun2use = var)
  colnames(VAR)[1]=paste0('var',chan)
  EPS = uscope.process.get.one.stat.allp(design = design.screen, data = data.screen3, col.of.interest=sprintf('%s_int_b5$',chan), info=c(chan,"ID","well",'plate'),fun2use = noise)
  colnames(EPS)[1]=paste0('eps',chan)
  MED = uscope.process.get.one.stat.allp(design = design.screen, data = data.screen3, col.of.interest=sprintf('%s_int_b5$',chan), info=c(chan,"ID","well",'plate'),fun2use = median)
  colnames(MED)[1]=paste0('med',chan)
  
  TMP = Reduce( function(x,y) merge(x,y,all=TRUE,sort=FALSE), list(AVG,VAR,EPS,MED) )
  
  NC = merge(NC,TMP,all=TRUE, sort=FALSE)
}


channel.cols = unlist(lapply(channels.order,function(x){ paste0(c("id.","avg","eps","var","med"),x) }))
idcol = paste0('id.',c("ID","well",'plate'))
STATS = NC[,c('SCREEN','PLATE','CHANNELS',idcol,'ncells',channel.cols)]

scinfo.cols = list(
  sel=c("ORF","Gene","UPROT","len","loc.ymd","loc.dtt","loc.h2o2","loc.starv","cyto","punc","over.tox","mrna.wang","ab.pax","ab.GFP.maya.ymd","ab.GFP.maya.dtt","ab.GFP.maya.h2o2","ab.GFP.maya.starv"),
  ab.prot=c("prot12","ab.ms.ratio","ab.ms.haploid","ab.ms.diploid","ab.apex.YPD","ab.apex.YMD","ab.GFP.YPD","ab.GFP.YMD","ab.western","ab.tap"),
  ab.mrna=c("mrna12","mrna.sage","mrna.HDA","cod.cai","cod.bias","cod.fop","cod.tAI","pest.1","pest.2"),
  compo=c("ribo.occu", "len2","diso05","diso1","diso2","sticky05","sticky1","sticky2","KRord","KRdes"),
  chrom=c("ord","chr","ch-start","ch-end","strand"),
  pheno=c("viable","pheno.giav.Var1","pheno.giav.Freq","pheno.growth.dudley.Var1","pheno.growth.dudley.Freq","pheno.resist.dudley.Var1","pheno.resist.dudley.Freq","pheno.growth.up.Var1","pheno.growth.up.Freq","pheno.growth.down.Var1","pheno.growth.down.Freq","pheno.growth.stab.Var1","pheno.growth.stab.Freq","pheno.resist.up.Var1","pheno.resist.up.Freq","pheno.resist.down.Var1","pheno.resist.down.Freq","pheno.resist.stab.Var1","pheno.resist.stab.Freq"),
  loc = c("Description","Control.Median","Control.STD","Control.Localization","DTT.Median","DTT.STD","DTT.FoldChange","Significance.DTT","DTT.Localization","H2O2.Median","H2O2.STD","H2O2.FoldChange","Significance.H2O2","H2O2.Localization","Starvation.Median","Starvation.STD","Starvation.FoldChange","Significance.Starvation","Starvation.Localization","Pup2.DAmP","MA3.Median","MA3.STD","MA3.FoldChange","Significance.MA3","MA3.Localization")
)

SC.screen = merge(x=STATS, y=scinfo[,scinfo.cols$sel], by.y="ORF", by.x="id.GFP", all.x=TRUE,sort=FALSE)
SC.screen$ORF = SC.screen$id.GFP

GENE = SC.screen$Gene
ORF  = SC.screen$ORF
# CORRECT GENE NAME (ONLY FOR VALID ORF - protein encoded in nuclear or mitochondrial chromosome)
norf = sort(table(ORF),decreasing = T)
nucORF = grepl(ORF, pattern="(Y[A-P][LR][0-9]{3}[WC])",ignore.case = T)
mitORF = grepl(ORF, pattern="(Q[0-9]+)",ignore.case = T)
validORF = nucORF | mitORF
SC.screen$NAME[validORF & !is.na(GENE)] = GENE[validORF & !is.na(GENE)]
SC.screen$NAME[validORF & is.na(GENE)] = ORF[validORF & is.na(GENE)]
SC.screen$Gene = NULL

# DATAFRAME WITH ALL INFOS FOR WEBSERVER
SC.screen.ord = SC.screen[,c('ORF','NAME','UPROT','len',colnames(STATS),"loc.ymd","loc.dtt","loc.h2o2","loc.starv","cyto","punc","over.tox","mrna.wang","ab.pax","ab.GFP.maya.ymd","ab.GFP.maya.dtt","ab.GFP.maya.h2o2","ab.GFP.maya.starv")]
save(SC.screen.ord,file = paste0(output.path,"/screen_",screen.name,"-","stats_",my.timestamp,".Rdata"))
## CSV FORMAT
write.table(file = paste0(design.screen$DIR.res,"/screen_",screen.name,"-stats.csv"),x = format(SC.screen.ord, trim = TRUE, digits = 3, nsmall = 4, width = 5,sci=TRUE),sep=';',quote = FALSE,row.names = FALSE)
## JSON FORMAT
library(rjson)
O2G = split(SC.screen.ord,f = SC.screen.ord$ORF)
O2G$Gene = O2G$NAME;
O2G.JSON = toJSON(O2G)
FILE=file(paste0(screen.path,"/",screen.name,"-by-ORF.JSON"), open = "w")
writeLines(text = gsub(O2G.JSON,pat='\\"\\},\\"',repl='\\"\\},\n\\"'), con = FILE)
close(FILE)
