###
### General functions to load function yeast data as well as to benchmark data
### 
library(org.Sc.sgd.db)
library(Biostrings)

##
## Return a matrix where yeast proteins are annotated according to a particular physical position
##
get.REF = function(){

  allGenes = read.csv("/data/elevy/70_R_Data/ALLGENES_R2.csv", stringsAsFactors=FALSE)
   tmp = read.csv("/data/elevy/70_R_Data/DHFR-LibraryOligosPosList.csv")
  allGenes$PCR = tmp$Fragment.size
  
  ## fakes the 16th plate
  fake.16 = cbind( paste( rep(NA, 384),1:384, sep=""), rep(NA, 384), rep(16,384), allGenes[ 1:384, c(4:6)], rep(0,384), rep(0,384), rep(0,384), rep(0,384),rep(0,384), rep(0,384), rep(0,384), rep(0,384) , rep(0,384) )
  colnames(fake.16) = colnames(allGenes)
  allGenes = rbind(allGenes, fake.16)
  imageJ.order = read.table("/data/elevy/62_ADE4/bin/imageJ.order", header=TRUE)
  new.order = order(imageJ.order$plate, imageJ.order$cadran, imageJ.order$pos)
  all.Genes.order = allGenes
  all.Genes.order[new.order,] = all.Genes.order
  REF = cbind(all.Genes.order,imageJ.order)
  REF[ which(REF[,1]==""),1]= paste( rep("Empty",4), 1:4, sep="")
  colnames(REF)[1]= c("ORF")
  REF$Desc = as.vector(unlist(sapply(as.list(REF[,1]), mget, env=org.Sc.sgdDESCRIPTION, ifnotfound=NA)))
  REF$row2 = rep(NA,6144)
  REF$col2 = rep(NA,6144)

  REF$row2[ REF$plate %in% c(1,3,9,11) ] =     REF$row[REF$plate %in% c(1,3,9,11)]*2
  REF$row2[ REF$plate %in% c(2,4,10,12) ] = 1+ (REF$row[REF$plate %in% c(2,4,10,12)]*2)

  REF$row2[ REF$plate %in% c(5,7,13,15) ] =     REF$row[REF$plate %in% c(5,7,13,15)]*2
  REF$row2[ REF$plate %in% c(6,8,14,16) ] = 1+ (REF$row[REF$plate %in% c(6,8,14,16)]*2)

  REF$col2[ REF$plate %in% c(1,3,9,11) ] =     REF$col[REF$plate %in% c(1,3,9,11)]*2
  REF$col2[ REF$plate %in% c(2,4,10,12) ] = 1+ (REF$col[REF$plate %in% c(2,4,10,12)]*2)

  REF$col2[ REF$plate %in% c(5,7,13,15) ] =     REF$col[REF$plate %in% c(5,7,13,15)]*2
  REF$col2[ REF$plate %in% c(6,8,14,16) ] = 1+ (REF$col[REF$plate %in% c(6,8,14,16)]*2)

  REF$P1536 = 0
  REF$P1536[ REF$plate %in% c(1,3,9,11) ] = 1
  REF$P1536[ REF$plate %in% c(2,4,10,12) ] = 2
  REF$P1536[ REF$plate %in% c(5,7,13,15) ] = 3
  REF$P1536[ REF$plate %in% c(6,8,14,16) ] = 4
  
  REF$ord=1:6144
  REF$Prot.Name = as.vector(unlist(sapply(as.list(REF[,1]), mget, env=org.Sc.sgdGENENAME, ifnotfound=NA)))
  REF$Prot.Name[is.na(REF$Prot.Name)]=REF$ORF[is.na(REF$Prot.Name)]
  REF$SGD = as.vector(unlist(sapply(as.list(REF[,1]), mget, env=org.Sc.sgdSGD, ifnotfound=NA)))
  
  return(REF)
}

##
## Return a matrix where yeast proteins are annotated according to a particular physical position
##
get.REF.full = function(){

  allGenes = read.csv("/data/elevy/70_R_Data/MichnickLablist_DHFRPCAcollectionList.csv", stringsAsFactors=FALSE)
  #tmp = read.csv("/data/elevy/70_R_Data/MichnickLablist_DHFRPCAcollectionList.csv")
  allGenes$PCR = tmp$Fragment.size
  
  ## fakes the 16th plate
  fake.16 = cbind( paste( rep(NA, 384),1:384, sep=""), rep(NA, 384), rep(16,384), allGenes[ 1:384, c(4:6)], rep(0,384), rep(0,384), rep(0,384), rep(0,384),rep(0,384), rep(0,384), rep(0,384), rep(0,384) , rep(0,384)  , rep(0,384)  , rep(0,384)  , rep(0,384)  , rep(0,384)  , rep(0,384)  , rep(0,384)  , rep(0,384)  , rep(0,384) )
  colnames(fake.16) = colnames(allGenes)
  allGenes = rbind(allGenes, fake.16)
  imageJ.order = read.table("/data/elevy/62_ADE4/bin/imageJ.order", header=TRUE)
  new.order = order(imageJ.order$plate, imageJ.order$cadran, imageJ.order$pos)
  all.Genes.order = allGenes
  all.Genes.order[new.order,] = all.Genes.order
  REF = cbind(all.Genes.order,imageJ.order)
  REF[ which(REF[,1]==""),1]= paste( rep("Empty",4), 1:4, sep="")
  colnames(REF)[1]= c("ORF")
  REF$Desc = as.vector(unlist(sapply(as.list(REF[,1]), mget, env=org.Sc.sgdDESCRIPTION, ifnotfound=NA)))
  REF$row2 = rep(NA,6144)
  REF$col2 = rep(NA,6144)

  REF$row2[ REF$plate %in% c(1,3,9,11) ] =     REF$row[REF$plate %in% c(1,3,9,11)]*2
  REF$row2[ REF$plate %in% c(2,4,10,12) ] = 1+ (REF$row[REF$plate %in% c(2,4,10,12)]*2)

  REF$row2[ REF$plate %in% c(5,7,13,15) ] =     REF$row[REF$plate %in% c(5,7,13,15)]*2
  REF$row2[ REF$plate %in% c(6,8,14,16) ] = 1+ (REF$row[REF$plate %in% c(6,8,14,16)]*2)

  REF$col2[ REF$plate %in% c(1,3,9,11) ] =     REF$col[REF$plate %in% c(1,3,9,11)]*2
  REF$col2[ REF$plate %in% c(2,4,10,12) ] = 1+ (REF$col[REF$plate %in% c(2,4,10,12)]*2)

  REF$col2[ REF$plate %in% c(5,7,13,15) ] =     REF$col[REF$plate %in% c(5,7,13,15)]*2
  REF$col2[ REF$plate %in% c(6,8,14,16) ] = 1+ (REF$col[REF$plate %in% c(6,8,14,16)]*2)

  REF$ord=1:6144
  REF$Prot.Name = as.vector(unlist(sapply(as.list(REF[,1]), mget, env=org.Sc.sgdGENENAME, ifnotfound=NA)))
  REF$SGD = as.vector(unlist(sapply(as.list(REF[,1]), mget, env=org.Sc.sgdSGD, ifnotfound=NA)))
  
  return(REF)
}


###
### Creates a matrix containing many descriptors of the yeast proteome.
###
get.proteome.table = function(REF=get.REF()){

  ORD = data.frame( ORF = as.vector(REF[,1]), ord=1:6144)
  
  SC.abund = read.csv("/data/elevy/70_R_Data/TABLES/sc_abund_christine.csv")[,1:16]
  colnames(SC.abund) = c("ORF","ab.apex.YPD","ab.apex.YMD","ab.GFP.YPD","ab.GFP.YMD","ab.western","ab.tap","mrna.sage","mrna.HDA","mrna.wang","mrna.av","cod.cai","cod.bias","cod.fop","pest.1","pest.2")
  SC.abund = SC.abund[,-c(11)]

  SC.tAI = read.table("/data/elevy/70_R_Data/TABLES/sc.tAI.values.txt", sep="\t", head=FALSE)
  colnames(SC.tAI) = c("ORF","cod.tAI")
  SC.tAI[,2]= round(SC.tAI[,2],3)
  
  SC.proteom = read.table("/data/elevy/48_subho_sc/data/Sc_expression_matrix.dat", sep="\t", header=TRUE)
  colnames(SC.proteom) = c("ORF","len","cod.tAI2","mrna.sde","ab.sde","mrna12","prot12","ribo.occu")
  SC.proteom = SC.proteom[,-c(3,4,5)]

  SC.addition = read.table("/data/elevy/70_R_Data/TABLES/sc.ABUND.paxDB.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)  
  colnames(SC.addition) = c("ORF","ab.pax")
  SC.addition[,2] = round(as.numeric(SC.addition[,2]),2)
  ## Subho has been using the rna HDA and the westerm for prots.
  
  SC.data = merge(SC.abund, SC.tAI, all.x=TRUE)
  SC = merge(SC.data,SC.proteom)
  SC = merge(SC,SC.addition, all.x=TRUE)

  SC.maya = read.csv("/data/elevy/70_R_Data/TABLES/sc.ab.GFP.maya.csv")[,c(1,3,5,6,9,10,13,14,17)];
  colnames(SC.maya) = c("ORF","ab.GFP.maya.ymd","loc.ymd","ab.GFP.maya.dtt","loc.dtt","ab.GFP.maya.h2o2","loc.h2o2","ab.GFP.maya.starv","loc.starv")
  #index.to.replace = which(SC.maya[, grep("ab", colnames(SC.maya)) ] < 0)
  #SC.maya[index.to.replace, grep("ab", colnames(SC.maya)) ] = NA
  SC = merge(SC,SC.maya, all.x=TRUE)

  
  SC.dosage = read.csv("/data/elevy/70_R_Data/TABLES/sc.DOSAGE.sopko.txt", header=TRUE)
  SC.dosage = data.frame(ORF=  SC.dosage[,1], over.tox= 1 )
  SC = merge(SC,SC.dosage, all.x=TRUE)
  SC$over.tox[is.na(SC$over.tox)]=0

  ## intensity.H = Haploid??
  SC.mann = read.csv("/data/elevy/70_R_Data/TABLES/sc.MannAbund.csv")
  SC.mann.ins = SC.mann[,c(1,34,38,39)]
  colnames(SC.mann.ins)=c("ORF","ab.ms.ratio","ab.ms.haploid","ab.ms.diploid")
  SC = merge(SC,SC.mann.ins, all.x=TRUE)
  #SC$over.tox[is.na(SC$over.tox)]=0
  
  SC.diso = read.table("/data/elevy/70_R_Data/TABLES/sc.diso.txt", header=FALSE)
  colnames(SC.diso) = c("ORF", "len2", "diso05","diso1","diso2","sticky05","sticky1","sticky2","KRord","KRdes")
  SC = merge(SC,SC.diso, all.x=TRUE)

  SC.genom = read.csv("/data/elevy/70_R_Data/TABLES/sc.GENOMFEATURES.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
  colnames(SC.genom) = c("SGDID","ORF","chr","ch-start","ch-end","strand","chr2")
  SC.genom = SC.genom[,-c(1,7)]
  SC.genom$chr = gsub("chr", "", SC.genom$chr)
  SC = merge(SC, SC.genom, all.x=TRUE, by=c(1))

  #SC.phenom = read.csv("/data/elevy/70_R_Data/TABLES/sc.multidrug.table.csv")[,-c(2,6,8,9,10:88)];
  #colnames(SC.phenom) = c("ORF","essential","slow.hete","slow.homo","num.cond")
  #SC = merge(SC, SC.phenom, all.x=TRUE, by=c(1))  

  ### There are a few lines with B" and the " need to be deleted
  SC.phenom = read.csv("/data/elevy/70_R_Data/TABLES/sc.genotypes.SGD.131008.tsv", sep="\t");
  colnames(SC.phenom) = c("SGD.num", "ORF","PROT","PROT.synonym","ORF.qualif","exp.type", "disturbance.type", "readout.type", "readout.outcome", "notsure", "strain", "drug", "media","remark","remark2","pmid","publi")
  SC.phenom = SC.phenom[, -c(1,3,4,5,10,14,15)]

  ### Viable / inviable 
  SC.phenom1 = data.frame(ORF = SC.phenom$ORF[ SC.phenom[,9]==12140549 & SC.phenom[,8]==""& SC.phenom[,7]==""],
      viable = as.vector(SC.phenom$readout.type[ SC.phenom[,9]==12140549 & SC.phenom[,8]==""& SC.phenom[,7]==""]),
      stringsAsFactors=FALSE)
  SC.phenom1$viable[which(SC.phenom1$viable=="Viable")]=1
  SC.phenom1$viable[which(SC.phenom1$viable=="Inviable")]=0
  SC.phenom1$viable[SC.phenom1$viable != 1 & SC.phenom1$viable != 0]= NA

  SC = merge(SC, SC.phenom1, all.x=TRUE, by=c(1))
  
  ### slow growth Giaver 2002
  SC.phenom2 = table(SC.phenom$ORF[SC.phenom[,9]==12140549 & SC.phenom[,8]=="" & SC.phenom[,5]=="decreased"])
  SC.phenom2 = SC.phenom2[SC.phenom2 > 0]
  SC.phenom2 = data.frame( ORF = names(SC.phenom2), pheno.giav = SC.phenom2, stringsAsFactors=FALSE)

  SC = merge(SC, SC.phenom2, all.x=TRUE, by=c(1))
  
  ### slow growth Dudley / carbon source
  SC.phenom3 = table(SC.phenom$ORF[SC.phenom$pmid==16729036 & SC.phenom[,4]== "Utilization of carbon source"])
  SC.phenom3 = SC.phenom3[SC.phenom3 > 0]
  SC.phenom3 = data.frame( ORF = names(SC.phenom3), pheno.growth.dudley = SC.phenom3, stringsAsFactors=FALSE)

  SC = merge(SC, SC.phenom3, all.x=TRUE, by=c(1))

  ### slow growth Dudley / drug resistance
  SC.phenom4 = table(SC.phenom$ORF[SC.phenom$pmid==16729036 & SC.phenom[,4]== "Resistance to chemicals"])
  SC.phenom4 = SC.phenom4[SC.phenom4 > 0]
  SC.phenom4 = data.frame( ORF = names(SC.phenom4), pheno.resist.dudley = SC.phenom4, stringsAsFactors=FALSE)

  SC = merge(SC, SC.phenom4, all.x=TRUE, by=c(1))

  ### slow growth all --> pheno.growth.less / pheno.growth.more / pheno.growth.same
  SC.pheno.growth.increase = table(SC.phenom$ORF[SC.phenom[,4]=="Competitive fitness" & SC.phenom[,5]=="increased" & grepl("S288", SC.phenom[,6])] )
  SC.pheno.growth.decrease = table(SC.phenom$ORF[SC.phenom[,4]=="Competitive fitness" & SC.phenom[,5]=="decreased" & grepl("S288", SC.phenom[,6])] )
  SC.pheno.growth.normal   = table(SC.phenom$ORF[SC.phenom[,4]=="Competitive fitness" & SC.phenom[,5]=="normal"    & grepl("S288", SC.phenom[,6])] )

  SC.pheno.growth.increase = SC.pheno.growth.increase[SC.pheno.growth.increase>0]
  SC.pheno.growth.decrease = SC.pheno.growth.decrease[SC.pheno.growth.decrease>0]
  SC.pheno.growth.normal   = SC.pheno.growth.normal  [SC.pheno.growth.normal  >0]

  SC.pheno.growth.increase = data.frame(ORF = names(SC.pheno.growth.increase),  pheno.growth.up   = SC.pheno.growth.increase, stringsAsFactors=FALSE )
  SC.pheno.growth.decrease = data.frame(ORF = names(SC.pheno.growth.decrease),  pheno.growth.down = SC.pheno.growth.decrease, stringsAsFactors=FALSE )
  SC.pheno.growth.normal   = data.frame(ORF = names(SC.pheno.growth.normal  ),  pheno.growth.stab = SC.pheno.growth.normal  , stringsAsFactors=FALSE )

  SC = merge(SC, SC.pheno.growth.increase, all.x=TRUE, by=c(1))
  SC = merge(SC, SC.pheno.growth.decrease, all.x=TRUE, by=c(1))
  SC = merge(SC, SC.pheno.growth.normal, all.x=TRUE, by=c(1))
  
  ### drug resistance --> pheno.resist.less / pheno.resist.more / pheno.resist.same
  SC.pheno.resist.increase = table(SC.phenom$ORF[SC.phenom[,4]=="Resistance to chemicals" & SC.phenom[,5]=="increased" & grepl("S288", SC.phenom[,6])] )
  SC.pheno.resist.decrease = table(SC.phenom$ORF[SC.phenom[,4]=="Resistance to chemicals" & SC.phenom[,5]=="decreased" & grepl("S288", SC.phenom[,6])] )
  SC.pheno.resist.normal   = table(SC.phenom$ORF[SC.phenom[,4]=="Resistance to chemicals" & SC.phenom[,5]=="normal"    & grepl("S288", SC.phenom[,6])] )

  SC.pheno.resist.increase = SC.pheno.resist.increase[SC.pheno.resist.increase >0]
  SC.pheno.resist.decrease = SC.pheno.resist.decrease[SC.pheno.resist.decrease >0]
  SC.pheno.resist.normal   = SC.pheno.resist.normal  [SC.pheno.resist.normal   >0]
  
  SC.pheno.resist.increase = data.frame(ORF = names(SC.pheno.resist.increase),  pheno.resist.up   = SC.pheno.resist.increase , stringsAsFactors=FALSE)
  SC.pheno.resist.decrease = data.frame(ORF = names(SC.pheno.resist.decrease),  pheno.resist.down = SC.pheno.resist.decrease , stringsAsFactors=FALSE)
  SC.pheno.resist.normal   = data.frame(ORF = names(SC.pheno.resist.normal  ),  pheno.resist.stab = SC.pheno.resist.normal   , stringsAsFactors=FALSE)

  SC = merge(SC, SC.pheno.resist.increase, all.x=TRUE, by=c(1))
  SC = merge(SC, SC.pheno.resist.decrease, all.x=TRUE, by=c(1))
  SC = merge(SC, SC.pheno.resist.normal, all.x=TRUE, by=c(1))

  ### homozygote hap
  ##SC = merge(SC, SC.phenom, all.x=TRUE, by=c(1))  
  
  SC.ord = merge(ORD, SC, all.x=TRUE, by=c(1))  
  SC = SC.ord[ order(SC.ord$ord),]  
  return(SC)
}
