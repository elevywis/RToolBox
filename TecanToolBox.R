library(gdata)
library(grofit)

Sys.setlocale('LC_ALL','C') ## text manipulation gave error messages with this

##
## Functions used to manipulate data created by the Tecan Infinite Reader
##
## There are also functions to feed the TecanEvo for cherry picking.
##                  
##

##
## Creates an experimental design for data acquired with the infinite reader
##
infinite.design = function(
    DIR = "/media/elusers/data/infinite_reader/Udi/codon_diploids_01_test/",
    repli = c("R1"),
    cond  = c("SD"),
    CH=c("OD","GFP","RFP"),
    nplates = 16,
    format=384,
    file.name = "growth_",
    direction = c()
    )  {
    
    design = list()
    design$DIR   = DIR
    design$CH   = CH
    design$repli = repli
    design$cond  = cond
    design$n.plates = nplates
    design$format  = format
    design$direction = direction
    design$n.repli = length(repli)
    design$n.cond  = length(cond)
    design$n.CH  = length(CH)
    design$file.name =file.name
    files = dir(DIR)
    files = files[grep(paste("^",file.name,sep=""), files, perl=TRUE)]

    design$modulo = design$n.plates * design$n.cond * design$n.repli
    design$n.timepoints = floor(length(files)/design$modulo)
    files = files[1: (design$modulo * design$n.timepoints ) ]
    design$n.files = length(files)
    design$files = files    
    return(design)
}


## Returns data created by the infinite reader.
## The data has the following form:
##
## - DATA$OD$1/2/3/4 ... ==>series of OD, and the name of each series is the time in minutes.
## - DATA$C1$1/2/3/4 ... ==>series of fluorescence reads for Chanel 1
##
get.infinite.data = function(DIR="/media/elusers/data/infinite_reader/elevy/test2/", pattern="^yeast.*07-2[67]-", types=c("OD","GFP"), R=2){

  files = dir(DIR)
  files = files[grep(pattern, files, perl=TRUE)]

  time = get.time.from.date2(files)

  N = 384
  T = length(time)
  
  all.res = list()
  for (i in 1:length(types) ){
    all.res[[i]] = matrix(NA, ncol = T, nrow=N)
    colnames(all.res[[i]]) = time
  }
  names(all.res) = types

  #types.index = c(17,50)
  types.index = c(18,51)

  index = 1
   
  for (each.file in files){

    print(each.file)
    
    tmp = read.xls(xls=paste(DIR,each.file,sep=""))
    
    for (i in 1:length(types) ){
      r   = as.numeric(as.matrix(tmp[types.index[i]:(types.index[i]+15),2:25]))
      all.res[[i]][,index] = r                 
    }

    index = index+1
  }
  return(all.res)
}



## Returns data created by the infinite reader.
## The data has the following form:
##
## - DATA$OD$1/2/3/4 ... ==>series of OD, and the name of each series is the time in minutes.
## - DATA$C1$1/2/3/4 ... ==>series of fluorescence reads for Chanel 1
##
load.infinite.data = function(EXP){

    res = list()    
    for ( i in 1:EXP$n.read ){
        res[[i]] = list()
        for ( j in 1:EXP$n.cond ){
            res[[i]][[j]] = list()
        }
        names(res[[i]]) = EXP$cond
    }
    names(res) = EXP$read

    ## processes exp condition
    ##
    index = 1
    #types.index = c(18,51)

    for( l in 1:EXP$n.plates ){

        for ( k in 1:EXP$n.repli ) {

            for ( j in 1:EXP$n.cond ){
            
                current.files = EXP$files[ seq( index, by=EXP$modulo, length=EXP$n.timepoints) ]
                tmp.res = list()
                
                ## First, builds a plate series (384 * timepoints
                for (file.i in 1:length(current.files) ){
                    each.file = current.files[file.i]                    
                    print(paste("Cond:", EXP$cond[j], " Repli:", EXP$repli[k], " Plate:", l, "file:",each.file))

                    tmp = read.xls(xls=paste(EXP$DIR,each.file,sep=""), stringsAsFactors=FALSE)
                    types.index = grep("Temperature",tmp[,2])+2
                    
                    for ( i in 1:EXP$n.read ){
                        
                        r = as.numeric(as.matrix(tmp[types.index[i]:(types.index[i]+15),2:25]))
                        
                        if(file.i==1){
                            tmp.res[[i]] = matrix( r , ncol=1)
                        } else {
                            tmp.res[[i]] = cbind(tmp.res[[i]], r)
                        }
                    }
                }
                ## Second, concatenate the series to other plates
                for ( i in 1:EXP$n.read ){
                    colnames(tmp.res[[i]]) = EXP$timepoints

                    #indexes = paste(rep(LETTERS[1:16],24), rep(gsub(" ","0",format(1:24, width=2)), rep(16,24)),sep="")        
                    if(length(EXP$direction) >= k & EXP$direction[k] == 1){
                        
                        row.index = paste(l, "-", paste(rep(LETTERS[16:1],24), rep(24:1, rep(16,24)),sep=""),sep="" )
                        rownames(tmp.res[[i]]) = row.index
                        tmp.res[[i]] = tmp.res[[i]][ order(row.index),]

                    } else {
                        row.index = paste(l, "-", paste(rep(LETTERS[1:16],24), rep(1:24, rep(16,24)),sep=""), sep="")
                        rownames(tmp.res[[i]]) = row.index
                        tmp.res[[i]] = tmp.res[[i]][ order(row.index),]
                        
                    }
                    
                    if( l==1){
                        res[[i]][[j]][[k]] = tmp.res[[i]]
                    } else {
                        res[[i]][[j]][[k]] = rbind( res[[i]][[j]][[k]], tmp.res[[i]])
                    }
                }
                index = index + 1
            }
        }
    }
    return(res)
}

## Computes the corresponding times in minutes from the following format
## 
##
add.time = function(design){

    data = design$files
    data = gsub(gsub(paste("^.*", design$file.name,sep=""), "", x=data, perl=TRUE), pattern=".xlsx", repl="")
    data = strsplit(data,"-| ")
    data = t( sapply(data,function(x){ return(as.numeric(x))} ) )
    
    time      = rep(0,NROW(data))
    
    for ( i in 1:(NROW(data)) ){
        
        delta = 1440 * ( data[i,3] - data[1,3]  ) +
            60 * ( data[i,4] - data[1,4] ) +
                1 * (  data[i,5] - data[1,5] )
        
        time[i] = delta
    }
    design$alltimes = time
    design$time = time[seq(1, len=KO.design$n.timepoints, by=KO.design$modulo)]
    return(design)
}

## Computes the corresponding times in minutes from the following format
##
##
get.time.from.date2 = function(files){

  data = files
  data = gsub(gsub("^.*001-", "", x=data, perl=TRUE), pattern=".xls", repl="")
  data = strsplit(data,"-")
  data = t( sapply(data,function(x){ return(as.numeric(x))} ) )
  
  time      = rep(0,NROW(data))
  
  for ( i in 1:(NROW(data)) ){

    delta = 1440 * ( data[i,3] - data[1,3]  ) +
      60 * ( data[i,4] - data[1,4] ) +
        1 * (  data[i,5] - data[1,5] )
    
    time[i] = delta
  }

  return(time)
}

##
## graphs the OD and adds the fluorescence
##
plot.a.strain = function(strain.no, data, cond=1, rep=1, main, design){
    ## First, plots the OD
    plot(design$time, data[[design$CH[1]]][[cond]][[rep]][strain.no,], "l", xlab="Time", ylab="OD")
    par(new=TRUE)
    plot(design$time, data[[design$CH[2]]][[cond]][[rep]][strain.no,], col="dark green", "l", axes=FALSE)
    axis(4)
    par(new=TRUE)
    plot(design$time, data[[design$CH[3]]][[cond]][[rep]][strain.no,], col="dark red","l", axes=FALSE)    
}

##
## graphs the OD and adds the fluorescence
##
plot.a.strain.duplex = function(strain.no, data, cond=1, rep=1, main, design){
    par(mai=c(1,1,1,1))
    ## First, plots the OD
    plot(design$time, data[[design$CH[1]]][[cond]][[rep]][strain.no,], "l", xlab="Time (minutes)", axes=FALSE, ylab="RFP", main=main)
    axis(1)
    par(new=TRUE)
    plot(design$time, data[[design$CH[2]]][[cond]][[rep]][strain.no,], col="yellow ", "l", axes=FALSE, xlab="", ylab="", lwd=3)
    axis(4)
    par(new=TRUE)
    plot(design$time, data[[design$CH[3]]][[cond]][[rep]][strain.no,], col="dark red","l", axes=FALSE, xlab="", ylab="", lwd=3)    
    axis(2)
    abline(h=c(1000,2500), lty=2, col="dark red")
    abline(h=c(1000,2500), lty=2, col="dark red")
}

##res = infer.fluo(L1=KO.data.ori

###
### This will fit a spline to each read-series and save parameters in a list that has the same structure as the original data
###
get.spline.params = function(data, design, smooth=0.4){

    my.res = data
    
    for( each.CH in names(data)){

        for( each.cond in names(data[[each.CH]])){

            for ( each.repli in 1:length(data[[each.CH]][[each.cond]]) ) {
                
                my.res[[each.CH]][[each.cond]][[each.repli]] = matrix(NA, nrow=NROW(data[[each.CH]][[each.cond]][[each.repli]]), ncol=6)
                colnames(my.res[[each.CH]][[each.cond]][[each.repli]])=c("min","max","rate.spline","lag.spline","rate.fit","lag.fit")
                
                for (each.strain in 1:NROW(data[[each.CH]][[each.cond]][[each.repli]]) ){
                    
                    fit= gcFitSpline(time=design$time, data = KO.data.ori[[each.CH]][[each.cond]][[each.repli]][each.strain,], control=grofit.control(smooth.gc=smooth))
                    
                    ## intersect, final value, rate-spline, lag-spline, rate-fit, lagfit
                    my.res[[each.CH]][[each.cond]][[each.repli]][each.strain,] = c(
                                                                    predict(fit$spline, c(design$time[1],design$time[length(design$time)] ))$y,
                                                                    fit[[6]][[2]],
                                                                    fit[[6]][[3]],
                                                                    fit[[7]][[2]],
                                                                    fit[[7]][[3]]
                                                                    )                    
                }
            }            
        }
    }
    return(my.res)
}

###
### This will fit a spline to each read-series and save parameters in a list that has the same structure as the original data
###
get.spline.values = function(data, design, smooth=0.4){

    my.res = data
    
    for( each.CH in names(data)){

        for( each.cond in names(data[[each.CH]])){

            for ( each.repli in 1:length(data[[each.CH]][[each.cond]]) ) {
                                
                for (each.strain in 1:NROW(data[[each.CH]][[each.cond]][[each.repli]]) ){
                    
                    fit= gcFitSpline(time=design$time, data = data[[each.CH]][[each.cond]][[each.repli]][each.strain,], control=grofit.control(smooth.gc=smooth))
                    my.res[[each.CH]][[each.cond]][[each.repli]][each.strain,] = 
                        predict(fit$spline, design$time)$y
                }
            }            
        }
    }
    return(my.res)
}

###
### This takes the original data and the "fitted splines" data, and substract the "background", which is the original value
### It also makes sure that no value is below 0, otherwise it's set to 0
###
substract.initial.value = function(data.ori, data.fit){
    
    for( each.CH in names(data.ori)){

        for( each.cond in names(data.ori[[each.CH]])){

            for ( each.repli in 1:length(data.ori[[each.CH]][[each.cond]]) ) {

                data.ori[[each.CH]][[each.cond]][[each.repli]] =
                    data.ori[[each.CH]][[each.cond]][[each.repli]] - data.fit[[each.CH]][[each.cond]][[each.repli]][1]

                data.ori[[each.CH]][[each.cond]][[each.repli]][ data.ori[[each.CH]][[each.cond]][[each.repli]]<0] = 0
            }
        }
    }
    return(data.ori)
}

##
## takes a list of two matrices:
##
## L1 --> rows = protein, cols = time
## L2 --> rows = protein, cols = time
##
## Fits a curve to both, L1 and L2, and returns values from L2 for particular L1 values.
##
infer.fluo = function(L1, L2, od.values=c(0.05,0.08, 0.1, 0.15, 0.2, 0.25, 0.27, 0.30, 0.32, 0.35), smooth=0.4, design){

    ### Very important to normalize the fluo?!
    #norma = rowMeans(L2[,1:4])
    #norma.mat = matrix( rep(norma, NCOL(L2)) , ncol=NCOL(L2), byrow=FALSE)
    #L2 = L2 - norma.mat
    #L2[L2<0]=0

    gc.time = design$time
    
    fits1 = apply(L1, 1, function(x) { gcFitSpline(time=gc.time, data = x, control=grofit.control(smooth.gc=smooth) )  } )
    fits2 = apply(L2, 1, function(x) { gcFitSpline(time=gc.time, data = x, control=grofit.control(smooth.gc=smooth) )  } )

    time.end   =   gc.time[ length(gc.time)]
    time.start =   gc.time[1]
  
    preds1 = sapply(1:length(fits1), function(x) {
        return (predict(fits1[[x]]$spline, seq(time.start, time.end, len= 1000)) )
    } )
    
  preds2 = sapply(1:length(fits2), function(x) {
      return (predict(fits2[[x]]$spline, seq(time.start, time.end, len= 1000)) )
  } )

  preds1x = preds1[1,1]$x
  preds1y = matrix( unlist(preds1[2,]), ncol=length(preds1x), nrow=length(fits2), byrow=TRUE)
  preds2x = preds2[1,1]$x
  preds2y = matrix( unlist(preds2[2,]), ncol=length(preds2x), nrow=length(fits2), byrow=TRUE)

  all.times = c()
  
  for( each.value in od.values){
  
    res = unlist(apply(preds1y, 1, function(x){
        
        pos = which.min(abs(x- each.value))
        if(length(pos) > 0){
            return(pos[length(pos)])
        } else {
            return(1)
        }
    }))
    all.times = cbind(all.times, res)
  }

  fluo.res = c()
  
  for (i in 1:NCOL(all.times)){

    res = sapply(1:NROW(L2), function(x){
      if(all.times[x,i] > 1){
        return(preds2y[x,all.times[x,i] ])
      } else {
        return(NA)
      }
    })

    fluo.res = cbind(fluo.res, res)
  }
  colnames(fluo.res) = od.values
  return(fluo.res)
}


indexes.2.positions = function(indexes=names(table(REF$Index.384)[1:10]), plate.def=get.384.by.col()){
    return( as.vector(unlist(sapply( indexes, index.2.position, plate.def=get.384.by.col() ) ) ))
}

index.2.position = function(x, plate.def=get.384.by.col() ){    
    return ( which(plate.def==x))
}

##
## This function takes a list of positions for a reference matrix, and creates a file that can
## be interpreted by the robot to pick specific positions.
##
## The reference matrix must contain the following columns (Name and CASE must match exactly!!), from which positions to pick will be determined:
##    - ORF       --> Gene ID
##    - Prot.Name --> Gene Name
##    - plate     --> plate info (1, 2 ,3 ...)
##    - Index.384 --> well info  (A01 ...)
##
TECAN.rearray.list = function(LIST = grepl("A", REF[,2]), REFMAT, FILENAME="~/array_test_rn.txt" ){

    if( sum( colnames(REFMAT) %in% c("ORF", "Prot.Name", "plate", "Index.384") )<4 ){
        print("Problem, the column names of the reference matrix should contain the following: ORF, Prot.Name, plate, Index.384")
        return(0)
    }
    
    sub.mat = REFMAT[which(LIST),]
    sub.mat$well = indexes.2.positions(indexes = sub.mat$Index.384)
    
    print(paste("There are ", length(table(sub.mat$plate))," plates to pick from"));

    new.order = order( sub.mat$plate, sub.mat$well)

    sub.mat = sub.mat[new.order,]

    line = 1
    plate.change = c(0,diff(sub.mat$plate))

    all.lines = c()
    all.print = c()
    
    dest.pos = 1;
    dest.plate = 1;
    FIRST = 1;
    
    while(line < (NROW(sub.mat))    ){        

        if(dest.pos == 385){
            dest.plate = dest.plate+1
            dest.pos   = 1
            FIRST = 1
        }
        
        if(FIRST){
            print.coord = c(dest.plate, dest.pos)
            FIRST=0
        } else {
            print.coord = c(dest.plate, dest.pos+1)
            FIRST=1
            dest.pos = dest.pos + 16
        }        

        all.print = rbind(all.print, print.coord)       
        picking.coordinates = c(sub.mat$plate[line])        
        next.batch = plate.change[line:(line+7)]
        
        ## if a plate.change exists in the next batch, I will only pick the N colonies
        ## that are on the current plate.
        if( sum(next.batch, na.rm=TRUE) > 0){
            
            INDEX.max = which(next.batch>0 & !is.na(next.batch)  )[1]            
            if(INDEX.max == 1){
                INDEX.max=8
            } else {
                INDEX.max=INDEX.max-1
            }
            
        } else {
            INDEX.max= max( which(next.batch==0) )
        }   
        for (line.tmp in 0:7){
            if( (line.tmp < INDEX.max) & (line+line.tmp) <= NROW(sub.mat) )  { ##  &  ((line + line.tmp) < 385)
                                        # push normal pickup coordinate
                picking.coordinates  = c(picking.coordinates,  sub.mat$well[ line + line.tmp])             
            } else {
                                        # push EMPTY pickup coordinate                
                picking.coordinates = c(picking.coordinates, 2000)
            }
        }
        all.lines = rbind(all.lines, picking.coordinates)

        ## coordinates are filled up.
        ## I will print them normally, and increment LINE only with what has been REALLY picked up (i.e., 
        line = line+INDEX.max
    }

    ## Concatenates the "source info" with the "target info".
    final.coord = cbind(all.lines, all.print)

    ## Last line added to tell the robot that it is the end.
    footer = c(max(sub.mat$plate)+1, rep(2000,8), dest.plate, 1)

    #print("Last line=")
    #print(footer)
    
    final.coord = rbind(final.coord, footer)

    ## Preparing file for writing
    header = paste(c("platesrc",paste("src",1:8, sep=""), "platedst","welldst"))

    final.coord.format = final.coord
    for (i in 1:NCOL(final.coord) ){
        final.coord.format[,i] = gsub(" ","",format(as.numeric(final.coord[,i]), nsmall=6))
    }
   
    
    final.coord.format = rbind(header,final.coord.format)
    
    write.table(file=FILENAME, x=final.coord.format, quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\r\n", sep=",")  

    ## Now takes care of the new definition file for this library.
    ##
    plate.table = table(sub.mat$plate)
    plate.index = as.numeric(names(table(sub.mat$plate)))
    K=1
    L=1
    final.REF = c()
    for (each.plate in plate.index){

        final.REF = rbind(final.REF, sub.mat[ which(sub.mat$plate == each.plate),])
        
        modulo = 8-( as.vector(plate.table[K]) %% 8)
        if(modulo < 8){
            for (n in 1:modulo){
                emptyline = final.REF[1,]
                emptyline[ colnames(emptyline) == "ORF" ]=paste( "empty", L,sep="")
                emptyline[ colnames(emptyline) == "Prot.Name" ]=paste( "empty", L,sep="")
                emptyline[ colnames(emptyline) == "plate" ]=each.plate
                emptyline[ ! colnames(emptyline) %in% c("plate","Prot.Name","ORF") ]=0                
                final.REF = rbind(final.REF, emptyline )
                L=L+1
            }
        }
        K=K+1
    }

    final.REF$pos      = rep(get.384.by.liha(),10)[1:NROW(final.REF)]
    final.REF$plateN   = rep(1:20,rep(384,20))[1:NROW(final.REF)]

    #nplates = as.integer(NROW(final.REF)/384)    
    #final.REF2 = final.REF[ 1:(nplates*384),]
    #final.REF2$plate = rep(1:(nplates), rep(384, nplates) )
    write.csv(file=paste(FILENAME, "_TARGET.csv",sep=""), x=final.REF)
    return(final.REF)
}



###
### Returns a list of positions corresponding to a 384 well plate
###
### --> goes by columns (A01, B01 ..... A02 ...)
###
get.384.by.col = function(){
    letters = rep(LETTERS[1:16],rep(16))
    numbers = gsub(" ","0",format(1:24,width=2))
    numbers = rep(numbers,rep(16,24))
    return ( paste( letters, numbers, sep="") )
}

###
### Returns a list of positions corresponding to a 96 well plate
###
### --> goes by columns (A01, B01 ..... A02 ...)
###
get.96.by.col = function(){
    letters = rep(LETTERS[1:8],rep(8))
    numbers = gsub(" ","0",format(1:12,width=2))
    numbers = rep(numbers,rep(8,12))
    return ( paste( letters, numbers, sep="") )
}

###
### Returns a list of positions corresponding to a 384 well plate
###
### --> goes by lines (A01, A02, ..... B01 ...)
###
get.384.by.row = function(){
    letters = rep(LETTERS[1:16],rep(24,16))
    numbers = rep(gsub(" ","0",format(1:24,width=2)),16)
    return ( paste( letters, numbers, sep="") )    
}

###
### Returns a list of positions corresponding to a 384 well plate
###
### --> goes by the order of the liquid handler (groups of 8 positions - positions are
###    *NON-CONSECUTIVE* and go every second well)
###
### [1] "A01" "C01" "E01" "G01" "I01" "K01" "M01" "O01" "B01" "D01" "F01" "H01"
### [13] "J01" "L01" "N01" "P01" "A02" "C02" "E02" "G02" "I02" "K02" "M02" "O02"
###
get.384.by.liha = function(){
    letters = rep(LETTERS[1:16],rep(16))
    numbers = gsub(" ","0",format(1:24,width=2))
    numbers = rep(numbers,rep(16,24))
    ind.384.let = paste( letters, numbers, sep="")
    ind.tmp = c( seq(1,16,by=2), seq(2,16,by=2) )    
    ind.all = rep(ind.tmp,24) + 16*rep(0:23,rep(16,24))    
    return ( ind.384.let[ind.all] )
}

###
### Returns a list of 96 positions matching one quadrant of a 384 well plate
###
get.384.byQuad = function(Q = 1){

    if(Q==1){
        numbers = seq(1,24,by=2)
        letters = LETTERS[ seq(1,16,by=2) ]
    }
    if(Q==2){
        numbers = seq(2,24,by=2) 
        letters = LETTERS[ seq(1,16,by=2) ]
    }
    if(Q==3){
        numbers = seq(1,24,by=2) 
        letters = LETTERS[ seq(2,16,by=2) ]
    }
    if(Q==4){
        numbers = seq(2,24,by=2) 
        letters = LETTERS[ seq(2,16,by=2) ]
    }

    letters.all = rep(letters,rep(12))
    numbers.all = gsub(" ","0",format(numbers,width=2))
    numbers.all = rep(numbers.all,rep(8,12))
    
    ind.384.let = paste( letters.all, numbers.all, sep="")
    return ( ind.384.let )
}
