###
### Set of general functions to work with multiwell plates
###


get.384.by.col = function(){
    letters = rep(LETTERS[1:16],rep(16))
    numbers = gsub(" ","0",format(1:24,width=2))
    numbers = rep(numbers,rep(16,24))
    return ( paste( letters, numbers, sep="") )
}

get.96.by.col = function(){
    letters = rep(LETTERS[1:8],rep(8))
    numbers = gsub(" ","0",format(1:12,width=2))
    numbers = rep(numbers,rep(8,12))
    return ( paste( letters, numbers, sep="") )
}

get.384.by.row = function(){
    letters = rep(LETTERS[1:16],rep(24,16))
    numbers = rep(gsub(" ","0",format(1:24,width=2)),16)
    return ( paste( letters, numbers, sep="") )    
}

get.96.by.row = function(){
    letters = rep(LETTERS[1:8],rep(12,8))
    numbers = rep(gsub(" ","0",format(1:12,width=2)),8)
    return ( paste( letters, numbers, sep="") )    
}

get.384.by.liha = function(){
    letters = rep(LETTERS[1:16],rep(16))
    numbers = gsub(" ","0",format(1:24,width=2))
    numbers = rep(numbers,rep(16,24))
    ind.384.let = paste( letters, numbers, sep="")
    ind.tmp = c( seq(1,16,by=2), seq(2,16,by=2) )    
    ind.all = rep(ind.tmp,24) + 16*rep(0:23,rep(16,24))    
    return ( ind.384.let[ind.all] )
}

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


##
## Takes a list of proteins and creates a file to reaaray them from a reference arrangement to a new arrangment.
##
## - LIST = list of TRUE / FALSE matching lines of the "plate definition"
##
## - REF  = Table containing plate/position of strains to be re-arrayed
##        --> Two important cols are *plate* and *well*)
##
## - FILENAME = name of file where re-arraying information is written
##
TECAN.rearray.list = function(LIST = grepl("A", REF[,2]), REF=TABLE, FILENAME="~/array_test_rn.txt" ){

    sub.mat = REF[which(LIST),]
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

    final.coord = cbind(all.lines, all.print)
    
    ## Preparing file for writing
    header = paste(c("platesrc",paste("src",1:8, sep=""), "platedst","welldst"))

    final.coord.format = final.coord
    for (i in 1:NCOL(final.coord) ){
        final.coord.format[,i] = gsub(" ","",format(as.numeric(final.coord[,i]), nsmall=6))
    }
    final.coord.format = rbind(header,final.coord.format)
    
    write.table(file=paste(FILENAME, ".tab",sep=""), x=final.coord.format, quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\r\n", sep=",")  

    ## Now takes care of the new definition file for this library.
    ##
    plate.table = table(sub.mat$Plate.Number)
    plate.index = as.numeric(names(table(sub.mat$Plate.Number)))
    K=1
    L=1
    final.REF = c()
    EMPTY.REF = cbind( paste( rep("empty", 384),1:384, sep=""), rep("empty", 384), rep(16,384), "", rep(0,384), rep(0,384), rep(0,384), rep(0,384),rep(0,384), rep(0,384), rep(0,384), rep(0,384) , rep(0,384) )
    for (each.plate in plate.index){

        final.REF = rbind(final.REF, sub.mat[ which(sub.mat$Plate.Number == each.plate),])
        
        modulo = 8-( as.vector(plate.table[K]) %% 8)
        if(modulo < 8){
            for (n in 1:modulo){
                
                final.REF = rbind(final.REF, c(paste( "empty", L,sep=""), paste( "empty", L,sep=""), each.plate, rep(0,25)) )
                L=L+1
            }
        }
        K=K+1
    }

    nplates = as.integer(NROW(final.REF)/1536)
    final.REF2 = final.REF[ 1:(nplates*1536),]
    final.REF2$plate = rep(1:(nplates*4), rep(384, nplates*4) )
    final.REF2$pos   = get.384.by.liha()
    write.csv(file=paste(FILENAME, "_TARGET.csv",sep=""), x=final.REF2)      
}
