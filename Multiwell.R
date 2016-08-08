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
