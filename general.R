##
## General functions that are neither related to plotting nor handling of real data
##

## Makes sure not to print to many lines in the buffer
options(max.print=4000)

######################################################################
## Color scales
######################################################################
BROWN2YELLOW = c()
for (white.tone in c(0,80)){
    BROWN2YELLOW  = c(BROWN2YELLOW, rgb( seq(150+white.tone,150+white.tone,len=5)/255, seq(60+white.tone,120+white.tone,l=5)/255, seq(30+white.tone,30+white.tone,l=5)/255)  )#, seq(0,0.95,l=20)), c( seq(1,0,l=20), seq(0,0,l=20)))
}

KROGAN.COL  = function(N=20){
  return(  rgb( c( seq(0,0,l=N), seq(0,0.95,l=N)), c( seq(0.7,0,l=N), seq(0,0.95,l=N)), c( seq(1,0,l=N), seq(0,0,l=N))) )
}

B2O.COL = colorRampPalette(colors=c("black","orange"))


######################################################################
## General Matrix/Vector functions 
######################################################################
## Will rotate a matrix so that the display by "image" looks the same.
rotate.matrix = function(M){
   return(t(M[nrow(M):1,,drop=FALSE]))
}



## MEAN function that tolerate strings to use in aggregates.
##
my.MEAN = function(x){

    if(is.double(x[1])){
        return(mean(x,na.rm=TRUE))
    } else {
        return(x[1])
    }
}
## MEDIAN function that tolerate strings to use in aggregates.
##
my.MEDIAN = function(x){

    if(is.double(x[1])){
        return(median(x,na.rm=TRUE))
    } else {
        return(x[1])
    }
}


######################################################################
## General Memory check functions 
######################################################################
free.memory = function(){
  gc()
}


memory.status = function(N=10){
  lsos(n=N)
}

# improved list of objects
.ls.objects <- function (pos = 1, pattern, order.by,
                        decreasing=FALSE, head=FALSE, n=5) {
    napply <- function(names, fn) sapply(names, function(x)
                                         fn(get(x, pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x)
                        as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.dim)
    names(out) <- c("Type", "Size", "Rows", "Columns")
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    if (head)
        out <- head(out, n)
    out
}
# shorthand
lsos <- function(..., n=10) {
    .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}

