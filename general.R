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
free.memory = function(){ gc() }
memory.status = function(N=10){ lsos(n=N) }
# improved list of objects
.ls.objects <- function (pos = 1, pattern, order.by, decreasing=FALSE, head=FALSE, n=5) {

    napply <- function(names, fn){ sapply(names, function(x){ fn(get(x, pos = pos)) }) }
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x){ as.numeric(dim(x))[1:2] }))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.dim)
    names(out) <- c("Type", "Size", "Rows", "Columns")
    if (!missing(order.by)){ out <- out[order(out[[order.by]], decreasing=decreasing), ] }
    if (head){ out <- head(out, n) }
    return(out)
}
# shorthand
lsos <- function(..., n=10) { .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n) }

#####
### TESTING AND COUNTING VALUES
######
# checks if value is even (0 can be considered as even or not)
is.even <- function(x,exclude0=TRUE){ 
  if( !missing(x) & !is.numeric(x) ){ stop("Need to provide a vector of numerical values") }
  res = (as.integer(x) %% 2 == 0)
  if(exclude0){ res[ x == 0 ] = FALSE }
  names(res) = x
  return(res)
}
# checks if value is odd
is.odd <- function(x){ 
  if( !missing(x) & !is.numeric(x) ){ stop("Need to provide a vector of numerical values") }
  res = (as.integer(x) %% 2 != 0)
  names(res) = x
  return(res)
}
# checks if vector can be interpreted as binary values (0/1 or 2 values)
is.binary =  function(x) {
  if( !missing(x) & !is.numeric(x) ){ 
    if( is.factor(x) & nlevels(x) == 2 ){ return(TRUE) }
    else if( is.factor(x) & nlevels(x) != 2 ){ return(FALSE) }
    else{ stop("Need to provide a vector of numerical values") }
  }
  if( length(unique(x)) == 2 ){ return(TRUE) }else{ return(FALSE)}
}

# Counts number of NA or not-NA values
sum.na <- function (x,notNA=FALSE){ if(is.null(x)){ return(NULL) }else{ sum(is.na(x) == !notNA) } }
# Counts number of not-NA values
sum.notNA = function(x){ sum.na(x,notNA = TRUE) }
# Counts number of unique values
Ulen = function(x){ length(unique(x)) }

#####
### TRIMMING UNWANTED CHARACTERS OR VALUES (also rounding values)
#####
# returns a vector without NaN 
trim.NaN <- function(d){ if( is.null(dim(d)) ){ d[!is.nan(d)] }else{ apply(d,2,trim.NaN) } }
# returns string w/o leading whitespace
trim.blanks.lead <- function (x)  sub("^\\s+", "", x)
# returns string w/o trailing whitespace
trim.blanks.trail <- function (x) sub("\\s+$", "", x)
# returns string w/o leading or trailing whitespace
trim.blanks <- function (x, everywhere = FALSE){ if(everywhere){ gsub("\\s+", "", x) }else{ gsub("^\\s+|\\s+$", "", x) } }

RoundUpToNearest = function(nb, roundto){
  if (roundto == 0){ return(nb) }
  else{ return( ceiling(nb / roundto) * roundto ) }
}

RoundDownToNearest = function(nb, roundto){
  if (roundto == 0){ return (nb) }
  else{ return( floor(nb / roundto) * roundto ) }
}

RoundToNearest = function(nb,roundto){
  a = RoundDownToNearest(nb,roundto);
  b = RoundUpToNearest(nb,roundto);
  
  if( abs(nb-a) < abs(nb-b) ){ return(a) }
  else if( abs(nb-a) > abs(nb-b) ){ return(b) }
  else{ return(nb) }
  
}

# Normalize number between 0 and 1
Norm01 = function(x){ scale(x,center=min(x),scale=diff(range(x)))[,1] }

#####
### GET VALUES IN SPECIFIC POSITIONS IN VECTOR
#####
# returns element before last 
b4last <- function(x) { head( tail(x,n=2) , n = 1) }
# returns vector w/o last element
alast <- function(x) { head(x, n = -1) }
# returns last element of a vector
last <- function(x) { tail(x, n = 1) }
# returns First and Last element of a vector
HandT <- function(x) { c(head(x,n=1),tail(x,n=1)) }

# returns top q% percent of values
ind.top.perc = function(x,my.p=seq(0,1,0.2)){
  q  = quantile(x,probs=my.p)
  return( x >= b4last(q) )
}

# returns down q% percent of values
ind.bottom.perc = function(x,my.p=seq(0,1,0.2)){
  q  = quantile(x,probs=my.p)
  return( x >= q[2] )
}

ind.perc = function(x,step=0.1){
  q  = quantile(x,probs=c(0,step,1))
  print(paste("CUT-OFF (top ",(100*(1-step)),"%) :",q[2]))
  return( x >= q[2] )
}

#####
### FORMAT VALUES INTO SPECIFIC STRINGS
#####
# return values in percentages with the percent symbol
percent = function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}

# return values in integer percentages with the percent symbol
percent0 = function(x, format = "f", ...) {
  percent(x,digits=0,format=format)
}

# return values with the chosen symbols for units
unit = function(x, digits = 0, format = "f", units,...) {
  if( missing(units) ){ print("You need to specify a string for the units") }
  paste0(formatC(x, format = format, digits = digits, ...), as.character(units) )
}

# return values based on their rank with the correct suffix
rank0 = function(x) {
  sapply( x, FUN = function(d){ 
    d.chr = as.character(d); ndigits = nchar(d.chr);
    last.digit = as.numeric(substr(d.chr,ndigits,ndigits))
    twolast.digit = as.numeric(substr(d.chr,ndigits-1,ndigits))
    if( !(last.digit %in% c(1,2,3)) | twolast.digit %in% c(11,12,13) ){ return(paste0(d,"th")) }
    if( last.digit == 1 & twolast.digit != 11){ return(paste0(d,"st")) }
    if( last.digit == 2 & twolast.digit != 12){ return(paste0(d,"nd")) }
    if( last.digit == 3 & twolast.digit != 13){ return(paste0(d,"rd")) }
  }
  )
}

# return values with the right prefix for their units
toUnits = function(x){
  if( x >= 1e+24 ){ return(c(pre="Y",x=format(x/1e+24,digit=1)))       }
  else if( x >= 1e+21  &  x < 1e+24 ){ return(c(pre="Z",x=format(x/1e21,digit=1)))   }
  else if( x >= 1e+18  &  x < 1e+21 ){ return(c(pre="E",x=format(x/1e18,digit=1)))   }
  else if( x >= 1e+15  &  x < 1e+18 ){ return(c(pre="P",x=format(x/1e15,digit=1)))   }
  else if( x >= 1e+12  &  x < 1e+15 ){ return(c(pre="T",x=format(x/1e12,digit=1)))   }
  else if( x >=  1e+9  &  x < 1e+12 ){ return(c(pre="G",x=format(x/1e9,digit=1)))    }
  else if( x >=  1e+6  &  x < 1e+9  ){ return(c(pre="M",x=format(x/1e6,digit=1)))    }
  else if( x >=  1e+3  &  x < 1e+6  ){ return(c(pre="k",x=format(x/1e3,digit=1)))    }
  else if( x >=  1e+2  &  x < 1e+3  ){ return(c(pre="h",x=format(x/1e2,digit=1)))    }
  else if( x >=  1e+1  &  x < 1e+2  ){ return(c(pre="da",x=format(x/1e1,digit=1)))   }
  else if( x <   1e+1  &  x > 1e-1  ){ return(c(pre=" ",x=format(x,digit=1)))        }
  else if( x <=  1e-1  &  x > 1e-2  ){ return(c(pre="d",x=format(x/1e-1,digit=1)))   }
  else if( x <=  1e-2  &  x > 1e-3  ){ return(c(pre="c",x=format(x/1e-2,digit=1)))   }
  else if( x <=  1e-3  &  x > 1e-6  ){ return(c(pre="m",x=format(x/1e-3,digit=1)))   }
  else if( x <=  1e-6  &  x > 1e-9  ){ return(c(pre="mu",x=format(x/1e-6,digit=1)))  }
  else if( x <=  1e-9  &  x > 1e-12 ){ return(c(pre="n",x=format(x/1e-9,digit=1)))  }
  else if( x <= 1e-12  &  x > 1e-15 ){ return(c(pre="p",x=format(x/1e-12,digit=1)))  }
  else if( x <= 1e-15  &  x > 1e-18 ){ return(c(pre="f",x=format(x/1e-15,digit=1)))  }
  else if( x <= 1e-18  &  x > 1e-21 ){ return(c(pre="a",x=format(x/1e-18,digit=1)))  }
  else if( x <= 1e-21  &  x > 1e-24 ){ return(c(pre="z",x=format(x/1e-21,digit=1)))  }
  else if( x <= 1e-24 ){ return(c(pre="y",x=format(x/1e-24,digit=1)))  }
}

# paste values by pairs separated by hypen (two consecutive values forms a pair)
paste.pair = function(vec,s="-"){
  new = c();
  for( i in 1:(length(vec)-1) ){
    new = c(new, paste(vec[i],vec[i+1],sep=s))
  }
  return(new)
}

