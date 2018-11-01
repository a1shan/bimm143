#Revisit function things
source("http://tinyurl.com/rescale-R")

rescale(1:10)
rescale(c(1,5,"sting"))

rescale2(c(1,5,"string"))
## build both_na function
#1. define examples
x <- c(1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)

#2. find NA (search on google)
is.na(x)
is.na(y)
is.na(x) & is.na(y)

#3. count # of both NA
sum(is.na(x))
sum(is.na(y))
sum(is.na(x) & is.na(y))

##put working snippet into funciton
both_na <- function(x,y) {
  sum(is.na(x) & is.na(y))
}

both_na(x,y)

## eejit proofing
x <- c(NA, NA, NA)
y1 <- c( 1, NA, NA)
y2 <- c( 1, NA, NA, NA)
both_na(x,y1)
both_na(x,y2)

both_na2 <-  function(x, y) {
  ## Check for NA elements in both input vectors and don't allow re-cycling 
  if(length(x) != length(y)) {
    stop("Input x and y should be vectors of the same length", call.=FALSE)
  }
  sum( is.na(x) & is.na(y) )
}

both_na2(x,y2)

##consolidate
both_na3 <- function(x, y) {
  ## Print some info on where NA's are as well as the number of them 
  if(length(x) != length(y)) {
    stop("Input x and y should be vectors of the same length", call.=FALSE)
  na.in.both <- ( is.na(x) & is.na(y) )
  na.number  <- sum(na.in.both)
  na.which   <- which(na.in.both)
}
  
both_na3(x,y2)


##find matching genes (intersect of two dfs)
# simplify as definitions
x <- df1$IDs
y <- df2$IDs