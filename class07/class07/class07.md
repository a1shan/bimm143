class07
================
Alyssa Shan
November 1, 2018

``` r
#Revisit function things

source("http://tinyurl.com/rescale-R")

rescale(1:10)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

``` r
#rescale(c(1,5,"sting"))

#rescale2(c(1,5,"string"))
## build both_na function
#1. define examples
x <- c(1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)

#2. find NA (search on google)
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
is.na(y)
```

    ## [1]  TRUE FALSE  TRUE FALSE FALSE

``` r
is.na(x) & is.na(y)
```

    ## [1] FALSE FALSE  TRUE FALSE FALSE

``` r
#3. count # of both NA
sum(is.na(x))
```

    ## [1] 2

``` r
sum(is.na(y))
```

    ## [1] 2

``` r
sum(is.na(x) & is.na(y))
```

    ## [1] 1

``` r
##put working snippet into funciton
both_na <- function(x,y) {
  sum(is.na(x) & is.na(y))
}

both_na(x,y)
```

    ## [1] 1

``` r
## eejit proofing
x <- c(NA, NA, NA)
y1 <- c( 1, NA, NA)
y2 <- c( 1, NA, NA, NA)
both_na(x,y1)
```

    ## [1] 2

``` r
both_na(x,y2)
```

    ## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
    ## shorter object length

    ## [1] 3

``` r
both_na2 <-  function(x, y) {
  ## Check for NA elements in both input vectors and don't allow re-cycling 
  if(length(x) != length(y)) {
    stop("Input x and y should be vectors of the same length", call.=FALSE)
  }
  sum( is.na(x) & is.na(y) )
}

#both_na2(x,y2)

##consolidate
#both_na3 <- function(x, y) {
  ## Print some info on where NA's are as well as the number of them 
  #if(length(x) != length(y)) {
    #stop("Input x and y should be vectors of the same length", call.=FALSE)
  #na.in.both <- ( is.na(x) & is.na(y) )
  #na.number  <- sum(na.in.both)
  #na.which   <- which(na.in.both)
#}
  
both_na3(x,y1)
```

    ## Found 2 NA's at position(s):2, 3

    ## $number
    ## [1] 2
    ## 
    ## $which
    ## [1] 2 3

``` r
##find matching genes (intersect of two dfs)
# simplify as definitions
x <- df1$IDs
y <- df2$IDs
```
