# Can you improve this analysis code?
library(bio3d)
s1 <- read.pdb("4AKE") # kinase with drug
s2 <- read.pdb("1AKE") # kinase no drug
s3 <- read.pdb("1E4Y") # kinase with drug
s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
s2.chainA <- trim.pdb(s2, chain="A", elety="CA")
s3.chainA <- trim.pdb(s3, chain="A", elety="CA")
s1.b <- s1.chainA$atom$b
s2.b <- s2.chainA$atom$b
s3.b <- s3.chainA$atom$b
plotb3(s1.b, sse=s1.chainA, typ="l", ylab="Bfactor")
plotb3(s2.b, sse=s2.chainA, typ="l", ylab="Bfactor")
plotb3(s3.b, sse=s3.chainA, typ="l", ylab="Bfactor")

#WORKSHEET: SIMPLIFY CODE
##SPECIFY TO ONE EXAMPLE
x = "4AKE"

##SIMPLIFY ROW 3-5 = READ FILE
s <- read.pdb(x)

##SIMPLIFY ROW 6-8 = TRIM FILE TO SPECIFIC CHAIN/REGION/RESIDUE
s.chainA <- trim.pdb(s, chain ="A", elety= "CA")

##SIMPLIFY ROW 9-11 = ASSIGN NAME TO BFACTOR OF RESIDUE
s.b <- s.chainA$atom$b

##SIMPLIFY ROW 12-14 = PLOT RESIDUE BFACTOR
plotb3(s.b, sse = s.chainA, typ= "l", ylab = "Bfactor")


#WRITE FUNCTION
myfunction <- function(x) {
  ##SIMPLIFIED CODES
  s <- read.pdb(x)
  s.chainA <- trim.pdb(s, chain ="A", elety= "CA")
  s.b <- s.chainA$atom$b
  plotb3(s.b, sse = s.chainA, typ= "l", ylab = "Bfactor")
  ##SHOW PLOT
  return(plotb3)
}

#TESTING CODE
library(bio3d)
x= "4AKE"
myfunction(x)