#' ---
#' title: "class07"
#' output: github_document
#' ---

#2A: line plot
# data input of weight of baby
weight <- read.table("bimm143_05_rstats/bimm143_05_rstats/weight_chart.txt", header = TRUE)
# data plot of weight of baby
plot(weight, typ="o", pch=15, cex=1.0, lwd=2, ylim=c(2,10), xlab="Age (months)", ylab="weight (kg)", main= "Growth of baby")

#2B: bar plot
#data input of feature counts
feature <- read.table("bimm143_05_rstats/bimm143_05_rstats/feature_counts.txt", sep="\t", header=TRUE)
#plot features
par(mar=c(3.1,11.0,4,2))  
barplot(feature$Count, horiz = T, ylab = "", names.arg = feature$Feature, main = "Freq of genetic features", las=1, xlim = c(0,80000))
   
#2C: histograms
hist( c(rnorm(10000), rnorm(10000)+4), breaks = 50)

#3A: color vectors
MF_count <- read.table("bimm143_05_rstats/bimm143_05_rstats/male_female_counts.txt", header=T, sep = "\t")
barplot(MF_count$Count, ylab = "Counts", names.arg = MF_count$Sample, las= 2, col = rainbow(nrow(MF_count)))

#3B: color values
genes <- read.delim("bimm143_05_rstats/bimm143_05_rstats/up_down_expression.txt")
par(mar=c(5,4,4,2))
plot(genes$Condition1, genes$Condition2, col= genes$State, xlab = "expression condition 1", ylab = "expression condition 2")

#3C: color dynamics
meth <- read.delim("bimm143_05_rstats/bimm143_05_rstats/expression_methylation.txt")
mycols <- densCols(meth$gene.meth, meth$expression)
plot(meth$gene.meth, meth$expression, col=mycols)

#narrow data
sig <- meth$expression > 0
mycols2 <- densCols(meth$gene.meth[sig], meth$expression[sig])
plot(meth$gene.meth[sig], meth$expression[sig], col=mycols2)

rescale <- function(x, na.rm=TRUE, plot=FALSE) {
  if(na.rm) {
    rng <-range(x, na.rm=TRUE)
  } else {
    rng <-range(x)
  }
  print("Hello")
  answer <- (x - rng[1]) / (rng[2] - rng[1])
  return(answer)
  print("is it me you are looking for?")
  if(plot) {
    plot(answer, typ="b", lwd=4)
  }
  print("I can see it in ...")
}

rescale(1:10)
