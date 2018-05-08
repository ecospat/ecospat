## ----load_library--------------------------------------------------------
library(ecospat)
citation("ecospat")

## ------------------------------------------------------------------------
data(ecospat.testData)
names(ecospat.testData)

## ------------------------------------------------------------------------
data(ecospat.testNiche.inv)
names(ecospat.testNiche.inv)

## ------------------------------------------------------------------------
data(ecospat.testNiche.nat)
names(ecospat.testNiche.nat)

## ------------------------------------------------------------------------
fpath <- system.file("extdata", "ecospat.testTree.tre", package="ecospat")
fpath
tree<-read.tree(fpath)
tree$tip.label

## ----tree----------------------------------------------------------------
plot(tree, cex=0.6)

## ----mantel_cor----------------------------------------------------------
ecospat.mantel.correlogram(dfvar=ecospat.testData[c(2:16)],colxy=1:2, n=100, 
                           colvar=3:7, max=1000, nclass=10, nperm=100)

## ------------------------------------------------------------------------
colvar <- ecospat.testData[c(4:8)]
x <- cor(colvar, method="pearson")
ecospat.npred (x, th=0.75)

## ------------------------------------------------------------------------
x <- cor(colvar, method="spearman")
ecospat.npred (x, th=0.75)

## ------------------------------------------------------------------------
x <- ecospat.testData[c(4:8)]
p<- x[1:90,] #A projection dataset.
ref<- x[91:300,] # A reference dataset

## ------------------------------------------------------------------------
ecospat.climan(ref,p)

## ------------------------------------------------------------------------
x <- ecospat.testData[c(2,3,4:8)]
proj<- x[1:90,] #A projection dataset.
cal<- x[91:300,] #A calibration dataset

## ------------------------------------------------------------------------
mess.object<-ecospat.mess (proj, cal, w="default")

