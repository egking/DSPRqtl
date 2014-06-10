## ----setup, include=FALSE------------------------------------------------
library("knitr")

## ----eval=FALSE----------------------------------------------------------
#  install.packages("DSPRqtl",
#                   repos = "http://wfitch.bio.uci.edu/R/",
#                   type = "source")

## ----eval=FALSE----------------------------------------------------------
#  install.packages("DSPRqtlDataA",
#                   repos = "http://wfitch.bio.uci.edu/R/",
#                   type = "source")

## ----eval=FALSE----------------------------------------------------------
#  install.packages("DSPRqtlDataB",
#                   repos = "http://wfitch.bio.uci.edu/R/",
#                   type = "source")

## ----eval=FALSE----------------------------------------------------------
#  install.packages("DSPRqtlDataA_2.0-1.tar.gz", repos = NULL, type = "source")
#  install.packages("DSPRqtlDataB_2.0-1.tar.gz", repos = NULL, type = "source")

## ----eval=FALSE----------------------------------------------------------
#  library(DSPRqtlDataA)

## ----eval=FALSE----------------------------------------------------------
#  library(DSPRqtlDataB)

## ----eval=FALSE----------------------------------------------------------
#  data(A_chromosome.arm_position.in.base.pairs)

## ----eval=FALSE----------------------------------------------------------
#  data(A_X_12000000)
#  # This gives a data.frame named A_X_12000000

## ------------------------------------------------------------------------
library(DSPRqtl)
data(positionlist_wgenetic)
# This gives a data.frame named poslist

## ------------------------------------------------------------------------
library(DSPRqtl)

## ------------------------------------------------------------------------
data(ADHdata)
# a data.frame named ADHdata

head(ADHdata)

## ----eval=FALSE----------------------------------------------------------
#  phenotype.dat <- read.table("/path/to/your/file.txt",
#                              header = TRUE)

## ----eval=FALSE----------------------------------------------------------
#  DSPRscan(model, design, phenotype.dat, id.col, batch, sex)

## ----eval=FALSE----------------------------------------------------------
#  data(ADHdata)
#  scan.results <- DSPRscan(adh ~ 1,
#                           design = "inbredA",
#                           phenotype.dat = ADHdata,
#                           id.col='patRIL')

## ------------------------------------------------------------------------
data(ADHscan)

## ------------------------------------------------------------------------
ADH.lod.scores <- ADHscan$LODscores

ADH.lod.scores[100:105, ]

## ----eval=FALSE----------------------------------------------------------
#  DSPRperm(model, design, phenotype.dat, id.col, batch, niter, alpha, sex)

## ----eval=FALSE----------------------------------------------------------
#  perm.test <- DSPRperm(adh ~ 1,
#                        design = "inbredA",
#                        phenotype.dat = ADHdata)

## ----eval=FALSE----------------------------------------------------------
#  quantile(perm.test$maxLODs, 1 - 0.01)

## ----eval=FALSE----------------------------------------------------------
#  DSPRpeaks(qtldat, method, threshold, LODdrop, BCIprob)

## ----echo=FALSE----------------------------------------------------------
load("peaks.rda")

## ----slowchunk, eval=FALSE-----------------------------------------------
#  peaks <- DSPRpeaks(ADHscan, threshold = 6.8, LODdrop = 2)

## ------------------------------------------------------------------------
peaks[[26]]

## ----dev="png", fig.width=6.5, fig.height=3.5, fig.cap="Output of `DSPRplot()`"----
DSPRplot(list(ADHscan), threshold=6.8)

## ----eval=FALSE----------------------------------------------------------
#  LocalInt(peakChr, peakPos, range, phenotype.dat, pheno.name, design)

## ------------------------------------------------------------------------
# The main QTL
main.peak <- peaks[[26]]
peakChr <- main.peak$peak$chr
peakPos <- main.peak$peak$Ppos

## ----eval=FALSE----------------------------------------------------------
#  peak.int <- LocalInt(peakChr,
#                       peakPos,
#                       phenotype.dat = ADHdata,
#                       pheno.name = "adh",
#                       design = "inbredA")

## ----eval=FALSE----------------------------------------------------------
#  DSPRgenos(design, phenotype.dat, id.col, output)

## ----eval=FALSE----------------------------------------------------------
#  myGenos <- DSPRgenos(design, phenotype.dat, id.col, output)
#  save(myGenos, file = "my file path")

## ----eval=FALSE----------------------------------------------------------
#  load(file="my file path")
#  # The object will load into the workspace with the same name as
#  # when save was used. In this example, this is myGenos.

