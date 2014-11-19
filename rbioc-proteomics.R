
## ----docver, echo=FALSE--------------------------------------------------
v <- system('git log -1 --format="%h [%ci]"', intern = TRUE)


## ----env0, message=FALSE, echo=FALSE, warning=FALSE----------------------
library("knitr")
opts_knit$set(error = FALSE)
library("BiocInstaller")
library("RforProteomics")


## ----, env, message=FALSE, echo=TRUE, warning=FALSE----------------------
library("mzR")
library("mzID")
library("MSnID")
library("MSGFplus")
library("MSnbase")
library("rpx")
library("MLInterfaces")
library("pRoloc")
library("pRolocdata")
library("rTANDEM")
library("MSGFplus")
library("MSGFgui")
library("rols")
library("hpar")


## ----r4pinstall, eval=FALSE----------------------------------------------
## library("BiocInstaller")
## biocLite("RforProteomics", dependencies = TRUE)


## ----, pk, echo=FALSE, warning=FALSE, cache=TRUE-------------------------
biocv <- as.character(biocVersion())
pp <- proteomicsPackages(biocv)
msp <- massSpectrometryPackages(biocv)
msdp <- massSpectrometryDataPackages(biocv)


## ----, pp, eval=FALSE----------------------------------------------------
## library("RforProteomics")
## pp <- proteomicsPackages()
## display(pp)


## ----, datatab, results='asis', echo=FALSE-------------------------------

datatab <-
    data.frame(Type = c("raw", "identification", "quantitation",
                   "peak lists", "other"),
               Format = c("mzML, mzXML, netCDF, mzData",
                   "mzIdentML", "mzQuantML", "mgf", "mzTab"),
               Package = c(
                   "[`mzR`](http://bioconductor.org/packages/release/bioc/html/mzR.html) (read)",
                   "`mzR` and [`mzID`](http://bioconductor.org/packages/release/bioc/html/mzID.html) (read)",
                   "",
                   "[`MSnbase`](http://bioconductor.org/packages/release/bioc/html/MSnbase.html) (read/write)", 
                   "[`MSnbase`](http://bioconductor.org/packages/release/bioc/html/MSnbase.html) (read/write)"))
library("knitr")
kable(datatab)


## ----, rpx---------------------------------------------------------------
library("rpx")
pxannounced()


## ----, pxd, cache=TRUE---------------------------------------------------
px <- PXDataset("PXD000001")
px
pxfiles(px)


## ----, pxvar, eval=FALSE-------------------------------------------------
## pxtax(px)
## pxurl(px)
## pxref(px)


## ----, pxget-------------------------------------------------------------
## mzf <- pxget(px, pxfiles(px)[6])
mzf <- "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML"
mzf


## ----pxd000561, cache=TRUE, eval=FALSE, echo=FALSE-----------------------
## library("rpx")
## hum <- PXDataset("PXD000561")
## hum
## humf <- pxfiles(hum)
## length(humf)
## table(sub("^.+\\.", "", humf))
## rawf <- grep("raw", humf, value = TRUE)
## table(sub("_.+$", "", rawf))


## ----, rawms-------------------------------------------------------------
library("mzR")
ms <- openMSfile(mzf)
ms


## ----, hd----------------------------------------------------------------
hd <- header(ms)
dim(hd)
names(hd)


## ----, ex_raw, fig.align='center'----------------------------------------
hd2 <- hd[hd$msLevel == 2, ]
i <- which.max(hd2$basePeakIntensity)
hd2[i, ]
pi <- peaks(ms, hd2[i, 1])
mz <- hd2[i, "basePeakMZ"]

par(mfrow = c(2, 2))
plot(pi, type = "h", main = paste("Acquisition", i))
plot(pi, type = "h", xlim = c(mz-0.5, mz+0.5))

pj <- peaks(ms, 100)
plot(pj, type = "l", main = paste("Acquisition", 100))
plot(pj, type = "l", xlim = c(536,540))


## ----, id, cache=TRUE----------------------------------------------------
library("mzID")
f <- dir(system.file("extdata", package = "RforProteomics"),
         pattern = "mzid", full.names=TRUE)
basename(f)
id <- mzID(f)
id


## ----, ex_id-------------------------------------------------------------
fid <- flatten(id)
x <- by(fid, fid$accession, function(x)
    c(unique(x$length),
      length(unique(x$pepseq)),
      mean(x$'ms-gf:specevalue')))
x <- data.frame(do.call(rbind, x))
colnames(x) <- c("plength", "npep", "eval")
x$bins <- cut(x$eval, summary(x$eval))
library("lattice")
xyplot(plength ~ npep | bins, data = x)


## ----mzrvsid, eval = FALSE-----------------------------------------------
## library("mzR")
## library("mzID")
## f <- dir(system.file("extdata", package = "RforProteomics"),
##          pattern = "mzid", full.names=TRUE)
## 
## system.time({
##     id0 <- mzID(f)
##     fid0 <- flatten(id0)
## })
## 
## head(fid0)
## 
## system.time({
##     id1 <- openIDfile(f)
##     fid1 <- psms(id1)
## })
## 
## head(fid1)


## ----, rtandem, eval=FALSE-----------------------------------------------
## library("rTANDEM")
## ?rtandem
## library("shinyTANDEM")
## ?shinyTANDEM


## ----msgfplus, eval=FALSE------------------------------------------------
## library("MSGFplus")
## parameters <- msgfPar(database = 'proteins.fasta',
##                       tolerance='20 ppm',
##                       instrument='TOF',
##                       enzyme='Lys-C')
## runMSGF(parameters, c('file1.mzML', 'file2.mzML'))


## ----msgfgui, eval=FALSE-------------------------------------------------
## library("MSGFgui")
## MSGFgui()


## ----, echo=FALSE--------------------------------------------------------
mzf <- "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML"


## ----ex_getfas-----------------------------------------------------------
## fas <- pxget(px, pxfiles(px)[8])
fas <- "erwinia_carotovora.fasta"


## ----ex_msgfcmd----------------------------------------------------------
msgf <- system.file(package = "MSGFplus", "MSGFPlus", "MSGFPlus.jar")
system(paste0("java -jar ", msgf))
cmd <- paste("java -jar", msgf, "-protocol 2 -inst 1 -s", mzf, "-d", fas)
cmd


## ----ex_msgfsys, eval=FALSE----------------------------------------------
## system(cmd)


## ----ex_msgfplus, eval=FALSE---------------------------------------------
## library("MSGFplus")
## msgfpar <- msgfPar(database = fas,
##                instrument = 'HighRes',
##                enzyme = 'Trypsin'm
##                protocol = 'iTRAQ')
## runMSGF(msgfpar, mzf)


## ----ex_msgfgui, eval=FALSE----------------------------------------------
## library("MSGFgui")
## MSGFgui()


## ----, msnid-------------------------------------------------------------
library("MSnID")
msnid <- MSnID(".")

PSMresults <- read.delim(system.file("extdata", "human_brain.txt",
                                     package="MSnID"),
                         stringsAsFactors=FALSE)
psms(msnid) <- PSMresults
show(msnid)


## ----msnidfilt-----------------------------------------------------------
msnid$msmsScore <- -log10(msnid$`MS.GF.SpecEValue`)
msnid$absParentMassErrorPPM <- abs(mass_measurement_error(msnid))

filtObj <- MSnIDFilter(msnid)
filtObj$absParentMassErrorPPM <- list(comparison="<", threshold=5.0)
filtObj$msmsScore <- list(comparison=">", threshold=8.0)
show(filtObj)
filtObj.grid <- optimize_filter(filtObj, msnid, fdr.max=0.01,
                                method="Grid", level="peptide",
                                n.iter=500)
show(filtObj.grid)
msnid <- apply_filter(msnid, filtObj.grid)
show(msnid)


## ----, msnset, echo=FALSE, fig.width = 5, fig.height = 7, fig.align='center'----
plot(NA, xlim = c(0, 5), ylim = c(0, 10), axes=FALSE, xlab = NA, ylab = NA)
rect(0, 0, 3, 1.9)
rect(0, 2, 3, 10)
rect(3.05, 2, 5, 10)

segments(seq(0, 3, length.out = 7),
         rep(0, 7),
         seq(0, 3, length.out = 7),
         rep(10, 7),
         lty = "dotted")

segments(rep(0, 50),
         seq(2, 10, length.out = 50),
         rep(5, 100),
         seq(2, 10, length.out = 50),
         lty = "dotted")

text(1.5, 1, "sample metadata", cex = 1.5)
text(1.5, 6, "assay data", cex = 1.5)
text(4, 6, "feature\nmetadata", cex = 1.5)


## ----, msnbase-----------------------------------------------------------
library("MSnbase")
rawFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
               full.name = TRUE, pattern = "mzXML$")
basename(rawFile)
msexp <- readMSData(rawFile, verbose = FALSE)
msexp


## ------------------------------------------------------------------------
length(msexp)
msexp[[2]]


## ----, addid-------------------------------------------------------------
fData(msexp)
## find path to a mzIdentML file
identFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                 full.name = TRUE, pattern = "dummyiTRAQ.mzid")
basename(identFile)
msexp <- addIdentificationData(msexp, identFile)
fData(msexp)


## ----, specplot----------------------------------------------------------
msexp[[1]]
plot(msexp[[1]], full=TRUE)


## ----, specplot2---------------------------------------------------------
msexp[1:3]
plot(msexp[1:3], full=TRUE)


## ------------------------------------------------------------------------
as(msexp[[1]], "data.frame")[100:105, ]


## ----, quanttab, echo=FALSE, results='asis'------------------------------

qtb <- matrix(c("XIC", "Counting", "SILAC, 15N", "iTRAQ, TMT"),
              nrow = 2, ncol = 2)
dimnames(qtb) <- list(
    'MS level' = c("MS1", "MS2"),
    'Quantitation' = c("Label-free", "Labelled"))

kable(qtb)



## ----, itraq4plot--------------------------------------------------------
plot(msexp[[1]], full=TRUE, reporters = iTRAQ4)


## ----, quantitraq--------------------------------------------------------
msset <- quantify(msexp, method = "trap", reporters = iTRAQ4, verbose=FALSE)
exprs(msset)
processingData(msset)


## ----, lfms2-------------------------------------------------------------
exprs(si <- quantify(msexp, method = "SIn"))     
exprs(saf <- quantify(msexp, method = "NSAF"))


## ----, mztab-------------------------------------------------------------
## mztf <- pxget(px, pxfiles(px)[2])
mztf <- "F063721.dat-mztab.txt"
(mzt <- readMzTabData(mztf, what = "PEP"))


## ----, readmsnset2-------------------------------------------------------
csv <- dir(system.file ("extdata" , package = "pRolocdata"),
           full.names = TRUE, pattern = "pr800866n_si_004-rep1.csv")
getEcols(csv, split = ",")
ecols <- 7:10
res <- readMSnSet2(csv, ecols)
head(exprs(res))
head(fData(res))


## ----, pure--------------------------------------------------------------
data(itraqdata)
qnt <- quantify(itraqdata, method = "trap",
                reporters = iTRAQ4, verbose = FALSE)
impurities <- matrix(c(0.929,0.059,0.002,0.000,
                       0.020,0.923,0.056,0.001,
                       0.000,0.030,0.924,0.045,
                       0.000,0.001,0.040,0.923),
                     nrow=4, byrow = TRUE)
## or, using makeImpuritiesMatrix()
## impurities <- makeImpuritiesMatrix(4)
qnt.crct <- purityCorrect(qnt, impurities)
processingData(qnt.crct)


## ----, pureplot----------------------------------------------------------

plot0 <- function(x, y, main = "") {
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    par(mar = c(4, 4, 1, 1))
    par(mfrow = c(2, 2))
    sx <- sampleNames(x)
    sy <- sampleNames(y)
    for (i in seq_len(ncol(x))) {
        plot(exprs(x)[, i], exprs(y)[, i], log = "xy",
             xlab = sx[i], ylab = sy[i])
        grid()
    }
}

plot0(qnt, qnt.crct)


## ----, norm, fig.align='center'------------------------------------------
qnt.crct.nrm <- normalise(qnt.crct,"quantiles")
plot0(qnt, qnt.crct.nrm)


## ----, comb--------------------------------------------------------------
## arbitraty grouping
g <- factor(c(rep(1, 25), rep(2, 15), rep(3, 15)))
prt <- combineFeatures(qnt.crct.nrm, groupBy = g, fun = "sum")
processingData(prt)


## ----impute--------------------------------------------------------------
set.seed(1)
qnt0 <- qnt
exprs(qnt0)[sample(prod(dim(qnt0)), 10)] <- NA
table(is.na(qnt0))
qnt00 <- filterNA(qnt0)
dim(qnt00)
qnt.imp <- impute(qnt0)
plot0(qnt, qnt.imp)


## ----, msmstest----------------------------------------------------------
library(msmsTests)
data(msms.dataset)
msms.dataset
e <- pp.msms.data(msms.dataset)
e
     
null.f <- "y~batch"
alt.f <- "y~treat+batch"
div <- apply(exprs(e),2,sum)
res <- msms.edgeR(e,alt.f,null.f,div=div,fnm="treat")
     
head(res)


## ----, ml----------------------------------------------------------------
library("MLInterfaces")
library("pRoloc")
library("pRolocdata")
data(dunkley2006)
traininds <- which(fData(dunkley2006)$markers != "unknown")
ans <- MLearn(markers ~ ., data = t(dunkley2006), knnI(k = 5), traininds)
ans


## ----clust---------------------------------------------------------------
kcl <- MLearn( ~ ., data = dunkley2006, kmeansI, centers = 12)
kcl
plot(kcl, exprs(dunkley2006))


## ----clust2--------------------------------------------------------------
hcl <- MLearn( ~ ., data = t(dunkley2006), hclustI(distFun =  dist, cutParm = list(k = 4)))
hcl
plot(hcl, exprs(t(dunkley2006)))


## ----nont, echo=FALSE, cache=TRUE----------------------------------------
library("rols")
nont <- nrow(ontologies())


## ----rols----------------------------------------------------------------
library("rols")
olsQuery("ESI", "MS")


## ----, si, echo=FALSE----------------------------------------------------
print(sessionInfo(), local = FALSE)

