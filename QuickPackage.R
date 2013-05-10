### R code from vignette source 'QuickPackage.Rnw'

###################################################
### code chunk number 1: knitr
###################################################
library("knitr")
options(width = 60)
opts_chunk$set(prompt = TRUE,
               comment = '',
               fig.align = 'center')


###################################################
### code chunk number 2: dummy
###################################################
## This is the dummy function that will return
## information about our 'QuickPackage'
qpf <- function() 
  packageDescription("QuickPackage")


###################################################
### code chunk number 3: pgkskel
###################################################
## This creates the package template
package.skeleton("QuickPackage", list = c("qpf"))


###################################################
### code chunk number 4: clean
###################################################
rm(list=ls())
system("rm ./QuickPackage_1.0.tar.gz")
system("rm ./QuickPackage/Read-and-delete-me")
system("rm ./QuickPackage/man/QuickPackage-package.Rd")
system("cp files/qpf.Rd QuickPackage/man/.")
system("cp files/DESCRIPTION QuickPackage/.")
system("R CMD build QuickPackage")
system("R CMD check QuickPackage_1.0.tar.gz")
system("R CMD INSTALL QuickPackage_1.0.tar.gz")


###################################################
### code chunk number 5: enjoy
###################################################
library("QuickPackage")
qpf()


