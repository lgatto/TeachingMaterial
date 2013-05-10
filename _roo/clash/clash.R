
# make package1
system("rm -rf package1")
package.skeleton("package1", code_files="files1.R")
system("rm -rf package1/man/*")
system("R CMD build package1")
install.packages("package1_1.0.tar.gz")


# now packge 2
system("rm -rf package2")
package.skeleton("package2", code_files="files2.R")
system("rm -rf package2/man/*")
system("R CMD build package2")
install.packages("package2_1.0.tar.gz")


cat("=====================================================\n")
library(package1)
library(package2)
