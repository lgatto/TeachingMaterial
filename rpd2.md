# R Package Development

## Introduction

Packages are the way to share R code in a structures, reproducible and
tractable way. Even if the intend is not to disseminate your code,
packaging it is worth it. Packages provide a mechanism for loading
optional code and attached documentation as needed.

- logically group your own functions 
- keep code and documentation together and consistent
- keep code and data together
- keep track of changes in code 
- summarise all packages used for a analysis (see `sessionInfo()`)
- make a reproducible research compendium (container for code, text,
  data as a means for distributing, managing and updating)
- optionally test your code
- ... project managment

### References

- [R packages](http://r-pkgs.had.co.nz/), by Hadley Wickham
- R Installation and Administration [R-admin], R Core team 
- Writing R Extensions [R-ext], R Core team 

> Use `help.start()` to access them from your local installation, or
> http://cran.r-project.org/manuals.html from the web.

### Terminology

A **package** is loaded from a **library** by the function `library()`. 
Thus a library is a directory containing installed packages.

> Calling `library("foo", lib.loc = "/path/to/bar")` loads the package
> (book) `foo` from the library `bar` located at `/path/to/bar`.

### Requirement


```r
library("devtools")
library("roxygen2")
```


### Course content

Basic workflow

1. Prepare R code
2. Create package directory: `mypackage`
3. **Build** the package tarball
4. **Check** the package
5. **Install** the package

Step 2 is done only once. Package developement cycles through 3 - 5.

Also 
- Wrinting package documentation
- Vignettes
- Testing packages
- Compiled code

## Preparing R code


```r
fn <- function() 
    message("I love R packages")    
```

## Package layout

We can use

- `package.skeleton("myRpackage", list = "fn")`
- `devtools::create("myRpackage")` also create an `.Rproj` file.
- Use the RStudio wizard: `New Project > New Directory > R Package`

```
myRpackage/
|-- DESCRIPTION
|-- NAMESPACE
|-- man
|   `-- fun.Rd
`-- R
    `-- fun.R
```

This is the *source package*. From this, we need to create the
*package tarball* (or *package bundle*), i.e. a compressed archive of
the source. We can also create *binary packages* for Windows and Mac.

## Package developement cycle

In the shell

```
R CMD build myPackage               ## creates myRpackage_1.0.tar.gz
R CMD check myPackage_1.0.tar.gz    ## create myRpackage.Rcheck
R CMD INSTALL myRpackage_1.0.tar.gz ## Installation in the default library
```

Using RStudio useful keyboard shortcuts for package authoring:

* Build and Reload Package: `Ctrl + Shift + B`
* Check Package: `Ctrl + Shift + E`
* Test Package: `Ctrl + Shift + T`

Using devtools:
* `devtools::build()` 
* `devtools::build(binary = TRUE)`
* `devtools::check()`
* `devtools::install()`

A shortcut when developing:

* `devtools::load_all()`

## Package metadata

The `DESCRIPTION` file

```
Package: myRpackage ## mandatory (*)
Type: Package ## optional, 'Package' is default type 
Title: What the package does (short line) ## *
Version: 1.0 ##  *
Date: 2013-05-10 ## release date of the current version 
Author: Who wrote it ## *
Maintainer: Who to complain to <yourfault@somewhere.net> ## *
Description: More about what it does (maybe more than one line) ## *
License: What license is it under? ## *
Depends: methods, Biostrings ## for e.g.
Imports: evd ## for e.g.
Suggests: BSgenome.Hsapiens.UCSC.hg19 ## for e.g.
Collate: 'DataClasses.R' 'read.R' ## for e.g.
```

Package dependencies:

* **Depends** A comma-separated list of package names (optionally with
  versions) which this package depends on.
* **Suggests** Packages that are not necessarily needed: used only in
  examples, tests or vignettes, loaded in the body of functions
* **Imports** Packages whose name spaces are imported from (as
  specified in the `NAMESPACE` file) which do not need to be attached
  to the search path.
* **Collate** Controls the collation order for the R code files in a
  package. If filed is present, all source files must be listed.

 Packages are attached to the search path with \Rfunction{library} or \Rfunction{require}. 
 

* **Attach** When a package is attached, then all of its dependencies
  (see `Depends` field in its `DESCRIPTION` file) are also
  attached. Such packages are part of the evaluation environment and
  will be searched.
	  
* **Load** One can also use the `Imports` field in the `NAMESPACE`
  file. Imported packages are loaded but are not attached: they do not
  appear on the search path and are available only to the package that
  imported them.

## `NAMESPACE`

Restricts the symbols that are exported and imports functionality from
other packages.  Only the exported symbols will have to be documented.

```
export(f, g) ## exports f and g 
exportPattern("^[^\\.]")
import(foo) ## imports all symbols from package foo
importFrom(foo, f, g) ## imports f and g from foo
```

It is possible to explicitely use symbol `s` from package `foo` with
`foo::s` or `foo:::s` if `s` is not exported.

## R code

Contains `source()`able R source code to be installed. Files must
start with an ASCII (lower or upper case) letter or digit and have one
of the extensions `.R`, `.S`, `.q`, `.r`, or `.s` (use `.R` or
`.r`). 

* General style guidelines and best practice apply.
* Any number of files in `R`.
* Any number of functions (methods, classes) in each source file.
* Order matters (somehow), as the files will be sourced in the
  alphanumeric order. If that doesn't fit, use the `collate` field in
  the `DESCRIPTION` files.

Example

```
## works fine without Collate field
AllGenerics.R       DataClasses.R
methods-ClassA.R    methods-ClassB.R
functions-ClassA.R  ...
```

`zzz.R` is generally used to define special functions used to
initialize (called after a package is loaded and attached) and clean
up (just before the package is detached).  See `help(".onLoad"))`,
`?.First.Lib` and `?.Last.Lib` for more details.

## Package sub-directories

* **vignettes** directory for vignettes in Sweave or R markdown format.
* **data** for R code, compressed tables (`.tab`, `.txt`, or .`csv`,
  see `?data` for the file formats) and binary R objects. Available
  with `data()`.
* **inst/docs** for additional documentation. That's also where the
  vignettes will be installed after compilation.
* **inst/extdata** directory for other data files, not belonging in `data`.
* **tests** code for unit tests (see
  [here](https://github.com/lgatto/R-debugging) and
  [here](https://github.com/lgatto/2016-02-25-adv-programming-EMBL)).
* **src** for compiled code (see the
  [rccpp](https://github.com/lgatto/rccpp/) material)
* **demo** for demo code (see `?demo`)

## Documentation

## Manual pages

Package functions, datasets, methods and classes are documented in
`Rd`, a LaTeX-like format. 

```
% File src/library/base/man/load.Rd
% Part of the R package, http://www.R-project.org
% Copyright 1995-2014 R Core Team
% Distributed under GPL 2 or later

\name{load}
\alias{load}
\title{Reload Saved Datasets}
\description{
  Reload datasets written with the function \code{save}.
}
\usage{
load(file, envir = parent.frame(), verbose = FALSE)
}
\arguments{
  \item{file}{a (readable binary-mode) \link{connection} or a character string
    giving the name of the file to load (when \link{tilde expansion}
    is done).}
  \item{envir}{the environment where the data should be loaded.}
  \item{verbose}{should item names be printed during loading?}
}
\details{
  \code{load} can load \R objects saved in the current or any earlier
  format.  It can read a compressed file (see \code{\link{save}})
  directly from a file or from a suitable connection (including a call
  to \code{\link{url}}).

[...]

\value{
  A character vector of the names of objects created, invisibly.
}
\section{Warning}{
  Saved \R objects are binary files, even those saved with
  \code{ascii = TRUE}, so ensure that they are transferred without
  conversion of end of line markers.  \code{load} tries to detect such a
  conversion and gives an informative error message.

[...]

\examples{
## save all data
xx <- pi # to ensure there is some data
save(list = ls(all = TRUE), file= "all.RData")
rm(xx)

## restore the saved values to the current environment
local({
   load("all.RData")
   ls()
})

xx <- exp(1:3)
## restore the saved values to the user's workspace
load("all.RData") ## which is here *equivalent* to
## load("all.RData", .GlobalEnv)
## This however annihilates all objects in .GlobalEnv with the same names !
xx # no longer exp(1:3)
rm(xx)
attach("all.RData") # safer and will warn about masked objects w/ same name in .GlobalEnv
ls(pos = 2)
##  also typically need to cleanup the search path:
detach("file:all.RData")

## clean up (the example):
unlink("all.RData")

\dontrun{
con <- url("http://some.where.net/R/data/example.rda")
## print the value to see what objects were created.
print(load(con))
close(con) # url() always opens the connection
}}
\keyword{file}
```

These R documentation files can then be converted into text, pdf or
html:


```r
help("load")
help("load", help_type = "html")
help("load", help_type = "pdf")
```

One can use `prompt`, `promptClass`, `promptMethods`, `promptPackage`,
`promptPackage`, `promptData` to generate `Rd` templates for
functions, classes, methods, packages and data.

### Use roxygen

The way documentation is managed in R packages separates the code from
the documentation, which makes it easier to adapt the latter when to
code is updated. In comes `roxygen2`, that allows developer to write
their documentation on top of their functions:

```
#' Reads sequences data in fasta and create \code{DnaSeq}
#' and \code{RnaSeq} instances. 
#'
#' This funtion reads DNA and RNA fasta files and generates
#' valid \code{"DnaSeq"} and \code{"RnaSeq"} instances.
#' 
#' @title Read fasta files.
#' @param infile  the name of the fasta file which the data are to be read from.
#' @return an instance of \code{DnaSeq} or \code{RnaSeq}.
#' @seealso \code{\linkS4class{GenericSeq}}, \code{\linkS4class{DnaSeq}} and \code{\linkS4class{RnaSeq}}.
#' @examples
#' f <- dir(system.file("extdata",package="sequences"),pattern="fasta",full.names=TRUE)
#' f
#' aa <- readFasta(f)
#' aa
#' @author Laurent Gatto \email{lg390@@cam.ac.uk}
#' @export
readFasta <- function(infile){
  lines <- readLines(infile)
  header <- grep("^>", lines)
  if (length(header)>1) {
    warning("Reading first sequence only.")
    lines <- lines[header[1]:(header[2]-1)]
    header <- header[1]
  }
  ##### (code cut for space reasons) #####
  if (validObject(newseq))
    return(newseq)
}
```

The roxygen code can then be parsed and converted to `Rd` using
`roxygen2::roxygenise` function. 

Note that the `roxygenise` function does more than produce
documentation (that part is handled by the `rd` roclet, set with
`roclete = "rd"`). It can also manage your `NAMESPACE` file and
`Collate` field.

Note: recently, support for markdown format has been added to roxygen.

### Vignettes

While manual pages are meant to be specific and technical, vignettes
are workflow-type documentation files that provide an overview and/or
a use-case demonstrating the package's functionality.

Vignettes can be written in Sweave format (`.Rnw` extension),
supporting R code chunks in LaTeX documents, or R markdown formart
(`.Rmd` extension) for R code and markdown.

The source document (`.Rnw` or `.Rmd`) can be 

* weaved into `tex` or `md` files respectively
* and converted into `pdf` or `html` (`Rnw` to `pdf` only, `Rmd` to
  either)

using the `utils::Sweave` (`Rnw` only) or `knitr::knit` functions.

Rstudio makes it very easy to write and compile `Rmd` documents
(independently of R packages). When inside a package, the documents
are stored in the `vignettes` directory and compiled/converted
automatically when the package is built.

Note: if you use `knitr` and `rmarkdown` for your vignette, you'll
have to add these dependencies in the `Suggests` field and specify
`VignetteBuilder: knitr` in the `DESCRIPTION` file.

## Additional files

* `.Rbuildignore` with a list of files/dirs to ignore when
  building. For example the `.Rproj` file.
* `.Rinstignore`with a list of files/dirs to ignore when installing.
* `CITATION` file (see `citation()` function)
* `README.Rmd`/`README.md` files if you use github.

## Distributing packages

* **CRAN** Read the CRAN Repository Policy
  (http://cran.r-project.org/web/packages/policies.html). Upload your
  `--as-cran` checked} `myPackage\_x.y.z.tar.gz` to
  `ftp://cran.R-project.org/incoming` or using
  `http://CRAN.R-project.org/submit.html`. Your package will be
  installable with `install.packages("myRpackage")`.
		
* **R-forge** Log in, register a project and wait for acceptance. Then
  commit you code to the svn repository. Your package will be
  installable with `install.packages` using
  `repos="http://R-Forge.R-project.org"`.


* **GitHub** (and **Bitbucket**) Great for development and promoting
  interaction and contributions. Unofficial. Autmatic checking
  possible through CI such as `travis-ci` for example. Packages can be
  installed with `devtools::install_github`
  (`devtools::install_bitbucket`).

* **Bioconductor** Make sure to satisfy submission criteria (pass
  `check` (and `BiocCheck`), have a vignette, use S4 if OO, make use
  of appropriate existing infrastructure, include a NEWS file, must
  **not** already be on CRAN, ...). Your package will then be reviewed
  on github publicly before acceptance. A svn (git very soon) account
  will then be created. Package will be installable with
  `biocLite("myPackage")`.


