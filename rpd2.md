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

Package developement cycles through 3 - 5.

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


```r
package.skeleton("myRpackage", list = "fn")
```

```
## Creating directories ...
```

```
## Creating DESCRIPTION ...
```

```
## Creating NAMESPACE ...
```

```
## Creating Read-and-delete-me ...
```

```
## Saving functions and data ...
```

```
## Making help files ...
```

```
## Done.
```

```
## Further steps are described in './myRpackage/Read-and-delete-me'.
```

```
myRpackage/
|-- DESCRIPTION
|-- NAMESPACE
|-- man
|   `-- myRpackage-package.Rd
|   `-- fun.Rd
|-- R
|   `-- fun.R
`-- Read-and-delete-me

2 directories, 6 files
```

We can use

- `package.skeleton` (only once!)
- 




