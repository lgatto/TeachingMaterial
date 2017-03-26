Object-oriented programming in R
================================

Object-oriented programming (OOP) can help us encapsulate complex data
structures into a simple interface.

Unlike other main-stream languages, R has multiple OOP systems that
co-exist in parallel:

-   S3 - used in base R
-   S4 - *de-facto* standard in Bioconductor
-   Reference classes (RC) - used in special use cases

S3 and S4 are similar, but S4 is more formal. RC is quite different and
is Java-like.

OOP concepts
------------

-   Abstraction - related data is stored and handled together
-   Polymorphism - the most appropriate function is called based on the
    object type (e.g various plot functions)
-   Inheritance - code reuse by hierarchy of more-to-less general object
    types

In R, everything has a class (i.e type):

    class("Hello")

    ## [1] "character"

    class(1:5)

    ## [1] "integer"

    class(2.1)

    ## [1] "numeric"

    class(table(c(1,4,5,3,2,1)))

    ## [1] "table"

    class(plot)

    ## [1] "function"

    x = 10
    class(x)

    ## [1] "numeric"

Sections
========

-   [S3 OOP](s3.md)
-   [S4 OOP](s4.md)
-   [RC OOP](rc.md)
