Reference classes OOP
=====================

Reference classes are a Java-like class system (no generics!)
implementation in R.

Reference classes use pass-by-reference logic that is different from the
usual pass-by-value logic R uses. This means when we assign `a=b` both
`a` and `b` point to the same object, and changing `a` will also modify
`b` and vica versa.

Creating classes
----------------

Similar syntax as S4, but now slots are called fields, and `setRefClass`
returns a constructor object.

    Seq <- setRefClass("Seq",
                       fields = list(
                         name = "character",
                         alphabet = "character",
                         sequence = "character"))

We then use the constuctor to create objects:

    a = Seq(name="My seq", alphabet=c("A", "T"), sequence="ATTTTAAAAAA")
    a

    ## Reference class object of class "Seq"
    ## Field "name":
    ## [1] "My seq"
    ## Field "alphabet":
    ## [1] "A" "T"
    ## Field "sequence":
    ## [1] "ATTTTAAAAAA"

Defining methods
----------------

Methods are defined either with the class, or later by modifying the
constructor object.

    Seq <- setRefClass("Seq",
                       fields = list(
                         name = "character",
                         alphabet = "character",
                         sequence = "character"),
                       methods = list(
                         comp = function() {
                           "Complements the (DNA) sequence" ## inline docs
                           sequence <<- chartr("ACGT","TGCA",.self$sequence)
                           name <<- paste(.self$name, "-- complemented")                     
                         }
                         ## there can be more, of course
                         ))

A couple of things to note here:

-   accessing fields and calling methods is done with the `$` operator.
-   the current object can be referred to in a method by the reserved
    name `.self`.
-   Changing fields of an object within methods needs to be done with
    the `<<-` operator.

Calling a method on the object will modify it:

    s <- Seq$new(name="foo", sequence="GATCATCA")
    s

    ## Reference class object of class "Seq"
    ## Field "name":
    ## [1] "foo"
    ## Field "alphabet":
    ## character(0)
    ## Field "sequence":
    ## [1] "GATCATCA"

    s$sequence

    ## [1] "GATCATCA"

    s$comp()
    s$sequence

    ## [1] "CTAGTAGT"

Objects are copied by reference and not value:

    s2 = s
    s$sequence

    ## [1] "CTAGTAGT"

    s2$sequence

    ## [1] "CTAGTAGT"

    s$sequence = "ATTTA"
    s$sequence

    ## [1] "ATTTA"

    s2$sequence

    ## [1] "ATTTA"

Summary
=======

Reference classes are suitable for objects that are dynamically tracked
by all the code: GUI components, read-only access to files (streams,
data bases), internet resources, editing facilities, ... Because they
are passed by reference they are *never* copied, which means they are
never duplicated in memory.

Because of their pass-by-reference semantics they are not appropriate as
a replacement for S3 and S4 as this would lead to a lot of unanticipated
side effects.
