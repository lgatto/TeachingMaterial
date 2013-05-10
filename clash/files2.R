setClass("Package2Class", 
	representation=representation(
		param = "numeric"
	)
)

setGeneric("score", function(x, ...) standardGeneric("score"))

setMethod("score", "Package2Class", function(x, ...){ -x@param })

c2 = new("Package2Class", param=10)
print(score(c2))
