setClass("Package1Class", 
	representation=representation(
		param = "numeric"
	)
)

setGeneric("score", function(object, ...) standardGeneric("score"))

setMethod("score", signature=c("object"="Package1Class"), 
	function(object, ...){ object@param^2 })

c1 = new("Package1Class", param=10)
print(score(c1))
