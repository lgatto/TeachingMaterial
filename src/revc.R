system("R CMD SHLIB rev.c")

obj <- structure(seq(0,1,0.1), myattr="important!")

dyn.load("rev.so")

.Call("rev", obj)
