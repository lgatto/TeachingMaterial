---
layout: page
title: minimal make
tagline: A minimal tutorial on make
---

I would argue that the most important tool for reproducible research
is not [Sweave](http://www.stat.uni-muenchen.de/~leisch/Sweave/) or 
[knitr](http://yihui.name/knitr/) but
*[GNU make](http://www.gnu.org/software/make)* (or
[rake](http://rake.rubyforge.org)).

Consider, for example, all of the files associated with a manuscript.
In the simplest case, I would have an [R](http://r-project.org)
script for each figure plus a [LaTeX](http://www.latex-project.org)
file for the main text.  And then a [BibTeX](http://www.bibtex.org)
file for the references.

Compiling the final PDF is a bit of work: 

- Run each R script through R to produce the relevant figure
- Run latex and then bibtex and then latex a couple of more times.

And the R scripts need to be run before latex is, and only if they've
changed.

### A simple example

[GNU make](http://www.gnu.org/software/make) makes this easy.  In your
directory for the manuscript, you create a text file called `Makefile`
that looks something like [the following](assets/Makefile) (here using
[pdflatex](http://www.tug.org/applications/pdftex/)).

    mypaper.pdf: mypaper.bib mypaper.tex Figs/fig1.pdf Figs/fig2.pdf
    	pdflatex mypaper
    	bibtex mypaper
    	pdflatex mypaper
    	pdflatex mypaper

    Figs/fig1.pdf: R/fig1.R
    	cd R;R CMD BATCH fig1.R fig1.Rout

    Figs/fig2.pdf: R/fig2.R
    	cd R;R CMD BATCH fig2.R fig2.Rout

Each batch of lines indicates a file to be created (the _target_), the files it
depends on (the _dependencies_), and then a set of commands needed to
construct the target from the dependent files.  Note that the lines
with the commands _must_ start with a tab character.  

One oddity: if you need to change directories to run a command, do
the `cd` on the same line as the related command.  The following
_*would not work*_:

    ### this doesn't work ###
    Figs/fig1.pdf: R/fig1.R
    	cd R
    	R CMD BATCH fig1.R fig1.Rout

Another great feature: in the example above, you'd only build
`fig1.pdf` when `fig1.R` changed.  And note that the dependencies
propagate.  If you change `fig1.R`, then `fig1.pdf` will change, and
so `mypaper.pdf` will be re-built.

### Using GNU make

You probably already have GNU make installed on your computer.  Type
`make --version` in a terminal/shell to see.

To use make:

- Go into the the directory for your project.
- Create the `Makefile` file.
- Every time you want to build the project, type `make`.
- In the example above, if you want to build `fig1.pdf` without
  building `mypaper.pdf`, just type `make fig1.pdf`.

### More complicated examples

There are complicated Makefiles all over the place.  Poke around
[github](http://github.com) and study them.

Here are some of my own examples:

- [Makefile](https://github.com/kbroman/ailProbPaper/blob/master/Makefile)
  for my [AIL probabilities paper](http://www.g3journal.org/content/2/2/199.long)

- [Makefile](https://github.com/kbroman/phyloQTLpaper/blob/master/Makefile)
  for
  [my phylo QTL paper](http://www.genetics.org/content/192/1/267.full)
  
- [Makefile](https://github.com/kbroman/preCCProbPaper/blob/master/Makefile)
  for my
  [pre-CC probabilities paper](http://www.genetics.org/content/190/2/403.full)

- [Makefile](https://github.com/kbroman/Talk_InteractiveGraphs1/blob/master/Makefile) 
  for a [talk on interactive graphs](http://www.biostat.wisc.edu/~kbroman/talks/InteractiveGraphs/).

- [Makefile](https://github.com/kbroman/Talk_FunQTL/blob/master/Makefile)
  for a [talk on QTL mapping for function-valued traits](http://www.biostat.wisc.edu/~kbroman/talks/FunQTL/).

And here are some examples from [Mike Bostock](http://bost.ocks.org/mike/):

- [Makefile](https://github.com/mbostock/world-atlas/blob/master/Makefile)
  for his [World Atlas](https://github.com/mbostock/world-atlas)

- [Makefile](https://github.com/mbostock/us-atlas/blob/master/Makefile)
  for his [U.S. Atlas](https://github.com/mbostock/us-atlas)

- [Makefile](https://github.com/mbostock/d3/blob/master/Makefile) for [D3](http://d3js.org/)

### Resources

- [GNU make webpage](http://www.gnu.org/software/make)

- [Official manual](http://www.gnu.org/software/make/manual/make.html)

- O'Reilly
  [Managing projects with GNU make](http://oreilly.com/catalog/make3/book/)
  book (part of the [Open Books project](http://oreilly.com/openbook/))

- [Software carpentry](http://software-carpentry.org/)'s [make tutorial](http://software-carpentry.org/4_0/make/index.html)

- [Mike Bostock](http://bost.ocks.org/mike/)'s &ldquo;[Why Use Make](http://bost.ocks.org/mike/make/)&rdquo;

Also see my [git/github guide](http://kbroman.github.io/github_tutorial).

