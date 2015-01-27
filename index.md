---
layout: page
title: minimal make
tagline: A minimal tutorial on make
description: A minimal tutorial on GNU Make, the most important tool for reproducible research.
---

I would argue that the most important tool for reproducible research
is not [Sweave](http://www.stat.uni-muenchen.de/~leisch/Sweave/) or
[knitr](http://yihui.name/knitr/) but
*[GNU make](http://www.gnu.org/software/make)*.

Consider, for example, all of the files associated with a manuscript.
In the simplest case, I would have an [R](http://r-project.org)
script for each figure plus a [LaTeX](http://www.latex-project.org)
file for the main text.  And then a [BibTeX](http://www.bibtex.org)
file for the references.

Compiling the final PDF is a bit of work:

- Run each R script through R to produce the relevant figure.
- Run latex and then bibtex and then latex a couple of more times.

And the R scripts need to be run before latex is, and only if they've
changed.

### A simple example

[GNU make](http://www.gnu.org/software/make) makes this easy.  In your
directory for the manuscript, you create a text file called `Makefile`
that looks something like [the following](examples/ex1/Makefile) (here using
[pdflatex](http://www.tug.org/applications/pdftex/)).

    mypaper.pdf: mypaper.bib mypaper.tex Figs/fig1.pdf Figs/fig2.pdf
        pdflatex mypaper
        bibtex mypaper
        pdflatex mypaper
        pdflatex mypaper

    Figs/fig1.pdf: R/fig1.R
        cd R;R CMD BATCH fig1.R

    Figs/fig2.pdf: R/fig2.R
        cd R;R CMD BATCH fig2.R

Each batch of lines indicates a file to be created (the _target_), the files it
depends on (the _prerequisites_), and then a set of commands needed to
construct the target from the dependent files.  Note that the lines
with the commands _must_ start with a **tab** character (**not spaces**).

Another great feature: in the example above, you'd only build
`fig1.pdf` when `fig1.R` changed.  And note that the dependencies
propagate.  If you change `fig1.R`, then `fig1.pdf` will change, and
so `mypaper.pdf` will be re-built.

**One oddity**: if you need to change directories to run a command, do
the `cd` on the same line as the related command.  The following
**_would not work_**:

    ### this doesn't work ###
    Figs/fig1.pdf: R/fig1.R
        cd R
        R CMD BATCH fig1.R

You can, however, use `\` for a continuation line, line so:

    ### this works ###
    Figs/fig1.pdf: R/fig1.R
        cd R;\
        R CMD BATCH fig1.R

Note that you still need to use the semicolon (`;`).

#### Using GNU make

You probably already have GNU make installed on your computer.  Type
`make --version` in a terminal/shell to see. (On Windows,
[go here to download make](http://gnuwin32.sourceforge.net/packages/make.htm).)

To use make:

- Go into the the directory for your project.
- Create the `Makefile` file.
- Every time you want to build the project, type `make`.
- In the example above, if you want to build `fig1.pdf` without
  building `mypaper.pdf`, just type `make fig1.pdf`.

### Frills

You can go along way with just simple make files as above, specifying
the target files, their dependencies, and the commands to create
them. But there are _a lot_ of frills you can add, to save some
typing.

Here are some of the options that I use. (See the
[make documentation](http://www.gnu.org/software/make/manual/make.html)
for further details.)

#### Variables

If you'll be repeating the same piece of code multiple times, you
might want to define a
[variable](http://www.gnu.org/software/make/manual/make.html#Using-Variables).

For example, you might want to run R with the flag `--vanilla`. You
could then define a variable `R_OPTS`:

    R_OPTS=--vanilla

You refer to this variable as `${R_OPTS}`, so in the R commands you
would use something like

    cd R;R CMD BATCH ${R_OPTS} fig1.R

An advantage of this is that you just need to type out the options you
want once; if you change your mind about the R options you want to
use, you just have to change them in the one place.

For example, I actually like to use the following:

    R_OPTS=--no-save --no-restore --no-init-file --no-site-file

This is like `--vanilla` but without `--no-environ` (which I need
because I use the `.Renviron` file to define `R_LIBS`, to say that I
have R packages defined in an alternative directory).


#### Automatic variables

There are a bunch of
[automatic variables](http://www.gnu.org/software/make/manual/make.html#Automatic-Variables)
that you can use to save yourself a lot of typing. Here are the ones
that I use most:

- `$@` &nbsp;&nbsp; the file name of the target
- `$<` &nbsp;&nbsp; the name of the first prerequisite (i.e., dependency)
- `$^` &nbsp;&nbsp; the names of all prerequisites (i.e., dependencies)
- `$(@D)` &nbsp;&nbsp; the directory part of the target
- `$(@F)` &nbsp;&nbsp; the file part of the target
- `$(<D)` &nbsp;&nbsp; the directory part of the first prerequisite (i.e., dependency)
- `$(<F)` &nbsp;&nbsp; the file part of the first prerequisite (i.e., dependency)

For example, in our simple example, we could simplify the lines

    Figs/fig1.pdf: R/fig1.R
        cd R;R CMD BATCH fig1.R

We could instead write

    Figs/fig1.pdf: R/fig1.R
        cd $(<D);R CMD BATCH $(<F)

The automatic variable `$(<D)` will take the value of the directory of
the first prerequisite, `R` in this case. `$(<F)` will take value of
the file part of the first prerequisite, `fig1.R` in this case.

Okay, that's not _really_ a simplification.
There doesn't seem to be much advantage to this, unless perhaps the
directory were an obnoxiously long string and we wanted to avoid having
to type it twice. The main advantage comes in the next section.

#### Pattern rules

If a number of files are to be built in the same way, you may want to
use a
[pattern rule](http://www.gnu.org/software/make/manual/make.html#Pattern-Rules).
The key idea is that you can use the symbol `%` as a wildcard, to be
expanded to any string of text.

For example, our two figures are being built in basically the same
way. We could simplify the example by including one set of lines
covering both `fig1.pdf` and `fig2.pdf`:

    Figs/%.pdf: R/%.R
        cd $(<D);R CMD BATCH $(<F)

This saves typing and makes the file easier to maintain and extend. If
you want to add a third figure, you just add it as another dependency
(i.e., prerequisite) for `mypaper.pdf`.

#### Our example, with the frills

Adding all of this together, here's what our example `Makefile`
[will look like](examples/ex2/Makefile).

    R_OPTS=--vanilla

    mypaper.pdf: mypaper.bib mypaper.tex Figs/fig1.pdf Figs/fig2.pdf
        pdflatex mypaper
        bibtex mypaper
        pdflatex mypaper
        pdflatex mypaper

    Figs/%.pdf: R/%.R
        cd $(<D);R CMD BATCH $(R_OPTS) $(<F)

The advantage of the added frills: less typing, and it's easier to
extend to include additional figures. The disadvantage: it's harder
for others who are less familiar with
[GNU Make](http://www.gnu.org/software/make/) to understand what it's
doing.

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

- [Makefile](https://github.com/kbroman/qtlcharts/blob/master/Makefile)
  for my [R/qtlcharts](http://kbroman.org/qtlcharts) package.

And here are some examples from [Mike Bostock](http://bost.ocks.org/mike/):

- [Makefile](https://github.com/mbostock/world-atlas/blob/master/Makefile)
  for his [World Atlas](https://github.com/mbostock/world-atlas)

- [Makefile](https://github.com/mbostock/us-atlas/blob/master/Makefile)
  for his [U.S. Atlas](https://github.com/mbostock/us-atlas)

- [Makefile](https://github.com/mbostock/d3/blob/master/Makefile) for [D3](http://d3js.org/)

Also look at the
[Makefile](https://github.com/yihui/knitr/blob/master/Makefile) for
[Yihui Xie](http://yihui.name/)'s [knitr](http://yihui.name/knitr/) package for [R](http://r-project.org).

Also of interest is [`maker`](https://github.com/ComputationalProteomicsUnit/maker), a `Makefile` for R package development.

### Resources

- [GNU make webpage](http://www.gnu.org/software/make)

- [Official manual](http://www.gnu.org/software/make/manual/make.html)

- O'Reilly
  [Managing projects with GNU make](http://oreilly.com/catalog/make3/book/)
  book (part of the [Open Books project](http://oreilly.com/openbook/))

- [Software carpentry](http://software-carpentry.org/)'s [make tutorial](http://software-carpentry.org/v4/make/index.html)

- [Mike Bostock](http://bost.ocks.org/mike/)'s &ldquo;[Why Use Make](http://bost.ocks.org/mike/make/)&rdquo;

- [GNU Make for reproducible data analysis](http://zmjones.com/make.html) by [Zachary Jones](http://zmjones.com/)

- [Makefiles for R/LaTeX projects](http://robjhyndman.com/hyndsight/makefiles/) by [Rob Hyndman](http://robjhyndman.com)

---

The source for this minimal tutorial is [on github](http://github.com/kbroman/minimal_make).

Also see my [git/github guide](http://kbroman.org/github_tutorial),
[knitr in a knutshell tutorial](http://kbroman.org/knitr_knutshell),
[R package primer](http://kbroman.org/pkg_primer),
[simple site tutorial](http://kbroman.org/simple_site),
and [initial steps towards reproducible research](http://kbroman.org/steps2rr).
