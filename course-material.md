# Best practices in bioinformatics research: open source software and reproducibility

- What does *open* mean?
- what should be, or must be open?
- What is *reproducibility*?
- Why do we want reproducibility?

## Open science

- open data
- open source
- open access

## Reproducible research


> D Knuth [Literate programming](http://literateprogramming.com/) is a
> methodology that combines a programming language with a
> documentation language, thereby making programs more robust, more
> portable, more easily maintained, and arguably more fun to write
> than programs that are written only in a high-level language.


> [D. Donoho](http://statweb.stanford.edu/~donoho/): An article about
> computational science in a scientific publication is not the
> scholarship itself, it is merely advertising of the scholarship. The
> actual scholarship is the complete software development environment
> and the complete set of instructions which generated the figures.


> R. Gentleman and D. Temple Land
> [Statistical Analyses and Reproducible Research](http://biostats.bepress.com/bioconductor/paper2/) 2004. We
> introduce the concept of a compendium as both a container for the
> different elements that make up the document and its computations
> (i.e. text, code, data, ...), and as a means for distributing,
> managing and updating the collection.

## Tools

#### R and Sweave/knitr

```
Rmd -> md -> html | pdf
```

```
Rnw -> tex -> pdf
```

#### Jupyter notebook

Previously known as IPython notebooks.

#### Other

`org-mode`


### Best practice 

Wilson G *et al.*
[Best Practices for Scientific Computing](http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001745)
2014.

1. Write programs for people, not computers.
   - A program should not require its readers to hold more than a handful of facts in memory at once.
   - Make names consistent, distinctive, and meaningful.
   - Make code style and formatting consistent.
2. Let the computer do the work.
   - Make the computer repeat tasks.
   - Save recent commands in a file for re-use.
   -  Use a build tool to automate workflows.
3. Make incremental changes.
   - Work in small steps with frequent feedback and course correction.
   - Use a version control system.
   - Put everything that has been created manually in version control.
4.Don't repeat yourself (or others).
   - Every piece of data must have a single authoritative representation in the system.
   - Modularize code rather than copying and pasting.
   - Re-use code instead of rewriting it.
5. Plan for mistakes.
   - Add assertions to programs to check their operation.
   - Use an off-the-shelf unit testing library.
   - Turn bugs into test cases.
   - Use a symbolic debugger.
6. Optimize software only after it works correctly.
   - Use a profiler to identify bottlenecks.
   - Write code in the highest-level language possible.
7. Document design and purpose, not mechanics.
   - Document interfaces and reasons, not implementations.
   - Refactor code in preference to explaining how it works.
   - Embed the documentation for a piece of software in that software.
8. Collaborate.
   - Use pre-merge code reviews.
   - Use pair programming when bringing someone new up to speed and when tackling particularly tricky problems.
   - Use an issue tracking tool.

## References



Gentleman, Robert and Temple Lang, Duncan, "Statistical Analyses and
Reproducible Research" (May 2004). Bioconductor Project Working
Papers. [Working Paper 2](http://biostats.bepress.com/bioconductor/paper2).


Wilson G, Aruliah DA, Brown CT, Chue Hong NP, Davis M, Guy RT, et
al. (2014) Best Practices for Scientific Computing. PLoS Biol 12(1):
e1001745. [doi:10.1371/journal.pbio.1001745](http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001745).



