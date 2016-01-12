# Best practices in bioinformatics research: open source software and reproducibility

![See, don't trust](./figs/CPcf52hWUAEDlF8.jpg)

#### Quick survey:

- What does *open* mean? What should be, or must be open? Is it that way?

- What is *reproducibility*? Why do we want reproducibility? Is it that way?

# Open Science

Science is supported by **public** funds. All outputs should this be
made publicly and openly available - to the research community and the
general public.

Who owns the research outputs?

**Open** means free from any fees (*free as in beer*) and free from
any re-use restrictions (*free as in speech).

#### Open data

The data, both in raw and processed form, and its annotation
(metadata) that accompanies a scientific paper (discovery) should be
freely/publicly available and should be free to re-use.

#### Open methodology

> An open methodology is simply one which has been described in
> sufficient detail to allow other researchers to repeat the work and
> apply it elsewhere (from Waton M., 2015)

Describe and release the process that lead raw data to processed data,
and how the processed has been further analysed to lead to results,
figures and conclusions.

- **[Open source](http://opensource.org/osd)**: the source code of
  your (scienfic) software should be openly released, so that others
  can read/understand/fix/re-use it. (see
  [software licences](http://opensource.org/licenses))

- **Open protocols**: openly share the hole sample processing, data
  acquisition.

#### Open access

Papers and supplementary material should be available to **read** and
**mine** by all, researchers and public.

There exist
[open access licenses](https://creativecommons.org/licenses/): CC-BY,
CC-BY-SA, ...

Also: open peer review, open education, ... 

See also 

- [Why Open Research?](http://whyopenresearch.org/)
- [The open research value proposition: How sharing can help researchers succeed](https://figshare.com/articles/The_open_research_value_proposition_How_sharing_can_help_researchers_succeed/1619902)

If you are interested in *All things open* (science, education, ...),
consider following/joining [OpenConCam](http://www.openconcam.org/)

# Reproducible research

> D Knuth: [Literate programming](http://literateprogramming.com/) is a
> methodology that combines a programming language with a
> documentation language, thereby making programs more robust, more
> portable, more easily maintained, and arguably more fun to write
> than programs that are written only in a high-level language.

<br />

> [D. Donoho](http://statweb.stanford.edu/~donoho/): An article about
> computational science in a scientific publication is not the
> scholarship itself, it is merely advertising of the scholarship. The
> actual scholarship is the complete software development environment
> and the complete set of instructions which generated the figures.

<br />

> R. Gentleman and D. Temple Land,
> [Statistical Analyses and Reproducible Research](http://biostats.bepress.com/bioconductor/paper2/) 2004. We
> introduce the concept of a compendium as both a container for the
> different elements that make up the document and its computations
> (i.e. text, code, data, ...), and as a means for distributing,
> managing and updating the collection.

## Nomenclature


- Reproducibility/reproduce


- Replication/replicate

Other terms: mechanical reproducibility, repeat/repeatability, re-use 


## Why 

> Open science isn't a movement, it's just (good) science. It's also
> the future. (from Watson M, 2015)

#### Open Science

- Moral argument
- Better science

But also

#### [Five selfish reasons to work reproducibly](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0850-7)

1. Reproducibility helps to avoid disaster
2. Reproducibility makes it easier to write papers
3. Reproducibility helps reviewers see it your way
4. Reproducibility enables continuity of your work
5. Reproducibility helps to build your reputation





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
4. Don't repeat yourself (or others).
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

## Conclusions

- Open, reproducibility

- [When will *open science* become simply *science*?](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0669-2)

- Traceability, transparency

- [This can not be an afterthought](http://www.ncbi.nlm.nih.gov/pubmed/26191404).

- Make your data and code suistainable

> A piece of software is being sustained if people are using it,
> fixing it, and improving it rather than replacing
> it [ref](http://software-carpentry.org/blog/2014/08/sustainability.html).

[The Software Sustainability Institute](http://software.ac.uk/)

## References

- Gentleman, Robert and Temple Lang, Duncan, "Statistical Analyses and
  Reproducible Research" (May 2004). Bioconductor Project Working
  Papers. [Working Paper 2](http://biostats.bepress.com/bioconductor/paper2).

- Wilson G, Aruliah DA, Brown CT, Chue Hong NP, Davis M, Guy RT, et
  al. (2014) Best Practices for Scientific Computing. PLoS Biol 12(1):
  e1001745. [doi:10.1371/journal.pbio.1001745](http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001745).

- Yihui Xie,
  [Dynamic Documents with R and knitr](https://github.com/yihui/knitr-book),
  Second Edition (Chapman & Hall/CRC The R Series) 2nd Edition, 2015.

- Victoria Stodden, Friedrich Leisch, Roger D. Peng,
  [Implementing Reproducible Research](https://osf.io/s9tya/) (Chapman
  & Hall/CRC: The R Series), 2014.

- Florian Markowetz, Five selfish reasons to work reproducibly. Genome
  Biology 2015 16:274
  [DOI:10.1186/s13059-015-0850-7](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0850-7).

- Leek JT, Peng RD. Opinion: Reproducible research can still be wrong:
  adopting a prevention approach. Proc Natl Acad Sci U S A. 2015 Feb
  10;112(6):1645-6. [doi:10.1073/pnas.1421412111](http://www.pnas.org/content/112/6/1645).


- Donoho DL. An invitation to reproducible computational
  research. Biostatistics. 2010
  Jul;11(3):385-8. [doi:10.1093/biostatistics/kxq028](http://biostatistics.oxfordjournals.org/content/11/3/385.long).

- Watson M. When will 'open science' become simply 'science'? Genome
  Biol. 2015 May
  19;16:101. [doi:10.1186/s13059-015-0669-2](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0669-2).
