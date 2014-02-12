---
layout: page
title: git/github with RStudio
---

[RStudio](http://www.rstudio.com/ide) is a popular development
environment for [R](http://www.r-project.org). 

RStudio has built-in facilities for [knitr](http://yihui.name/knitr/)
(the superb successor to [Sweave](http://www.stat.uni-muenchen.de/~leisch/Sweave/)) and [Rmarkdown](http://www.rstudio.com/ide/docs/r_markdown) (a
variant of
[Markdown](http://daringfireball.net/projects/markdown/syntax)), which
are highly recommended for writing data analysis reports.

And RStudio has built-in facilities for [git](http://git-scm.com/) and
[GitHub](http://www.github.com).  See
[this document on how to use version control with RStudio](http://www.rstudio.com/ide/docs/version_control/overview). It's
pretty straightforward.

Basically, in RStudio you want to create a Project
([described further here](http://www.rstudio.com/ide/docs/using/projects)),
which is basically a directory with some special files to describe
project-specific RStudio options. This Project will be your git
repository. Or you can easily make a current git repository into an
RStudio Project. And then you can use items in the menu bar to commit
changes and push or pull from GitHub.

**Next**: [Other resources](resources.html)
