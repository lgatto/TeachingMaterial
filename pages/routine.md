---
layout: page
title: Routine use of git and github
---

The routine use of [git](http://git-scm.com) involves just a few commands:
principally `add`, `commit`, and `push`, but also `status` and
`diff`. 

You can deal with git and github via a [GUI](http://mac.github.com/),
but I prefer the command line, and so that's all I'll discuss.

### Add and commit

After you've made some small modifications to your project and
checked that they work, use `git add` to indicate that they're ready.
  
    $ git add R/modified.R man/modified.Rd

Then use `git commit` to add the modifications to the repository.

    $ git commit

A text editor (e.g., [emacs](http://www.gnu.org/software/emacs)) will
open; add a short message describing the changes.

To abandon your commit, exit the editor without adding text.

Note that `git add` is used to add completely new files as well as to
&ldquo;add&rdquo; modifications to files that already exist in the
repository.  

The commit message should be short (40 or 60 characters) so it's easy
to read in a list.  For a more complex commit, write an initial line
that is short and gives the overall idea, followed by as many lines as
you want giving the details.

People tend to write commit messages in the present rather than past tense
(eg, &ldquo;Fix such and such&rdquo; rather than &ldquo;Fixed such and
such&rdquo;).

For a one-line commit message, you can skip the text editor business
and just type

    $ git commit -m "Fix such and such"

### Add everything

If you want to commit all of the modifications you've made, without
having to explicitly &ldquo;add&rdquo; each file, you can skip the
separate `add` and `commit` commands and just type

    $ git commit -a

I try to avoid this, as it can lead to mistakes (committing more
modifications than intended).

### Push to [github](http://github.com)

To push committed changes to github, type

    $ git push

You don't need to do this every time. Do it after you've completed a batch
of changes that you're thoroughly happy with and before you move on to
something else.

Once you've pushed a commit, it's hard to take it away.  If you've
not pushed it yet, you _can_ go back and scrap it and not have it be part
of your project's history.

But if you move on to something else without having pushed the
changes, they may not get to github for months.


### Status

You've made some changes to a project, but you're not sure what.  Type

    git status
    
It'll give you a list of files that have been changed, plus new
files that haven't been formally added.


### Diff

Exactly what changes have you made?  Type

    git diff
    
Or to see your changes to a particular file, type

    git diff R/modified.R

It'll show you which lines have been added and which have been
deleted.


### How often to commit?

I prefer to do many small commits, each for a set of related changes:

- Think of something that needs to be fixed, or a feature to add.
- Do the work.
- Test that it is okay.
- Add and commit.

Look at others' projects on github, to see what they do and what sort
of commit messages they write.

### What to commit?

Don't include files that are derived from other files in the repository.  (Are you using
[make](http://www.gnu.org/software/make/) or
[rake](http://rake.rubyforge.org/)?  You should be!  See my
[make tutorial](kbroman.github.io/minimal_make).)

For example, for a [LaTeX](http://www.latex-project.org/) manuscript,
I wouldn't include all the .log, .dvi, .aux, etc., files.  And if I
have R code to generate a figure, I'll include the R code but not the
figure.

Be careful about committing binary files, or really big files.  Git
works best with text files (like source code), as you can see
just the lines that were changed.  A new copy of a file will get added to
the repository every time you change it. For small text files, that's
no big deal; for big images, you'll get a bloated repository.

And once you've committed a big file to your repository, it's there
forever, even if you use `git rm` to remove it later.

For big data files that are changing, you'll want to track a
text-based version (not .xls!), and you may want to make a
fully separate git repository for the data.

### .gitignore

The various files in your project directory that you're not tracking
in git should be indicated in a `.gitignore` file.

You don't _have_ to have a `.gitignore` file, but if you don't, those
files will show up every time you type `git status`.

Each subdirectory can have its own `.gitignore` file, too.  

Also, you can have a global such in your home directory; I use
`~/.gitignore_global`, which contains:

    *~
    .*~
    .DS_Store
    .Rhistory
    .RData

You have to tell git about the global `.gitignore` file:

    $ git config --global core.excludesfile ~/.gitignore_global

**Next**: [Start a new repository](init.html)
