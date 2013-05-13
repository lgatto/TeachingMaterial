---
layout: page
title: Exploring code and its history
---

For my own repositories, I like to use [gitx](http://gitx.frim.nl/) on
my Mac to explore the commit history.

But, of course, you can also just poke around
[github](http://github.com) itself.

And you can do a lot of things from the command line, say to look at
the combined changes to a file that were made over the last year.
See the
[github page about `git diff`](http://learn.github.com/p/diff.html).

### Tags

Each commit is assigned a &ldquo;hash tag&rdquo; which is a unique
sequence of letters and numbers, like
[4d9fe11e56652bd19e19e28eac3906f09d5a3074](https://github.com/kbroman/github_tutorial/commit/4d9fe11e56652bd19e19e28eac3906f09d5a3074).
When you refer to these hash tags, you can just use an initial substr,
like `4d9fe`, that is unique to your repository.

For my [R/qtl package](http://github.com/kbroman/qtl), I like to tag
particular commits by the version number of the package, then I can
use my assigned tag in place of the less memorable hash tag.

To assign a tag, use something like

    git tag -a -m "Tagging version 1.28-5" 1.28-5

To push the tags to github, you need to use

    git push --tags

To delete a tag, use

    git tag -d 1.28-5

and then you need to remove the tag from github

    git push origin :refs/tags/1.28-5

### Uses of diff

To see all of the changes since the last commit, type

    git diff

To see all of the changes since a given commit, type

    git diff [commit]

Where in place of `[commit]` you use the initial part of a hash tag,
like `9f4668c`, or a tag you've created, like `1.28-1`.

To see all of the changes to a given file since a given commit, type

    git diff [commit] [file]

For example,

    git diff 1.22-21 R/scanone.R

To see all of the changes between two revisions, type something like

    git diff 1.22-21 1.23-16

If you use [gitx](http://gitx.frim.nl/), you can use it to view the
differences using a pipe:

    git diff 1.22-21 1.23-16 | gitx
