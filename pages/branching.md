---
layout: page
title: Branching and merging
description: Git/GitHub Guide: branching and merging
---

I touched on just a few things about git.  Get yourself going with git and
github and then start looking at some of the many
[resources](resources.html).

In particular, pay attention to: [Branching](http://git-scm.com/book/en/Git-Branching-Basic-Branching-and-Merging):
Create a separate _branch_ to develop a feature (or work on a
bug) without disturbing the _master_ branch.  If it works out, you
can merge it back into the master; if it doesn't, you can trash it.
Branching is super easy, so for big projects, you should probably do it more
often than not.

To create a branch called `new_feature`:

    $ git branch new_feature

Then &ldquo;check it out&rdquo;:

    $ git checkout new_feature
    
Make various modifications, and then add and commit.

To go back to the master branch, check it out:

    $ git checkout master
    
To push the branch to github, use this:

    $ git push origin new_feature

If you make changes to the master branch, you'll want to merge them
into your exploratory one:

    $ git checkout new_feature
    $ git merge master
    
If you're satisfied with your changes in the exploratory branch, merge
them into the master:

    $ git checkout master
    $ git merge new_feature
    
If you're done with the branch and want to delete it:

    $ git branch -d new_feature

But if you pushed it to github, it will still exist there.  This is
how to delete the branch from github:

    $ git push origin --delete new_feature
    
After pulling from github, use the following to get access to a branch
that is only on github:

    $ git checkout -b new_feature origin/new_feature
    
If you want to pull a particular branch from a collaborator's
repository, do this:

    $ git checkout new_feature
    $ git pull myfriend new_feature

One final point: note that
[`git pull`](https://www.kernel.org/pub/software/scm/git/docs/git-pull.html)
is really doing a
[`git fetch`](https://www.kernel.org/pub/software/scm/git/docs/git-fetch.html)
followed by a
[`git merge`](https://www.kernel.org/pub/software/scm/git/docs/git-merge.html).

**Next**: [Git/github with RStudio](rstudio.html)
