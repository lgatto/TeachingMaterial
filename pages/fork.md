---
layout: page
title: Contribute to someone's repository
---

Say you want to contribute changes to someone else's repository (eg,
[this one](http://github.com/kbroman/github_tutorial)).

- Go to the repository on github.  (Say it's by `myfriend`, and is
  called `the_repo`, then you'll find it at `http://github.com/myfriend/the_repo`.

- Click the &ldquo;Fork&rdquo; button at the top right.

- You'll now have your own copy of that repository in your github account.

- Open a terminal/shell.

- Type

    ````
    $ git clone git@github.com:username/the_repo.git
    ````

- You'll now have a local copy of _your version_ of that repository.
- Add a connection to the original owner's repository.

    ````
    $ git remote add myfriend git://github.com/myfriend/the_repo.git
    ````

- Note the distinction between `git@github...` in the first case and
  `git://github...` in the second case.  I'm not sure why these need
  to be the way they are, but that's what works for me.

- Make changes to files.

- `git add` and `git commit` those changes

- `git push` them back to [github](http://github.com).  These will go
  to _your version_ of the repository.

- Go to _your version_ of the repository on github.

- Click &ldquo;Pull Request&rdquo; button at the top.

- Note that your friend's repository will be on the left and _your
  repository_ will be on the right.

- Give a short explaination of the changes and click the &ldquo;Send
  pull request&rdquo; button.


### Pulling others' changes

Before you make further changes to the repository, you should check
that your version is up to date relative to your friend's version.

Go into the directory for the project and type:

    $ git pull myfriend master

This will pull down and merge all of the changes that your friend has
made.

Now push them back to your github repository.

    $ git push
