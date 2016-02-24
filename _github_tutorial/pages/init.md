---
layout: page
title: Start a new git repository
description: Creating a new git repository
---

Your first instinct, when you start to do something new, should be
`git init`.  You're starting to write a new paper, you're writing a
bit of code to do a computer simulation, you're mucking around with some
new data ... _anything_: think `git init`.

### A new repo from scratch

Say you've just got some data from a collaborator and are about to
start exploring it.

- Create a directory to contain the project.
- Go into the new directory.
- Type `git init`.
- Write some code.
- Type `git add` to add the files (see the
  [typical use page](routine.html)).
- Type `git commit`.

The first file to create (and add and commit) is probably a ReadMe
file, either as plain text or with
[Markdown](http://daringfireball.net/projects/markdown/), describing
the project.

Markdown allows you to add a bit of text markup, like
[hyperlinks](http://en.wikipedia.org/wiki/Hyperlink),
**bold**/_italics_, or to indicate code with a `monospace
font`. Markdown is easily converted to html for viewing in a web
browser, and GitHub will do this for you automatically.



### A new repo from an existing project

Say you've got an existing project that you want to start tracking
with git.

- Go into the directory containing the project.
- Type `git init`.
- Type `git add` to add all of the relevant files.
- You'll probably want to create a `.gitignore` file right away, to
  indicate all of the files you don't want to track.  Use `git add
  .gitignore`, too.
- Type `git commit`.


### Connect it to github

You've now got a local git repository.  You can use git locally, like
that, if you want.  But if you want the thing to have a home on github, do
the following.

- Go to [github](http://github.com).
- Log in to your account.
- Click the [new repository](https://github.com/new) button in the
top-right.  You'll have an option there to initialize the repository with a README
file, but I don't.
- Click the &ldquo;Create repository&rdquo; button.

Now, follow the second set of instructions, &ldquo;Push an existing
repository...&rdquo;

    $ git remote add origin git@github.com:username/new_repo
    $ git push -u origin master

Actually, the first line of the instructions will say

    $ git remote add origin https://github.com/username/new_repo

But I use `git@github.com:username/new_repo` rather than `https://github.com/username/new_repo`, as the
former is for use with
[ssh](http://en.wikipedia.org/wiki/Secure_Shell) (if you set up ssh as
I mentioned in "[Your first time](first_time.html)", then you won't
have to type your password every time you push things to github). If
you use the latter construction, you'll have to type your github
password every time you push to github.

**Next**: [Contribute to someone's repository](fork.html)
