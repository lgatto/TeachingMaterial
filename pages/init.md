---
layout: page
title: Start a new git repository
---

Your first instinct, when you start to do something new, should be
`git init`.  You're starting to write a new paper, you're writing a
bit of code to do a computer simulation, you mucking around with some
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
[markdown](http://daringfireball.net/projects/markdown/), describing
the project.


### A new repo from an existing project

Say you've got an existing project on which you've not been using
version control, and you want to start.

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

    $ git remote add origin git@github.com:username/new_repo.git
    $ git push -u origin master

Actually, the first line of the instructions will say

    $ git remote add origin https://github.com/username/new_repo.git
    
But I've found that the `https://github...` bit doesn't work and I
need to use `git@github...`  I don't really understand the difference.
