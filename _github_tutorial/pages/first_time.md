---
layout: page
title: Your first time with git and github
---

If you've never used git or github before, there are a bunch of things
that you need to do.  It's
[very well explained on github](http://help.github.com/articles/set-up-git),
but repeated here for completeness.

- Get a [github](http://github.com) account.
- Download and install [git](http://git-scm.com/downloads).
- Set up git with your user name and email.

  - Open a terminal/shell and type:
  
    ````
    $ git config --global user.name "Your name here"
    ````

    ````
    $ git config --global user.email "your_email@example.com"
    ````

- Set up ssh on your computer.  I like
  [Roger Peng](http://www.biostat.jhsph.edu/~rpeng)'s
  [guide to setting up password-less logins](http://www.biostat.jhsph.edu/bit/nopassword.html).
  Also see [github's guide to generating SSH keys](http://help.github.com/articles/generating-ssh-keys).

  - Look to see if you have files `~/.ssh/id_rsa` and
  `~/.ssh/id_rsa.pub`.
  - If not, create such public/private keys: Open a terminal/shell and type:
  
    ````
    $ ssh-keygen -t rsa -C "your_email@example.com"
    ````

  - Copy your public key (the contents of the newly-created
    `id_rsa.pub` file) into your clipboard.  On a Mac, in the terminal/shell, type:
  
    ````
    $ pbcopy < ~/.ssh/id_rsa.pub
    ````

- Paste your ssh public key into your github account settings.

  - Go to your github [Account Settings](http://github.com/settings/profile)
  - Click &ldquo;[SSH Keys](http://github.com/settings/ssh)&rdquo; on the left.
  - Click &ldquo;Add SSH Key&rdquo; on the right.
  - Add a label (like &ldquo;My laptop&rdquo;) and paste the public
    key into the big text box.
  - In a terminal/shell, type the following to test it:
  
    ````
    $ ssh -T git@github.com
    ````
    
  - If it says something like the following, it worked:
  
    ````
    Hi username! You've successfully authenticated, but Github does
    not provide shell access.
    ````
