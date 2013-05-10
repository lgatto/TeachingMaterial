---
layout: page
title: Handling merge conflicts
---

One of the best features of [git](http://git-scm.com) is its ability
to easily merge multiple changes by different people.  

Say you and a
friend have both made changes to the same file at the same time.  
When you pull your friend's changes, git will often be able to combine
them without any problem.  

Sometimes, though, after you do

    $ git pull myfriend master

You'll get a message like

    Auto-merging README.md
    CONFLICT (content): Merge conflict in README.md
    Automatic merge failed; fix conflicts and then commit the result.

If you open the offending file in a text editor, you'll find an indication
of the bits that are different, something like this:

    <<<<<<< HEAD
    A line in my file.
    =======
    A line in my friend's file
    >>>>>>> 031389f2cd2acde08e32f0beb084b2f7c3257fff

Edit the bits from `<<<<<<<` to `>>>>>>>`, to make the file just as
you want it

Then do `git add`, `git commit`, and `git push`.
