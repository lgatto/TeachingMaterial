---
layout: page
title: Amend the last commit message
description: How to amend the last git commit message
---

&ldquo;Oops!  That last commit message was screwed up.
[How do I edit an incorrect commit message?](http://stackoverflow.com/questions/179123/how-do-i-edit-an-incorrect-commit-message-in-git)&rdquo;

This usually happens to me when I intended to use just `git commit`
but typed `git commit -a` and committed a whole bunch more stuff
that I hadn't mentioned in the commit message.

It's easy to fix just the message for the last commit:

    $ git commit --amend -m "New commit message"

Or leave off the `-m "New commit message"` and type the message in the
text editor that opens.

**Next**: [Exploring code and its history](exploring_code.html)
