## Conclusions and further information

We've seen how we can use version control to,

* Keep track of changes like a lab notebook for code and documents.
* Roll back changes to any point in the history of changes to our files - "undo" and "redo" for files.
* Back up our entire history of changes in various locations.
* Work on our files from multiple locations.
* Identify and resolve conflicts when the same file is edited within two repositories without losing any work.
* Collaboratively work on code or documents or any other files.

Now, consider again our initial questions,

* If someone asks you, "can I have the code you used to create the data that you graphed in your conference paper?" You can use your version control logs and tags to immediately retrieve the exact version of the code that you used.
* If someone tells you, "er, your laptop's just been stolen!" ideally you'll have lost no work as you push it regularly to a remote repository.
* If you're working with a colleague on a journal paper who storms into your office and shouts, "you've just deleted my analysis section" you can just ask them to retrieve the previous version from the repository.
* If you're working with colleagues on a code and you find that a function you wrote has been rewritten and you want to know why, you could just look through the logs at the commit messages and, hopefully, these'll explain why it was changed.

Version control serves as a log book for your software and documents, ideas you've explored, fixes you've made, refactorings you've done, false paths you've explored - what was changed, who by, when and why - with a powerful undo and redo feature!

It also allows you to work with others on a project, whether that be writing code or papers, down to the level of individual files, without the risk of overwriting and losing each others work, and being able to record and understand who changed what, when, and why. 

> "If you are not using version control then, whatever else you may be doing with a computer, you are not doing science" -- Greg Wilson

### Find out more...

* K. Ram  (2013) "git can facilitate greater reproducibility and increased transparency in science", Source Code for Biology and Medicine 2013, 8:7 doi:[10.1186/1751-0473-8-7](http://dx.doi.org/10.1186/1751-0473-8-7) - survey of the range of ways in which version control can help research.
* [Visual Git Reference](http://marklodato.github.com/visual-git-guide/index-en.html) - pictorial representations of what Git commands do.
* [Pro Git](http://git-scm.com/book) - the "official" online Git book.
* [Version control by example](http://www.ericsink.com/vcbe/) - an acclaimed online book on version control by Eric Sink.
* [Git commit policies](http://osteele.com/posts/2008/05/commit-policies) - images on what Git commands to with reference to the working directory, staging area, local and remote repositories.
* [Gitolite](https://github.com/sitaramc/gitolite) - a way for you to host
your own multi-user Git repositories. Your collaborators send you their public SSH keys then they can pull and push from/to the repositories.
* G. Wilson, D. A. Aruliah, C. T. Brown, N. P. Chue Hong, M. Davis, R. T. Guy, S. H. D. Haddock, K. Huff, I. M. Mitchell, M. Plumbley, B. Waugh, E. P. White, P. Wilson (2012) "[Best Practices for Scientific Computing](http://arxiv.org/abs/1210.0530)", arXiv:1210.0530 [cs.MS].

Previous: [Collaborating with our colleagues](3_Collaboration.md) Next: [Git hints and tips](5_HintsAndTips.md)
