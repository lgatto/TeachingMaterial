## Collaborating with our colleagues

### Collaboration: pulling changes from a remote repository

Now when we have a remote repository, we can share it and collaborate with others (and we can also work from multiple locations: for example from a laptop and a desktorp in the lab). But how do we get the latest changes? One way is simply to clone the repository every time but this is inefficient, especially if our repository is very large. So, git allows us to get the latest changes down from a repository. 

Find a partner to work in a pair. We will work, in pairs, on one of the Makefiles used yesterday. One person should start with creating a repository on their machine. Leave the  current local repository and let us clone the `rbc` repository from Github  but this time specify the local directory name,

    $ git clone https://github.com/lgatto/rbc.git 2014-01-Cambridge
    Cloning into '2014-01-Cambridge'...

Now create an empty repository on Github and push the cloned repository.

The second person can now clone that repository from Github on their own machine. In order to enable the second person to push back to the repository, its owner must add that person's Github name to 'Collaborators' (under 'Settings' on the right hand side).  

Now both of you should have the same version of the repository on your machines. One person should search for the Makefile with the typo and fix it. Then commit the change and push it to the remote repository on Github. 

The other person can now *fetch* the changes from the remote repository,

    $ git fetch

We can now see what the differences are by doing,

    $ git diff origin/master

which compares our current, `master` branch, with an `origin/master` branch which is the name of the `master` branch in `origin` which is the alias for our cloned repository, the one on GitHub.

We can then *merge* these changes into our current repository, which merges the branches together,

    $ git merge origin/master

And then we can check that we have our changes,

    $ cat Makefile
    $ git log

As a short-hand, we can do a Git *pull* which does a *fetch* then a *merge*,

	$ git pull

   
### Collaboration: conflicts and how to resolve them

Let's do some further work on the Makefile. Try working independently on the same file (and possibly on the same part of the file) and push it to the remote repository. Inevitably, one of you will be the first one to push. When the other person tries to push, they should get a message similar to the following:

    $ git push
    To https://github.com/lgatto/rbc.git
 	! [rejected]        master -> master (non-fast-forward)
	error: failed to push some refs to 'https://github.com/lgatto/rbc.git'
	hint: Updates were rejected because the tip of your current branch is behind
	hint: its remote counterpart. Merge the remote changes (e.g. 'git pull')
	hint: before pushing again.
	hint: See the 'Note about fast-forwards' in 'git push --help' for details.


The push fails, as the other person needs to pull down the most recent changes in the repository. Before pushing we should always pull, so let's do that...

    $ git pull

If you actually edited the same part of the Makefile, you may get a conflict:

    Auto-merging Makefile
    CONFLICT (content): Merge conflict in Makefile
    Automatic merge failed; fix conflicts and then commit the result.

As we saw earlier, with the fetch and merge, a pull pulls down changes from the repository and tries to merge these. It does this on a file-by-file basis, merging files line by line. We get a *conflict* when if a file has changes that affect the same lines and those changes can't be seamlessly merged. If we look at the status,

    $ git status

we can see that our file is listed as `Unmerged` and if we look at the file that was modified, we will see the markers for the HEAD and the current commit, just like in the previous example. 

The mark-up shows us the parts of the file causing the conflict and the versions they come from. We now need to manually edit the file to *resolve* the conflict. Just like we did when we had to deal with the conflict when we were merging the branches.

We edit the file. Then commit our changes e.g. Now if we push,

    $ git push

All goes well. If we now go to GitHub and click on the Overview tab we can see where our repository diverged and came together again.

This is where version control proves itself better than DropBox or GoogleDrive, this ability to merge text files line-by-line and highlight the conflicts between them, so no work is ever lost.

Previous: [Working from multiple locations with a remote repository](2_Remote.md) Next: [Conclusions and further information](4_Conclusion.md)
