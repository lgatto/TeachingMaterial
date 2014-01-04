## Working from multiple locations with a remote repository

We're going to set up a remote repository that we can use from multiple locations. The remote repository can also be shared with colleagues, if we want to.

### GitHub

[GitHub](http://GitHub.com) is a web service which allows users to set up  their private and public source code Git repositories. It provides tools for browsing, collaborating on and documenting code. Your organisation may also offer support for hosting Git repositories - ask your local system administrator. GitHub, like other services such as [Launchpad](https://launchpad.net), [Bitucket](https://bitbucket.org),
[GoogleCode](http://code.google.com), and [SourceForge](http://sourceforge.net) provides a wealth of resources to support projects including:

* Time histories changes to repositories
* Commit-triggered e-mails
* Browsing code from within a web browser, with syntax highlighting
* Software release management
* Issue (ticket) and bug tracking
* Download
* Varying permissions for various groups of users
* Other service hooks e.g. to Twitter.

**Note**  GitHub's free repositories have public licences **by default**. If you don't want to share (in the most liberal sense) your stuff with the world and you want to use GitHub (instead of GitHub), you will need to pay for the private GitHub repositories (GitHub offers up to 5 free private repositories, if you are an academic - but do check this information as T&C may change).

### Get an account

Setting up GitHub at first requires a user name and password. Let's now  [create a free one](https://GitHub.com). 

### Create a new repository

Now, we can create a repository on GitHub,

* Log in to [GitHub](https://GitHub.com/)
* Click on the Create icon on the top right
* Enter Repository name: "2014-01-Cambridge"
* For the purpose of this exercise we'll create a public repository
* Make sure the Initialize this repository with a README is *unselected*
* Click Create Repository

You'll get a page with new information about your repository. We already have our local repository and we will be *pushing* it to GitHub. 

    git remote add origin https://github.com/USERNAME/2014-01-Cambridge.git
    git push -u origin master

This sets up an alias, `origin`, to correspond to the URL of our new repository on GitHub.

Now copy and paste the second line,

    $ git push -u origin master
    Counting objects: 38, done.
    Delta compression using up to 4 threads.
    Compressing objects: 100% (30/30), done.
    Writing objects: 100% (38/38), 3.59 KiB, done.
    Total 38 (delta 9), reused 0 (delta 0)
    To https://github.com/USERNAME/2014-01-Cambridge.git
    * [new branch]      master -> master
    Branch master set up to track remote branch master from origin.


This *pushes* our `master` branch to the remote repository, named via the alias `origin` and creates a new `master` branch in the remote repository.

Now, on GitHub, we should see our code and click the Commits tab we should see our complete history of commits.  

Our local repository is now available on GitHub. So, anywhere we can access GitHub, we can access our repository.

### Cloning a remote repository

Now, let's do something drastic!

    $ cd ..
    $ rm -rf add_nums

We've just wiped our local repository! But, as we've a copy on GitHub we can just copy, or *clone* that,

    $ git clone https://github.com/USERNAME/2014-01-Cambridge.git
    Cloning into '2014-01-Cambridge'...
    Password for 'https://USERNAME@GitHub.com':
    remote: Counting objects: 38, done.
    remote: Compressing objects: 100% (21/21), done.
    remote: Total 38 (delta 9), reused 38 (delta 9)
    Unpacking objects: 100% (38/38), done.

In GitHub, there is also an option for called *forking* which essentially allows to do the same thing but forking is done in the web service (that is not using the command line and creating a repository on your local file system).

Now, if we change into `2014-01-Cambridge` we can see that we have our repository,

    $ cd 2014-01-Cambridge
    $ git log

and we can see our Git configuration files too,

    $ ls -A

But where is the `add_nums` directory, you might ask? `add_nums` was the directory that held our local repository but was not a part of it.

### Push changes to a remote repository

We can use our cloned repository just as if it was a local repository so let's make some changes to our files and commit these.

Having done that, how do we send our changes back to the remote repository? We can do this by *pushing* our changes,

    $ git push

If we now check our GitHub page we should be able to see our new changes under the Commit tab.

If we created a new branch to develop a new feature and then we want to push it into GitHub:

	$ git push origin feature2 

This will work assumig that 'origin' is still an alias for our remote repository in GitHub. The feature2 branch should now be created in our GitHub repository.

### Collaboration: pulling changes from a remote repository

Now when we have a remote repository, we can share it and collaborate with others (and we can also work from multiple locations: for example from a laptop and a desktorp in the lab). But how do we get the latest changes? One way is simply to clone the repository every time but this is inefficient, especially if our repository is very large. So, Git allows us to get the latest changes down from a repository. Let's assume that we work on our 

First, let us leave our current local repository,

    $ cd ..
    $ ls
    2014-01-Cambridge

And let us clone our repository again, but this time specify the local directory name,

    $ git clone https://github.com/USERNAME/2014-01-Cambridge.git bootcamp2
    Cloning into 'bootcamp2'...


So we now have two clones of our repository,

    $ ls
    $ bootcamp bootcamp2

Let's pretend these clones are on two separate machines! So we have 3 versions of our repository - our two local versions, on our separate machines (we're still pretending!) and one on GitHub. So let's go into one of our clones, make some changes, commit these and push these to GitHub:

    $ cd bootcamp
    $ nano add_numb.py
    $ git add add_numb.py
    $ git commit -m "Added some more comments" add_numb.py
    $ git push

Now let's change to our other repository and *fetch* the changes from our remote repository,

    $ cd ../bootcamp2
    $ git fetch

We can now see what the differences are by doing,

    $ git diff origin/master

which compares our current, `master` branch, with an `origin/master` branch which is the name of the `master` branch in `origin` which is the alias for our cloned repository, the one on GitHub.

We can then *merge* these changes into our current repository, which merges the branches together,

    $ git merge origin/master

And then we can check that we have our changes,

    $ cat add_numb.py
    $ git log

As a short-hand, we can do a Git *pull* which does a *fetch* then a *merge*,

    $ nano add_numb.py
    $ git add add_numb.py
    $ git commit -m "Added credits" add_numb.py
    $ git push
    $ cd ..
    $ cd ../2014-01-Cambridge
    $ git pull

And then check that we have our changes,

    $ cat add_numb.py
    $ git log

### Collaboration: conflicts and how to resolve them

Let's continue to pretend that our two local, cloned, repositories are hosted on two different machines, and make some changes to our file, and push these to GitHub:

    $ nano add_numb.py
    $ git add add_numb.py
    $ git commit -m "Credits - added our names" add_numb.py
    $ git push

Now let us suppose, at a later, date, we use our other repository and we want to change the credits.

    $ cd ../bootcamp2
    $ nano add_numb.py
    $ git add add_numb.py
    $ git commit -m "Changed the first author" add_numb.py
    $ git push

Our push fails, as we've not yet pulled down our changes from our remote repository. Before pushing we should always pull, so let's do that...

    $ git pull

and we get:

    Auto-merging add_numb.py
    CONFLICT (content): Merge conflict in add_numb.py
    Automatic merge failed; fix conflicts and then commit the result.

As we saw earlier, with the fetch and merge, a pull pulls down changes from the repository and tries to merge these. It does this on a file-by-file basis, merging files line by line. We get a *conflict* when if a file has changes that affect the same lines and those changes can't be seamlessly merged. If we look at the status,

    $ git status

we can see that our file is listed as `Unmerged` and if we look at `add_numb.py`, we may see something like,

    <<<<<<< HEAD 
    """ This short program just adds numbers. Developed by John and Aleksandra."""
    =======
    """ This short program just adds numbers. Developed by Aleksandra and John."""
    >>>>>>> 71d34decd32124ea809e50cfbb7da8e3e354ac26 

The mark-up shows us the parts of the file causing the conflict and the versions they come from. We now need to manually edit the file to *resolve* the conflict. Just like we did when we had to deal with the conflict when we were merging the branches.

We edit the file. Then commit our changes e.g. Now if we push,

    $ git push

All goes well. If we now go to GitHub and click on the Overview tab we can see where our repository diverged and came together again.

This is where version control proves itself better than DropBox or GoogleDrive, this ability to merge text files line-by-line and highlight the conflicts between them, so no work is ever lost.

### Remote repository and branching

    $ git checkout -b feature2

To push that branch to the remote repository:

    $ git push origin feature2

To list all branches (local and remote):

    $ git branch -a

To delete the remote branch:

    $ git push origin :feature2 


## The story so far...

So far, we've now seen how we can,

* Host our private repository on GitHub
* Copy, or clone, our remote repository onto a local machine
* Make changes in a local repository and push these to a remote repository
* Fetch and merge, or pull, changes from a remote repository into our local repository
* Identify and resolve conflicts when the same file is edited within two repositories
* Create and delete remote branches

For the next exercise, please pair up (or work in groups of 3) and work on a shared repository in GitHub.

Previous: [Tracking our changes with a local repository](1_Local.md) Next: [Exercise - collaborating](3_Collaboration.md)
