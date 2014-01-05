# Make and Makefiles

- Make is an automatated build system, designed to avoid costly
  recomputation.

- *make* examines a *Makefile*, which contains a set of rules describing
  dependencies among files.

- A rule is run if the *target* is older than any of its *dependencies*.

- Older: compare creation time of files.

- Example:

        res.txt: param1.dat param2.dat
			simulation param1.dat param2.dat > res1.dat
			post-process res1.dat > res.txt

- Commands to be run should be indented with a TAB.  (Check your
  editor settings, esp R studio.)


## Makefile conventions

   - PHONY targets: denote actions; ignore filenames with same
     name. PHONY targets are always out of date, and so always run.

        |------------+---------------------------------------|
        | command    | action                                |
        |------------+---------------------------------------|
        | make       | check first rule                      |
        | make all   | rebuild everything                    |
        | make clean | remove files that can be rebuilt      |
        | touch file | update timestamp, preserving contents |
        |------------+---------------------------------------|


## References

[http://linuxdevcenter.com/pub/a/linux/2002/01/31/make_intro.html](http://linuxdevcenter.com/pub/a/linux/2002/01/31/make_intro.html)

[http://www.gnu.org/software/make/](http://www.gnu.org/software/make/)


