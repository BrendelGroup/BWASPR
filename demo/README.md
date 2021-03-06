# BWASPR Demo

Analyses in this directory are based on the BWASPR package, available on github:

	https://github.com/brendelgroup/BWASPR

As documented on the github page, BWASPR processes (BWASP-generated) *.mcalls
files.
Input specification is via two files: the *.dat data file and the *.par
parameter file, which are read by the BWASPR function setup_BWASPR().

Although there are of course several ways of running R-based workflows, the data
studies in this directory were conveniently produced by Rscript.
The workflow is represented in the file

	Rscript.BWASPR

which coulde be invoked as follows (putting output into OUT-arg1):

	Rscript --default-packages=methods,datasets,utils,grDevices,graphics,stats  Rscript.BWASPR arg1 OUT-arg1

where arg1 is the prefix of a configuration file arg1.conf.
The command will source(arg1.conf) within Rscript.BWASPR, then run the workflow
according to the specifications in arg1.conf.

For convenience and testing, use

	./xcheck arg1

which will call the Rscript command and put output in a directory NEW-arg1 (and
compare with an existing output directory OUT-arg1, if it exists, a useful check
on updated code).

What will be run is determined by the logicals RUNsomething in the *.conf file.
The default

	RUNload    <- FALSE

	RUNcms     <- TRUE
	RUNpwc     <- TRUE
	RUNcrl     <- TRUE

	RUNrepcms  <- TRUE
	RUNrepcrl  <- TRUE

	RUNmmp     <- TRUE
	RUNacs     <- TRUE
	RUNrnk     <- TRUE
	RUNmrpr    <- TRUE

	RUNdmt     <- TRUE
	RUNdmsg    <- TRUE
	RUNdmgdtls <- TRUE
	RUNogl     <- TRUE

	RUNsave    <- FALSE

would run all implemented analyses from scratch.
To omit certain steps (e.g., to re-run only some of the analyses with changed
parameters), change TRUE to FALSE.
RUNload=TRUE will load a previously saved arg1.RData file into the R workspace
before executing other analyses steps.
To generate an arg1.RData file from the current run, set RUNsave=TRUE.

To check the mechanics of running BWASPR, you can use the small data set in
[../inst/exdata](../inst/extdata).
Run

	./xdoit sample OUT-sample

and take a look at the 0README files in the output directory and its
subdirectories.

If you decided to forgo the R package installation on your machine but have
[Singularity](https://www.sylabs.io/docs/) installed, then the following will do
nicely:

	singularity pull --name bwaspr.simg shub://BrendelGroup/BWASPR
	singularity exec -e -B `pwd` bwaspr.simg ./xcheck sample

(using _xcheck_ here in case you ran _xdoit_ and want to compare local versus
singularity output).

## From sample to real data:
This demo is designed to show you the mechanics of the workflow on a small
data set.
So, how do you get a real data set?
Back to [BWASP](https://github.com/brendelgroup/BWASP), which tells you how to
generate the required _*.mcalls_ and _GFF3DIR_ files.
