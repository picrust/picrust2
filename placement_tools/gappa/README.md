![gappa](/doc/logo/logo_readme.png?raw=true "gappa")

<!-- Genesis Applications for Phylogenetic Placement Analysis -->

<!-- [![Build Status](https://travis-ci.org/lczech/gappa.svg?branch=master)](https://travis-ci.org/lczech/gappa) -->
[![License](https://img.shields.io/badge/license-GPLv3-blue.svg)](http://www.gnu.org/licenses/gpl.html)
![Language](https://img.shields.io/badge/language-C%2B%2B11-lightgrey.svg)
<!-- ![Language](https://img.shields.io/badge/language-python-lightgrey.svg)-->

Features
-------------------

gappa is a collection of commands for working with phylogenetic data.
Its main focus are evolutionary placements of short environmental sequences on a reference phylogenetic tree.
Such data is typically produced by tools like [EPA-ng](https://github.com/Pbdas/epa-ng),
[RAxML-EPA](http://sco.h-its.org/exelixis/web/software/epa/index.html) or
[pplacer](http://matsen.fhcrc.org/pplacer/) and usually stored in
[jplace](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0031009) files.
<!-- It however also offers some commands for working with data such as sequences or trees. -->

Many commands in gappa are implementations of our novel methods.<!-- described here ... and here ... -->
At the same time, it offers some commands that are also implemented in the excellent
[guppy](http://matsen.github.io/pplacer/generated_rst/guppy.html) tool.
However, being written in C++, our gappa is much faster and needs less memory for most of the tasks.

Setup
-------------------

To run gappa on your machine, simply get it, and build it:

~~~.sh
git clone --recursive https://github.com/lczech/gappa.git
cd gappa
make
~~~

You can also use the green "Clone and download" button to get the source.
Still, call `make` in the main directory to build everything.

Requirements:

 *  [Make](https://www.gnu.org/software/make/) and [CMake](https://cmake.org/) 2.8.7 or higher.
 *  A fairly up-to-date C++11 compiler, e.g.,
    [clang++ 3.6](http://clang.llvm.org/) or [GCC 4.9](https://gcc.gnu.org/), or higher.

After building, the executable is stored in the `bin` directory, and used as follows.

Command Line Interface
-------------------

gappa is used via its command line interface, with subcommands for each task.
The subcommands have the general structure:

    gappa <module> <subcommand> <options>

<!-- The modules are simply a way of organizing the subcommands,
and have no [deeper meaning](https://en.wikipedia.org/wiki/42_%28answer%29). -->

For a list of all subcommands and their documentation,
see [the Wiki pages](https://github.com/lczech/gappa/wiki).

Behind the scenes
-------------------

gappa is short for Genesis Applications for Phylogenetic Placement Analysis.
This is because most of the work of gappa is actually performed by our [genesis](https://github.com/lczech/genesis) library.
See there if you are interested in the implementation details.
