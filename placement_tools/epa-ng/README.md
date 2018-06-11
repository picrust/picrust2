[![Gitter chat](https://badges.gitter.im/Pbdas/epa-ng.png)](https://gitter.im/epa-ng/Lobby)

# EPA-ng

1. **[Introduction](#introduction)**
2. **[Build Instructions](#build-instructions)**
3. **[Usage](#usage)**
4. **[Citing EPA-ng](#citing-epa-ng)**

## IMPORTANT
The by far easiest way to clone this repository is to use the following command
```
git clone --recursive https://github.com/Pbdas/epa.git
```

If that is not an option (perhaps you downloaded the zip file), you can fix the submodules by running

```
git submodule update --init --recursive
```

### SUPPORT

The most reliable way to get in touch with us is to head over to the [Phylogenetic Placement Google Group](https://groups.google.com/forum/#!forum/phylogenetic-placement). You can also search its history, or the hostory of the  [RAxML Google Group](https://groups.google.com/forum/#!forum/raxml) for your particular question.

Alternatively I've created a [gitter chat room](https://gitter.im/epa-ng/Lobby) where I can usually be found during office hours.

### DISCLAIMER

This tool is still in an active, *beta* phase of development. Suggestions, bug reports and constructive comments are more than encuraged! Please do so in the [google group](https://groups.google.com/forum/#!categories/phylogenetic-placement/epa-ng).

## Introduction

`EPA-ng` is a complete rewrite of the [Evolutionary Placement Algorithm (EPA)](https://academic.oup.com/sysbio/article/60/3/291/1667010/Performance-Accuracy-and-Web-Server-for), previously implemented in [RAxML](https://github.com/stamatak/standard-RAxML). It uses [libpll](https://github.com/xflouris/libpll-2) and [pll-modules](https://github.com/ddarriba/pll-modules) to perform maximum likelihood-based phylogenetic placement of genetic sequences on a user-supplied reference tree and alignment.

### What can EPA-ng do?

- do phylogenetic placement using the **GTR+GAMMA ML model only** (for now)
- take as input **separated reference and query alignment files**, in the **fasta** format (for now)
- handle **DNA** and **Amino Acid** data
- distributed computing suitable for the **cluster**
- **prepare inputs** for the cluster:
  - convert query fasta file into a random access, binary encoded file called a `bfast`-file
- output the placement results in the [jplace format](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0031009) ready for downstream analysis by frameworks such as [genesis](https://github.com/lczech/genesis)

## Build Instructions

After cloning with the above command, ensure the following packages are installed or otherwise available (relevant modules loaded on your cluster):

```
sudo apt-get install autotools-dev libtool flex bison cmake automake autoconf
```

Once these dependencies are available, you need to ensure that your compiler is recent enough, as EPA-ng is built using C++14 features. The minimum required versions are as follows:

| Compiler | Min. Version |
| - | - |
| gcc | 4.9.2  |
| clang  | 3.8  |
| icc | 16 |

Any one of these compilers will be sufficient. gcc is the most wide spread, and current versions of Ubuntu have gcc versions exceeding the minimum.

Now it's time to build the program.

```
make
```

Thats it! The executable should now be located in the `epa-ng/bin/` folder.

### Apple

Supported in theory, but currently work in progress. In principle same procedure as under Linux.

### Windows

No support yet, will gaugue interest first.

## Usage

`EPA-ng` is used from the command line, as the main use-case is processing large amounts of data using a supercomputing cluster.

Here is a list of the most basic arguments you will use:

| Flag | Long Flag | Meaning |
| - | - | - |
| -s | --ref-msa | reference MSA (fasta)  |
| -t | --tree | reference Tree (newick)  |
| -q | --query | query sequences (fasta or [bfast](#converting-the-query-file)) |
| -w | --outdir | output directory (default: current directory) |
|  | --model | [model parameter specification](#setting-the-model-parameters) |
| -T | --threads | number of threads to use |

For a full overview of command line options either run `EPA-ng` with no input, or with the flag `-h` (or `--help`).

### Basic

On a single computer, an example execution might look like this:

```
epa-ng --ref-msa $REF_MSA --tree $TREE --query $QRY_MSA --outdir $OUT
```

Note that this will use as many threads as specified by the environment variable `OMP_NUM_THREADS`.
Usually this defaults to the number of cores.
Note however, that no speedup is to be expected from hyperthreads, meaning the number of threads should be set to the number of physical cores.

#### Setting the Model Parameters
As of version 0.2.0, GTRGAMMA model parameters have to be specified explicitly.
There are currently two ways of doing this:
Either specify a raxml-ng-style model descriptor, like so:
```
epa-ng <...> --model GTR{0.7/1.8/1.2/0.6/3.0/1.0}+FU{0.25/0.23/0.30/0.22}+G4{0.47}
```

... or pass a RAxML 8.x RAxML_info-file to the program, where the info file was generated from a call to RAxML option `-f e`:
```
raxmlHPC-AVX -f e -s $REF_MSA -t $TREE -n file -m GTRGAMMAX
(generates a RAxML_info file with model parameters, called RAxML_info.file)
epa-ng <...> --model RAxML_info.file
```


### Advanced

Overview of advanced features:

| Flag | Long Flag | Meaning |
| - | - | - |
| -g | --dyn-heur | use dynamic [preplacement heuristic](#configuring-the-heuristic-preplacement) (default)  |
| -G | --fix-heur | use fixed [preplacement heuristic](#configuring-the-heuristic-preplacement) |
|  | --no-heur | disable [preplacement heuristic](#configuring-the-heuristic-preplacement) |
|  | --no-pre-mask | disable [premasking](#premasking) |
| -c | --bfast | [convert query fasta to binary format](#converting-the-query-file) |

The description of basic cluster usage starts [here](#running-on-the-cluster)

#### Configuring the Heuristic Preplacement

By default, `EPA-ng` performs placement of a sequence in two stages: first selecting promising branches quickly (preplacement), then evaluating the selected branches in greater detail.

`EPA-ng` currently offers three ways of selecting these candidates.

The default is the *accumulated threshold* method, in which branches are added to the set of candidates until the sum of their LWR exceed a user specified threshold.
The flag controlling this mode is `-g` (or `--dyn-heur`), with a default setting of `0.99999`, corresponding to a covered likelihood weight of 99.999%.

The second mode functions identically to the candidate selection mode in the original implementation of the `EPA` in `RAxML`.
Here again the branches are sorted by the LWR of the placement of a sequence.
Then, the top x% of the total number of branches are selected into the set of candidates.
Like in `RAxML`, this behavior is controlled via the `-G` (or `--fix-heur`) flag.

The third mode works identically to the *baseball heuristic* from [pplacer](http://matsen.github.io/pplacer/), with default settings (--strike-box 3.0, --max-strikes 6, --max-pitches 40) and is enabled using the `--baseball-heur` flag.

Lastly, to disable the preplacement completely, you can simply supply the `--no-heur` flag.
Be warned however: doing so will be significantly more computationally demanding.
Our advice is to use the heuristic, as it sacrifices only insignificant amounts of accuracy for greatly improved speed.

#### Premasking

By default, `EPA-ng` enables *premasking*, which works similarily to the same option in [pplacer](http://matsen.github.io/pplacer/):
If a site of the alignment is all-gaps in either the reference OR query alignment, throw it out.
Further, for each query sequence, ignore the leading, and trailing gap columns (this is where we differ from pplacer, as they ignore ALL query gap columns).

This reduces both runtime and memory footprint greatly, depending on the data.
For short read data, the impact will be massive, as typically query alignments will be mostly all-gap.

### Cluster usage

Before using the cluster version of `EPA-ng`, the input files must be preprocessed.
In general, these preprocessed files also enable more streamlined re-runs of experiments and are reccomended for use, always.

#### Converting the query file

You may also explicitly convert the input query fasta file to our internal fasta format.
This format is binary encoded (reducing the size by half) and randomly accessible.
Again, using this format is highly reccomended, and required to use the MPI version.

To convert the fasta file, simply run the program with the query file specified thusly:

```
epa-ng --bfast query.fasta --outdir $OUT
```

This will produce a file called `query.fasta.bfast` in the specified output directory.

#### Running on the cluster

To use distributed parallelism in `EPA-ng`, first we must re-compile the program with MPI enabled.
This requires a version of MPI to be loaded/installed on your system.
The only additional requirement `EPA-ng` has, is that the compiler that is loaded in conjunction with MPI satisfies the [minimum version requirements](#build-instructions).
Often this can be assured by the order in which the relevant modules are loaded on the cluster: first MPI, then the compiler.
However we reccomend you contact your support team should this cause issues for you.

The actual compilation is very straight-forward:
```
make clean && make EPA_HYBRID=1
```

This will attempt to compile the program with both MPI and OpenMP, as the most efficient way to run the program is to map one MPI rank per node (good alternative: one rank per socket!), each rank starting as many threads as there are *physical* cores.

In your job submission script, you can then call the program in a highly similar way to before:

```
mpirun epa-ng --ref-msa $REF_MSA --tree $TREE -q query.fasta.bin -w ./some/output/dir
```

**Note** that using the binary format is not strictly required, however it is highly reccomended to increase parallel efficiency.


## Citing EPA-ng

For now, please cite [this preprint](https://www.biorxiv.org/content/early/2018/04/17/291658) of the paper when using EPA-ng.
