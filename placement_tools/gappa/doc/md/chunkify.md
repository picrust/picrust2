## Description

The command is one of the steps of our data preprocessing pipeline for phylogenetic placements as described [here](https://www.biorxiv.org/content/early/2018/04/11/299792). It takes one or more fasta files as input, e.g., each representing an environmental sample. It then writes out numbered chunks files of equal size, containing the unique sequences of the input. For each input file, it also writes an abundance map file, which stores the per-sequence abundances in the input. In order to identify unique sequences, it uses a hash value of the sequence data, which is also assigned as a new name to the sequences in the chunks.

![Chunkify and Unchunkify Workflow.](https://github.com/lczech/gappa/blob/master/doc/png/workflow_chunkify.png?raw=true)

The produced chunk files are intended to be used with [phylogenetic placement](../wiki/Phylogenetic-Placement) next (after potentially aligning them first to the reference). Using chunks of equal size ensures relatively stable run times for each chunk, so that large datasets can be processed efficiently on a computer cluster. Furthermore, as the chunks only contain unique sequences, compute time is further reduced.

After finishing phylogentic placement, the [unchunkify](../wiki/Subcommand:-unchunkify) command then takes the per-chunk placement files as well as the abundance map files produced here, and creates placement files for each of the original input files, with all abundances and original sequences names restored. Thus, the combination of these two commands achieves the same effect as placing each input file separately, but lowers computational cost and maximizes load balancing.

## Details

The memory usage depends on the number of unique sequences in the input data.
For example, we used a test dataset with 1,170 fasta files (31.5 GB)
containing 182,556,655 sequences, thereof 104,947,033 unique.
Using `--hash-function SHA1` and `--threads 1`,
the program ran for 40 min and used 4 GB of memory on a laptop computer.
The run time can be reduced by using multiple threads,
at the cost of slightly more memory usage.
As most of the sequences only appeared with an abundance of 1 per file,
we also tried the option `--min-abundance 2`.
In that case, 82% of the sequences were filtered due to low abundance,
resulting in a run time of 12 min and 440 MB memory usage with a single thread.

### `--min-abundance`

The `--min-abundance` option is meant for fasta files which were already deduplicated, that is, which do not contain duplicate sequences, but have unique sequences with their abundances annotated in the sequence label itself.
Such files are for example produced by the [vsearch](https://github.com/torognes/vsearch) command `--derep_fulllength`.
If the `--min-abundance` option is used with a value greater than 1 for files that do not contain abundance information in the sequence labels, all sequences will be filtered out, because their abundance is then considered to be 1.

The abundance information in sequence labels can be specified in two ways:

  * Using the format `[;]size=123[;]`, appearing anywhere in the label. The semicoli are optional.
  * Appended via underscore: `name_123`. In this case, the number has to be the last in
    the label, that is, no other text may follow.
