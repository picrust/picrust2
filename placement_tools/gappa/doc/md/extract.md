## Description

The command extracts the queries that are placed in specified clades of the reference tree and writes per-clade `jplace` files. It is mainly intended for the multilevel placement approach as explained [here](https://www.biorxiv.org/content/early/2018/04/11/299792). It can of course also be used for other purposes where one is interested in just working with placements in a specific clade.

![Multilevel placement extraction workflow.](https://github.com/lczech/gappa/blob/master/doc/png/workflow_multilevel.png?raw=true)

The command takes one or more jplace files as input, as well as a file describing clades of the reference tree used in the jplace files.
It then finds all placements in those clades and writes per-clade placement files, each of them containing only those placements that had more of their mass (likelihood weights) in that clade than specified by `--threshold`.
Furthermore, two special clades are produced: `basal_branches`, which collects all placements that have their mass on branches that do not belong to any clade, as well as `uncertain`, which collects placements where no clade (including the basal clade) have more than the threshold amount of the mass in them.

Furthermore, if a set of fasta sequence files is provided, the command also creates per-clade fasta files, containing the sequences corresponding to the placements of the jplace files.
This of course necessitates that the sequences are named the same as the placements - which is given if the placement files are simply the result of placing the sequences on a reference tree.

## Details

The algorithm assigns a clade to each of the branches of the reference tree (either one of the specified ones, or the basal clade).
For terminal branches (leaves), the assigned clade is simply as specified in the `--clade-list-file`.
Inner branches are assigned to a clade if all leaves on one side of the split that is induced by the branch belong to the same clade.
In other words, all branches of a subtree that contains only taxa from one clade are assigned to that clade.
See the figure below for an example.

### `--clade-list-file`

This file describes which taxa of the reference tree are considered to belong to which clade.
Each line of the file needs to contain a taxon name of the tree, and the name of the clade it belongs to, separated by a tab:

    AF401522_Carchesium_polypinum	Alveolata
    X56165_Tetrahymena_thermophila	Alveolata
    X03772_Paramecium_tetraurelia	Alveolata
    ...

Not all taxa of the reference tree have to be part of the file; all missing ones are simply considered to be part of the special basal clade.

### `--color-tree-file`

If provided with an output file name, an svg file is written that shows which branches of the tree were assigned to which clade:

![Tree with clades marked in color.](https://github.com/lczech/gappa/blob/master/doc/png/prepare_extract.png?raw=true)

This is useful to verify the process and to make sure that the correct branches were selected.
In the figure, the basal branches are gray, while three exemplary clades are marked in color.

The behavior of selecting branches so that their subtrees are monophyletic with respect to a clade is visible here as well: For example, the green clade is split into two subtrees and a few single branches.
