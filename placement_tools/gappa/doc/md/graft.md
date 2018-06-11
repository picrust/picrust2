## Description

The command takes the reference tree of the provided placefile(s), and for each pquery, it attaches a new leaf node to the tree, positioned according to its proximal length and pendant length. The resulting tree is useful to get an overview of the distribution of placements. It is mainly intended to view a few placements. For large samples, it might be a bit cluttered.

Similar trees are produced by [RAxML-EPA](http://sco.h-its.org/exelixis/web/software/epa/index.html), where the file is called `RAxML_labelledTree`, and by the [`guppy tog` command](http://matsen.github.io/pplacer/generated_rst/guppy_tog.html#guppy-tog). Both programs differ in the exact way the the placements are added as edges. To control this behaviour, use the `--fully-resolve` parameter.

## Details

The provided `jplace` files are processed individually, producing a `newick` tree for each of them.
They are named like the input files, but replace the file extension by `.newick`.

### Without `--fully-resolve`

If `--fully-resolve` is not provided (default), all placements at one edge are collected as children of one central base edge:

![Multifurcating grafted tree.](https://github.com/lczech/genesis/blob/master/doc/img/placement/labelled_tree_multifurcating.png?raw=true)

This method is similar to the way [RAxML-EPA](http://sco.h-its.org/exelixis/web/software/epa/index.html) produces a grafted tree, which is there called "labelled tree".

The base edge is positioned on the original edge at the average `proximal_length` of the placements. The base edge has a multifurcation if there are more than two placements on the edge.

The pendant length of the placements is used to calculate the branch length of the new placement edges. This calculation subtracts the shortest pendant length of the placements on the edge, so that the base edge is maximally "moved" towards the placement edges. This also implies that at least one of the placement edges has branch length == 0.0. Furthermore, the placements are sorted by their pendant length.

Using this method, the new nodes of the resulting tree are easier to distinguish and collapse, as all placements are collected as children branching off from the base edge. However, this comes at the cost of losing the detailled information of the proximal length of the placements. If you want to keep this information, use `--fully-resolve` instead.

### With `--fully-resolve`

If `--fully-resolve` is provided, all placements per branch are turned into single leaf nodes:

![Fully resolved grafted tree.](https://github.com/lczech/genesis/blob/master/doc/img/placement/labelled_tree_fully_resolved.png?raw=true)

This method is similar to the way [`guppy tog`]((http://matsen.github.io/pplacer/generated_rst/guppy_tog.html#guppy-tog)) produces a grafted tree.

The original edge is splitted into separate parts where each placement edge is attached. The branch lengths between those parts are calculated using the proximal length of the placements, while the branch lengths of the placement edges use their pendant length.

Using this method gives maximum information, but results in a more crowded tree. The new placement edges are "sorted" along the original edge by their proximal length. For this reason in the example image above, "Query 2" is closer to "Node A" then "Query 1": it has a higher proximal length. This information was lost in the multifurcating tree shown before.

### Further Details

For edges that contain only a single placement (or none at all), both versions of fully resolve behave the same. In this case, the placement is simply attached using its proximal length and pendant length.

Pqueries with multiple names are treated as if each name is a separate placement, i.e., for each of them, a new (identical) edge is added to the Tree. If using `--fully-resolve`, this results in a branch length of 0.0 between the nodes of those placements.

### `--name-prefix`

Specify a prefix to be added to all new leaf nodes (the ones that represent placements). This is useful if a pquery name also occurs as a name in the original tree. By default, empty. In order to get the same naming as grafted trees as produced by RAxML, use `QUERY___`.
