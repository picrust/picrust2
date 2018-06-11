## Description

The command runs a Phylogenetic k-means clustering on a set of `jplace` files (called samples). The aim is to group samples that are similar to each other regarding the Phylogenetic KR distance. This is for example useful to find structure in a set of samples from different locations or points in time.

## Details

### Values for `k`

It is often not obvious what the "natural" number of clusters of a set of samples is. To this end, it makes sense to try different values for `k` and explore how the clustering changes. Then, techniques like the Elbow Method can be used to estimate a reasonable number of clusters. See below for more on that.

To this end, the option `--k` accepts multiple values, separated by commas, as well as ranges of numbers, specified via a dash. This is similar to how specific pages can be selected in common software before printing.

Example: `--k 1-6,10,15`

### Output Format

For each specified `k`, the result of the clustering is written to an assignment table, which lists for each sample the cluster number it was grouped into, as we as the distance (Phylogenetic KR distance) from the sample to the centroid of the cluster. The cluster numbers are zero based, and thus span the range `[0, k-1]`.

### Centroid Trees

If furthermore an output tree format is specified (via one of the `---write-...-tree` options), the centroids of each cluster are visualized as mass trees. That is, the average mass distribution of all samples that were assigned to a cluster is calculated and visualized on the tree. This is useful to explore what each cluster represents - that is, how the samples were clustered.

### Multiple `k` and Overview File

If multiple values for `k` are specified (see above), the option `--write-overview-file` can be used to write an overview table that lists for each value of `k` the average distance and variance from each sample to its assigned cluster centroid. This table can directly be visualized to create plots such as the Elbow Method.
