## Description

The command takes a set of `jplace` files, and calculates and visualizes the Edge Dispersion per edge of the reference tree. The files need to have the same reference tree.

Edge Dispersion is explained and evaluated in detail in our article (in preparation). The following figure and its caption are an example adapted from this article:

<br>

![Dispersion Trees.](https://github.com/lczech/gappa/blob/master/doc/png/analyze_dispersion.png?raw=true)

<br>

The command is useful as a first exploratory tool to detect placement heterogeneity across samples. Subfigure (a) shows the standard deviation of the edge masses, without any further processing. One outlier (marked with an arrow) dominates the variances, which hides the values on most other edges. Thus, in subfigure (b), we used logarithmic scaling, which reveals more details on the edges with lower placement mass variance.
Subfigure (c) shows the Index of Dispersion of the edge masses, that is, the variance normalized by the mean. That means, edges with a higher number of placements can also have a higher variance. The subfigure again uses a logarithmic scale because of the outlier. The subfigure reveals more details on edges that exhibit a lower variance, which are shown in medium green colors. Lastly, subfigure (d) shows the variance of edge imbalances (instead of edge masses), and thus reveals information about whole clades of the tree.

## Details

By default, the command creates dispersion trees using all valid combinations of variants of the method. The following two options change this behavior.

### Edge Masses and Imbalances (`--edge-values`)

Controls whether to use masses or imbalances. By default, trees using both of them are crated. Using masses highlights the dispersion on single edges, while using imbalances considers whole clades. See the article for details on the differences between these two variants.

### Dispersion Method (`--method`)

Controls which method of dispersion is used for the visualization. By default, all valid ones are used, that is, trees for each of them are created.

When using edge masses (see `--edge-values`), the per-branch values can be scaled and normalized in different ways: Simple variance (`var`) or standard deviation (`sd`), coefficient of variation (`cv`, that is, standard deviation divided by mean), or variance to mean ratio (`vmr`, also called the Index of Dispersion), and the logarithmically scaled versions of these (`var-log`, `sd-log`, `cv-log`, and `vmr-log`).

When using edge imbalances however, only the variance and standard deviation are valid methods. This is because imbalances are not zero-based values, so dividing by mean is not a reasonable operation.

### Normalization (`--mass-norm`)

As the command is meant to show differences in a set of `jplace` samples files, it is important how those are normalized. Thus, the option is required.

If using `--mass-norm relative`, each sample (that is, each input `jplace` file) is normalized to unit mass 1.0, so that they all contribute equally to the result. Hence, the dispersion is measured relatively. That is, a branch exhibits a high dispersion if samples differ in the relative amount of placements on that branch (or in the clade, for imbalances) compared to the other placements in that sample.

On the other hand, if `--mass-norm absolute` is specified, the samples are not normalized. Thus, dispersion is measured absolutely. Branches then exhibit a high dispersion, if samples differ in the absolute number of placements on that branch (or clade). This can vastly differ from the normalized result, as the dispersion then depends on the total number of pqueries in each sample - which in turn depend on things like amplification bias, rarefaction, and other factors that can change the total number of sequences per sample.

The decision whether to use relative or absolute abundances depends on the use case and what each sample represents. See our article for details.

<!--
Example to run both:
${GAPPA} analyze dispersion --jplace-path ${SAMPLES} --write-svg-tree --svg-tree-ladderize --out-dir ${BASEDIR}/dispersion/ --tree-file-prefix disp_rel_ --mass-norm relative
${GAPPA} analyze dispersion --jplace-path ${SAMPLES} --write-svg-tree --svg-tree-ladderize --out-dir ${BASEDIR}/dispersion/ --tree-file-prefix disp_abs_ --mass-norm absolute

Caveat: imbalances basically only makes sense with relative masses... should mention that!
-->
