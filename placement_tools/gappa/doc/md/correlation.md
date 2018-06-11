## Description

The command takes a set of `jplace` files (called samples), as well as a table containing metadata features for each sample. It then calculates and visualizes the Edge Correlation with the metadata features per edge of the reference tree. The files need to have the same reference tree.

Edge Correlation is explained and evaluated in detail in our article (in preparation). The following figure and its caption are an example adapted from this article:

<br>

![Correlation Trees.](https://github.com/lczech/gappa/blob/master/doc/png/analyze_correlation.png?raw=true)

<br>

All subfigures show red edges or red paths at the clade on the left hand side of the tree. This indicates that presence of placements in this clade is anti-correlated with the used metadata feature. On the other hand, blue and green edges, which indicate positive correlation, are spread throughout the tree the same way in all subfigures. The extent of correlation is larger for Spearman’s Coefficient, indicating that the correlation is
monotonic, but not strictly linear.

## Details

By default, the command creates correlation trees for all valid metadata features, using all variants of the method. In the following, we first explain how to specify the metadata, and then how to change the default behavior.

### Metadata Features (`--metadata-file`)

The metadata features are specified in a comma separated table file (`.csv`). The first row needs to contain the feature names, which are used as file names for the output files. The first column needs to contain the file names of the `jplace` files (samples) without extension.

Example:

```
File,Temperature,Salinity Sensor,Oxygen Sensor
ERR562588,19.85,36.32,221.47
ERR562558,23.83,37.49,n/a
ERR562591,26.23,36.62,199.94
ERR562643,21.44,37.89,207.79
ERR562637,26.64,35.36,189.81
```

This table specifies three types of metadata for five files `ERR562588.jplace`, `ERR562558.jplace`, etc. Note the `n/a` value in the last column. Any non-numerical value is interpreted as missing data, and is simply left out when calculating the correlation. That is, the last column only uses four data points.

### Features Selection (`--metadata-fields`)

When specifying a comma-separated list of column headers of the meatadata table, only these features are used. Otherwise, all numerical columns are used, and trees for all for all of them are created.

Example: In order to only use the first two features of the above table, specify `--metadata-fields "Temperature,Salinity Sensor"` with the command. Note the double quotes, which are necessary here, as one of the feature names contains a space.

### Edge Masses and Imbalances (`--edge-values`)

Controls whether to use masses or imbalances. By default, trees using both of them are crated. Using masses highlights the correlation of single edges, while using imbalances considers whole clades. See the article for details on the differences between these two variants.

### Correlation Method (`--method`)

Controls which method of correlation is used for the visualization. By default, Pearsons and Spearmans are used, that is, trees for each of them are created.

### Normalization (`--mass-norm`)

As the command is meant to show differences in a set of `jplace` samples files, it is important how those are normalized. Thus, the option is required.

If using `--mass-norm relative`, each sample (that is, each input `jplace` file) is normalized to unit mass 1.0, so that they all contribute equally to the result. Hence, the correlation is measured relatively. That is, a branch exhibits a high correlation with a metadata feature depending on the relative amount of placements on that branch (or in the clade, for imbalances) compared to the other placements in that sample.

On the other hand, if `--mass-norm absolute` is specified, the samples are not normalized. Thus, correlation is measured absolutely. Branches then exhibit a high correlation (or anti-correlation) with a metadata feature depending on the absolute number of placements on that branch (or clade). This can vastly differ from the normalized result, as the values then depends on the total number of pqueries in each sample - which in turn depend on things like amplification bias, rarefaction, and other factors that can change the total number of sequences per sample.

The decision whether to use relative or absolute abundances depends on the use case and what each sample represents. See our article for details.
