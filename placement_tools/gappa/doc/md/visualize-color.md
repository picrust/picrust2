## Description

The command takes one or more jplace files as input and visualizes the distribution of placements on the branches of the tree. It uses color coding to show how much placement mass there is per branch.

![Placements visualized by per-branch colors.](https://github.com/lczech/gappa/blob/master/doc/png/analyze_visualize_color.png?raw=true)

**Important remark:**
If multiple jplace files are provided as input, their combined placements are visualized. It is then critical to correctly set the `--mass-norm` option. If set to `absolute`, no normalization is performed per jplace file - thus, absolute abundances are shown. However, if set to `relative`, the placement mass in each input file is normalized to unit mass 1.0 first, thus showing relative abundances.
