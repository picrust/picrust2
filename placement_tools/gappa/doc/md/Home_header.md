gappa is a collection of commands for working with phylogenetic data.
Its main focus are evolutionary placements of short environmental sequences on a reference phylogenetic tree.
See [Phylogenetic Placement](../wiki/Phylogenetic-Placement) for an introduction describing a typical pipeline.
<!-- It however also offers some commands for working with data such as sequences or trees. -->

Many commands in gappa are implementations of our novel methods.
At the same time, it offers some commands that are also implemented in the excellent
[guppy](http://matsen.github.io/pplacer/generated_rst/guppy.html) tool.
However, being written in C++, our gappa is much faster and needs less memory for most of the tasks.

## Command Line Interface

gappa is used via its command line interface, with subcommands for each task.
The subcommands have the general structure:

    gappa <module> <subcommand> <options>

The modules are simply a way of organizing the subcommands.
