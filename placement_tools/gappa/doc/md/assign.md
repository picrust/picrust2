## Description

The command takes one or more `jplace` files as input and assigns the likelihood weights of each placement to a taxonomic rank, then prints a high level profile of how the total likelihood weight is distributed on the taxonomy (specified by the `--taxon-file`). To achieve this, the command operates in three phases.

First, the tree found in the `jplace` input is labelled according to the information found in the [`taxon-file`](#taxonomic-label-file--taxon-file), beginning at the tip nodes.
Inner nodes of the tree are labelled according to the most specific taxonomic rank shared by both its child nodes (i.e. if one child is labelled "A;B;C", and the other is labelled "A;B;D", the parent node will be labelled "A;B"). The resulting labelled tree is printed to file as the first intermediate result (filename `labelled_tree`).

Second, the algorithm goes through each placement in the `jplace` input and assigns its likelihood weight to one or more taxonomic ranks, according to the [specified strategy](#distribution-ratio---distribution-ratio). This triggers the second intermediate result, which is a file containing the per-query results of this assignment (filename `per_query_assign`).

Third, the command summarizes these assignments by collapsing them into one tabulation, showing information about the total distribution of likelihood weight across the taxonomy ([example](#final-output), filename `profile.csv`).

## Details

### jplace File(s) (`--jplace-path`)

List of jplace files or directories to process. For directories, only files with the extension `.jplace` are processed. When multiple files are specified, the command treats them all as one collective input (as opposed to processing each independently).
All files must have compatible trees, i.e. same topology and tip labels that are congruent with the assignment in the `taxon-file`.
We further strongly recommend that also the branch lengths are identical to facilitate better comparability, however differences in branch lengths will not cause the command to fail.

### Taxonomic Label File (`--taxon-file`)

This file is used to assign machine-readable taxonomic paths to taxa (tips) of the reference tree (which is taken from the `jplace` input).
The format is as follows. Each line assigns a taxonomic path to one taxon, and contains two columns: the taxon label as it appears in the tree, followed by the semicolon-separated taxonomic path. The two columns are separated by a tab character.

#### Example:

```
Seal    Eukaryota;Animalia;Chordata;Mammalia;Carnivora;Phocoiae
Whale   Eukaryota;Animalia;Chordata;Mammalia;Cetartiodactyla;
Mouse   Eukaryota;Animalia;Chordata;Mammalia;Rodentia;Muridae
Human   Eukaryota;Animalia;Chordata;Mammalia;Primates;Homonidae
Chicken Eukaryota;Animalia;Chordata;Amphibia;Galliformes;Phasianidae
Frog    Eukaryota;Animalia;Chordata;Amphibia;Anura;Dendrobatidae
Loach   Eukaryota;Animalia;Chordata;Amphibia;Anura;Rhacophoridae
Cow     Eukaryota;Animalia;Chordata;Mammalia;Artiodactyla;Bovidae
```

### Distribution Ratio (`--distribution-ratio`)

This option controls the strategy by which the likelihood weight of a placement is assigned to the taxonomic labels associated with the placement branch.
If the option is omitted, the strategy is to use an automatic ratio, meaning that the ratio is calculated from the attachment point of the placement on the placement branch, as specified by the `distal` (or `proximal`) length field of the placement.

#### Automatic Ratio Example:

In this example, the edge drawn as the line between `p` and `d` represent the placement branch in the reference tree, while `q` represents the attached query sequence. Where the line from `q` meets the reference branch is the attachment point, defined by the `distal length`. As in this case the distal length is about 1/3rd of the total branch length, the automatic distribution ratio will result in 2/3rds of the likelihood weight associated with this placement to contribute to the taxonomic label of the node `d`, while the label associated with node `p` receives 1/3rd of the likelihood weight.

```
          q
          |
          |
p---------------d
           <--->
              |
         distal length
```

#### Fixed Ratio Example:

When `--distribution-ratio` is used with a fixed specified value (between 0.0 and 1.0), that value defines how the likelihood weight is split between the two labels.
It defines the fraction of the likelihood weight that will contribute to the label at node `p` in the example below, where `p` is the _proximal node_, meaning the node that is closer to the (virtual) root of the tree.

```
--distribution-ratio 0.25

      q
      |
      |
p-----------d
▲ 0.25      ▲ 0.75
 \         /
  \       /
   \     /
LWR: 1.0
```

This is for example useful to produce assignments similar to [Sativa](https://github.com/amkozlov/sativa), which uses a fixed ratio of 0.49 (in order to break ties).
We however usually recommend using the default automatic ratio (that is, do not specify the `--distribution-ratio` options), as this takes more phylogenetic information into account.

### Final Output

The final output of the command is a tabulation of the distribution of the total likelihood weight across the taxonomy, which is written to the file `profile.csv`.

The meaning of the column headers are:

 - `LWR`: likelihood weight that was assigned to this exact taxonomic path
 - `fract`: `LWR` divided by the global total likelihood weight
 - `aLWR`: accumulated likelihood weights that were assigned either to this taxonomic path or any taxonomic path below this
 - `afract`: `aLWR` divided by the global total likelihood weight
 - `taxopath`: the taxonomic path

#### Example

```
LWR     fract   aLWR    afract  taxopath
0       0       2       1       Eukaryota
0       0       2       1       Eukaryota;Animalia
0.49    0.245   2       1       Eukaryota;Animalia;Chordata
0       0       1       0.5     Eukaryota;Animalia;Chordata;Amphibia
0.49    0.245   1       0.5     Eukaryota;Animalia;Chordata;Amphibia;Anura
0.51    0.255   0.51    0.255   Eukaryota;Animalia;Chordata;Amphibia;Anura;Rhacophoridae
0       0       0.51    0.255   Eukaryota;Animalia;Chordata;Mammalia
0       0       0.51    0.255   Eukaryota;Animalia;Chordata;Mammalia;Rodentia
0.51    0.255   0.51    0.255   Eukaryota;Animalia;Chordata;Mammalia;Rodentia;Muridae
```

### Filtering (`--sub-taxopath`)

Additionally, the tabulated output may be filtered (constrained) to include only a part of the taxonomy.
This behavior is regulated via the `--sub-taxopath` option, and the result is printed to a file called `profile_filtered.csv` in the specified output directory.

For this output, the normalization of the accumulated LWR and fraction columns only takes the specified subtaxonomy into account.
This option is hence useful to get a more detailed view at a specific part of the taxonomy.

#### Example `--sub-taxopath "Eukaryota;Animalia;Chordata;Amphibia"`:

```
LWR     fract   aLWR    afract  taxopath
0.49    0.49    1       1       Anura
0.51    0.51    0.51    0.51    Anura;Rhacophoridae
```
