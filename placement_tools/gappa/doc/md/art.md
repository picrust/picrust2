## Description

Given a set of sequences and a fitting taxonomy, the command produces consensus sequences representing taxonomic clades, according to our ART method as described [here](https://www.biorxiv.org/content/early/2018/04/11/299792).
The main inputs are `--sequence-file` and `--taxonomy-file`, which provide the input data, as well as the `--target-size` of how many consensus sequences to build.

![ART workflow.](https://github.com/lczech/gappa/blob/master/doc/png/workflow_art.png?raw=true)

After running the command, the resulting set of sequences can be used to infer a reference tree using any tree inference program.

## Details

### `--taxonomy-file`

The taxonomy file needs to contain a list of the taxa used for the taxonomic expansion algorithm. Each line of the file lists a semicolon-separated taxonomic clade. Everything after the first tab is ignored.

Example:

```
Eukaryota;	4	domain		
Eukaryota;Amoebozoa;	4052	kingdom		119
Eukaryota;Amoebozoa;Myxogastria;	4094	phylum		119
Eukaryota;Amoebozoa;Myxogastria;Amaurochaete;	4095	genus		119
Eukaryota;Amoebozoa;Myxogastria;Badhamia;	4096	genus		119
Eukaryota;Amoebozoa;Myxogastria;Brefeldia;	4097	genus		119
Eukaryota;Amoebozoa;Myxogastria;Comatricha;	4098	genus		119
...
```

### `--sequence-file`

The sequence file needs to be in fasta format, and contain sequences that are labelled with the taxonomic path that they belong to. This taxonomic path can either be the whole label, or everything after the first whitespace (space or tab). This allows to have sequences with unique identifiers as the first part of the label.

For example, sequences in the [Silva](https://www.arb-silva.de/) database are labelled like this:

    >AY842031.1.1855 Eukaryota;Amoebozoa;Myxogastria;Amaurochaete;Amaurochaete comata
    >JQ031957.1.4380 Eukaryota;Amoebozoa;Myxogastria;Brefeldia;Brefeldia maxima

In the example, the sequences first contain a unique identifier, followed by a space and the taxonomic path the sequence belongs to. The path contains an additional taxonomic level which is not present in the database. If this occurs, the last level is assumed to be species level, and removed from the path. The resulting taxonommic path is part of the taxonomy, and hence the sequence can be used.

### `--sub-taxonomy`

If provided with a semicolon-separated taxonomic path (e.g., `Eukaryota;Amoebozoa;`), only this subclade is used for the algorithm. That is, the algorithm behaves as if the `--taxonomy-file` and `--sequence-file` only contained the taxa and sequences of the provided clade.
