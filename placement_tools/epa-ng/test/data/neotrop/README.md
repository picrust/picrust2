## Description

This directory contains exemplary data from the *neotrop* data set mentioned in the paper. It was obtained from [1].

Setting up the directory:
```
tar xzvf neotrop.tar.gz
```
This dataset contains a reference tree (`tree.newick`) and reference alignment (`reference.fasta`).

Two query files are provided: one with 10,000 aligned query sequences (`query_10k.fasta`, 46MB) and one with 100,000 aligned query sequences (`query_100k.fasta`, 457MB).

Additinally I have provided a file called `infofile`. This is one way of passing model parameters to EPA-ng, as these parameters are _not_ re-inferred in the program itself.
To create this file I ran RAxML-8 with the `-f e` option which computes and outputs the model parameters in the `RAxML_info.*` file.

## Example

```
mkdir result
../../../bin/epa-ng --tree tree.newick --ref-msa reference.fasta --query query_10k.fasta --model infofile --outdir result
```

## References

[1] Mahé, F., de Vargas, C., Bass, D., Czech, L., Stamatakis, A., Lara, E., Singer, D., Mayor, J., Bunge, J., Sernaker, S. and Siemensmeyer, T., 2017. Parasites dominate hyperdiverse soil protist communities in Neotropical rainforests. Nature ecology & evolution, 1(4), p.0091.
