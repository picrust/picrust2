## Description

The command reverses the effects of the [chunkify](../wiki/Subcommand:-chunkify) command (see there for details on the workflow). That is, it takes the abundance map files and the per-chunk placement files as input, and creates a placement file for each of the original input sequences files, with all abundances and sequences names correctly restored. The command is thus one of the steps of our data preprocessing pipeline for phylogenetic placements as described [here](https://www.biorxiv.org/content/early/2018/04/11/299792).

The easiest way to input the placement files to the command is the `--jplace-path` option, which takes a list of files or a directory containing `.jplace` files. This option works in all cases, and can even handle cases where sequences were moved around between chunks, or chunks that were merged later, and so on. It simply uses the hash names of the sequences to identify them.

## Details

For large datasets, using the `--jplace-path` option might need too much memory, as all files have to be scanned for the sequence hash names first. This is necessary if the jplace files do not correspond exactly to the chunk files. However, if each jplace file was created from one chunk file, there is no need to scan for hashes in other files. Thus, we offer two memory- and time-saving alternatives:

### `--chunk-list-file`

The option takes a file, which needs to contain one jplace file path per line, in the order of the original chunks. For example, let's say the original sequence files were split into 13 chunks `chunk_0.fasta` to `chunk_12.fasta` by the `chunkify` command. Each of them was then placed on the reference tree, producing 13 jplace files. Then, the list file could look like this:

```
/path/to/chunk_0/result_0.jplace
/path/to/chunk_1/result_1.jplace
/path/to/chunk_2/result_2.jplace
/path/to/chunk_3/result_3.jplace
/path/to/chunk_4/result_4.jplace
/path/to/chunk_5/result_5.jplace
/path/to/chunk_6/result_6.jplace
/path/to/chunk_7/result_7.jplace
/path/to/chunk_8/result_8.jplace
/path/to/chunk_9/result_9.jplace
/path/to/chunk_10/result_10.jplace
/path/to/chunk_11/result_11.jplace
/path/to/chunk_12/result_12.jplace
```

That is, each line contains a path, in the original order of the chunks. Then, in order to create the placement entry for a sequence, the number `n` of the chunk in which the sequence was "chunkified" is used to find the correct jplace file by using the file in the `n`-th line of the list.

<!-- Theoretically, the same effect could be achieved by providing the files in the correct order to the `--jplace-path` option, but this is error-prone and less easy to trace later. Thus, we did not implement this way. -->

### `--chunk-file-expression`

Alternatively, if the naming of the per-chunk jplace files is as straight forward as above, that is, the file names are just numbered, it is also possible to use an expression instead of the list file, where the `@` character is used as a placeholder for the number:

    --chunk-file-expression /path/to/chunk_@/result_@.jplace

This has the same effect as using the list file.
