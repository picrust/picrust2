#!/usr/bin/env Rscript

### This script will calculate nearest sequenced taxon index (NSTI) between study 
### sequences in a tree and the nearest reference sequence.

# Read in quietly to avoid outputting "Loading required package: Rcpp" to stderr.
library(castor, quietly = TRUE)

Args <- commandArgs(TRUE)

# Read in command-line arguments.
full_tree <- read_tree(file=Args[1], check_label_uniqueness = TRUE)

# Read in a file with each known tip, 1 per line.
known_tips <- read.table(Args[2], header=FALSE, stringsAsFactors = FALSE)$V1

# Read in path to output file.
output_path <- Args[3]

# Determine which tips are unknown.
unknown_tips_index <- which(! full_tree$tip.label %in% known_tips)
unknown_tips <- full_tree$tip.label[unknown_tips_index]
all_tip_range <- 1:length(full_tree$tip.label)
known_tip_range <- which(! full_tree$tip.label %in% unknown_tips)

nsti_values <- find_nearest_tips(full_tree,
                                 target_tips=known_tip_range,
                                 check_input=TRUE)$nearest_distance_per_tip[unknown_tips_index]
nsti_genomes <- find_nearest_tips(full_tree,
                                 target_tips=known_tip_range,
                                 check_input=TRUE)$nearest_tip_per_tip[unknown_tips_index]
print(known_tips, nsti_genomes)
nsti_genomes = known_tips[nsti_genomes]

# Make dataframe of study sequences (unknown tips) and nsti values as 2nd column.
write.table(x = data.frame("sequence" = unknown_tips, "metadata_NSTI" = nsti_values, "closest_reference_genome" = nsti_genomes),
            file = output_path,
            sep="\t",
            quote = FALSE,
            row.names=FALSE)
