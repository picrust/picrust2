#!/usr/bin/env Rscript
#./castor_asr_pic.R  <input: reference_tree_with_internal_nodes_labelled>
#                    <input: functional_matrix_file>
#                    <output: count_table>
#                    <output_CI_flag: boolean>
#                    <output: prob_table>

# Read in quietly to avoid outputting "Loading required package: Rcpp" to stderr.
library(castor, quietly = TRUE)

# Read in ape as well for "read.tree" function.
library(ape)

Args <- commandArgs(TRUE)

# Read in command-line arguments.
pruned_tree <- read.tree(Args[1])
trait_values <- read.delim(Args[2], check.names=FALSE, row.names=1)
count_out_file <- Args[3]
ci_set <- as.logical(Args[4])

# If ci_set == TRUE then make sure that Args[5] exists.
if(ci_set == TRUE) {
  if(is.na(Args[5])) {
    stop("If CI intervals set then need to specify an outfile!")
  }
}

ci_out_file <- Args[5]

#If tips in tree contain '_' then read.tree() places single quotes around these tip labels.
#This then causes sorting errors below since the rownames are different between the trait table and the tree.
#Fix this by putting quotes around any labels in the trait table that have a '_'.
for(i in grep("_",rownames(trait_values))) {
 rownames(trait_values)[i] <- paste("'", rownames(trait_values)[i], "'", sep="")
}

# Order the trait table to match the tree tip labels.
trait_values_ordered <- trait_values[pruned_tree$tip.label , , drop=FALSE]

# Run asr_independent_contrasts on each trait.
reconstructions <- apply(trait_values_ordered, 2, asr_independent_contrasts, tree=pruned_tree, weighted=TRUE, check_input=FALSE, include_CI=ci_set)

# Pull out only the node predictions.
just_ancestral_states <- lapply(1:length(reconstructions), function(x) reconstructions[[x]]$ancestral_states)
names(just_ancestral_states) <- names(reconstructions)

# Reformat the list into a matrix.
ancestral_states_matrix <- do.call(cbind, just_ancestral_states)

# Relabel the node names (ones created internally by castor) with the actual node labels in the tree.
ancestral_states_matrix <- cbind(pruned_tree$node.label, ancestral_states_matrix)

# Set first column name to be "nodes".
colnames(ancestral_states_matrix)[1] <- 'nodes'

# Write ancestral node predictions to file.
write.table(data.frame(ancestral_states_matrix, check.names=FALSE), file=count_out_file, row.names=FALSE, quote=FALSE, sep="\t")

if(ci_set) {
  # Get lower and upper 95% CI based on output radius.
  ci <- lapply(1:length(reconstructions),
               function(x) paste(round(reconstructions[[x]]$ancestral_states - reconstructions[[x]]$CI95, digits=4), 
                                 round(reconstructions[[x]]$ancestral_states + reconstructions[[x]]$CI95, digits=4),
                                 sep="|")
  )
  
  names(ci) <- names(reconstructions)                                                                                
  ci_matrix <- do.call(cbind, ci)
  ci_matrix <- cbind(pruned_tree$node.label, ci_matrix)
  
  # Set column name to be "nodes".
  colnames(ci_matrix)[1] <- 'nodes'

  # Output the data to file.
  write.table(data.frame(ci_matrix, check.names=FALSE), file=ci_out_file, row.names=FALSE, quote=FALSE, sep="\t")
}
