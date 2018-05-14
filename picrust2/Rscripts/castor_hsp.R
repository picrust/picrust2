#!/usr/bin/env Rscript

# Read in quietly to avoid outputting "Loading required package: Rcpp" to stderr.
library(castor, quietly = TRUE)

# Read in ape as well for "read.tree" function.
library(ape)

# Load parallel package to run over multiple cores.
library(parallel)

Args <- commandArgs(TRUE)

# Read in command-line arguments.
full_tree <- read.tree(Args[1])
trait_values <- read.delim(Args[2], check.names=FALSE, row.names=1)
hsp_method <- Args[3]
calc_nsti <- as.logical(Args[4])
calc_ci <- as.logical(Args[5])
check_input_set <- as.logical(Args[6])
num_cores <- Args[7]
predict_outfile <- Args[8]
ci_outfile <- Args[9]
rds_outfile <- Args[10]
seed_setting <- Args[11]
write_rds <- as.logical(Args[12])

# Set random seed if integer specified.
if(seed_setting != "None") {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(as.integer(seed_setting))
}


# Function to get CIs for certain HSP methods.
ci_95_states2values <- function(state_probs) {
  
  if(ncol(state_probs) > 1) {
    state_prob_cumsum <- t(apply(state_probs, 1, cumsum))
  } else {
    state_prob_cumsum <- state_probs
  }
  
  ci_5 <- apply(state_prob_cumsum, 1, function(x) { min(which(x >= 0.05)) })
  ci_95 <- apply(state_prob_cumsum, 1, function(x) { min(which(x >= 0.95)) })
  
  return(c(ci_5, ci_95))
}


# Function to get HSP state probabilities for study (i.e. "unknown" tips only).
# Adds rownames of sequences and colnames of counts. 
# Also remove columns that are all zeros (no probability of that state).
get_sorted_prob <- function(in_likelihood, study_tips_i, tree_tips) {
  
  # Subet to study sequences only and set as rownames.
  tmp_lik <- in_likelihood[study_tips_i, , drop=FALSE]
  rownames(tmp_lik) <- tree_tips[study_tips_i]
  
  # Set column names to be 0 to max num of counts.
  colnames(tmp_lik) <- c(0:(ncol(tmp_lik)-1))
  
  # Remove columns that are 0 across all sequences.
  col2remove <- which(colSums(tmp_lik) == 0)
  if(length(col2remove) > 0) {
    tmp_lik <- tmp_lik[, -col2remove, drop=FALSE]
  }
  
  return(tmp_lik)
  
}


# Order the trait table to match the tree tip labels. Set all tips without a value to be NA.
unknown_tips_index <- which(! full_tree$tip.label %in% rownames(trait_values))
unknown_tips <- full_tree$tip.label[unknown_tips_index]

unknown_df <- as.data.frame(matrix(NA,
                                   nrow=length(unknown_tips),
                                   ncol=ncol(trait_values)))

rownames(unknown_df) = unknown_tips
colnames(unknown_df) = colnames(trait_values)

trait_values_ordered <- rbind(trait_values, unknown_df)

trait_values_ordered <- trait_values_ordered[full_tree$tip.label, , drop=FALSE]

num_tip <- nrow(trait_values_ordered)

if (hsp_method == "pic" | hsp_method == "scp" | hsp_method == "subtree_average") {
  
  if (hsp_method == "pic") {
    predict_out <- mclapply(trait_values_ordered,
                            hsp_independent_contrasts,
                            tree=full_tree,
                            weighted=TRUE,
                            check_input=check_input_set,
                            mc.cores = num_cores)
    
  } else if (hsp_method == "scp") {
    
    predict_out <- mclapply(trait_values_ordered,
                            hsp_squared_change_parsimony,
                            tree=full_tree,
                            weighted=TRUE,
                            check_input=check_input_set,
                            mc.cores = num_cores)
    
  } else if (hsp_method == "subtree_average") {
    
    predict_out <- mclapply(trait_values_ordered,
                            hsp_subtree_averaging,
                            tree = full_tree,
                            check_input = check_input_set,
                            mc.cores = num_cores)
  }
  
  predicted_values <- mclapply(predict_out, function(x) { x$states[unknown_tips_index] }, mc.cores = num_cores)
  
} else if(hsp_method == "emp_prob" | hsp_method == "mp") {
  
  # Add 1 to all input counts because because traits states need to start at 1.
  trait_states_mapped <- trait_values_ordered + 1
  
  if (hsp_method == "emp_prob") {
    
    hsp_out_models <- mclapply(trait_states_mapped,
                               hsp_empirical_probabilities,
                               tree = full_tree,
                               check_input = check_input_set,
                               mc.cores = num_cores)
    
  } else if (hsp_method == "mp") {
    
    hsp_out_models <- mclapply(trait_states_mapped,
                               hsp_max_parsimony,
                               tree = full_tree,
                               check_input = check_input_set,
                               transition_costs = "proportional",
                               weight_by_scenarios = TRUE,
                               mc.cores = num_cores)
    
  }
  
  # Get subset of likelihood matrices for previously unknown tips only and output RDS file.
  num_unknown <- length(unknown_tips)
  
  hsp_out_models_unknown_lik <- mclapply(names(hsp_out_models), 
                                         function(x) { get_sorted_prob(hsp_out_models[[x]]$likelihoods,
                                                                       study_tips_i=unknown_tips_index, 
                                                                       tree_tips=full_tree$tip.label)},
                                         mc.cores = num_cores)
  
  names(hsp_out_models_unknown_lik) <- names(hsp_out_models)
  
  # Get state with highest probability in each case.
  predicted_values <- mclapply(hsp_out_models_unknown_lik,
                               function(x) { colnames(x)[max.col(x)] },
                               mc.cores = num_cores)
  
  # Save RDS object if option set.  
  if(write_rds) {
    saveRDS(object = hsp_out_models_unknown_lik, file = rds_outfile)
  }
  
  # If calc_ci set then figure out what the assigned trait would be at the 95% CI and output resulting matrix.
  if(calc_ci) {
    
    ci_values <- data.frame(mclapply(hsp_out_models_unknown_lik,
                                     function(x) { ci_95_states2values(x) },
                                     mc.cores = num_cores),
                            check.names = FALSE)
    
    colnames(ci_values) <- names(hsp_out_models_unknown_lik)
    
    ci_values_ci_5 <- ci_values[1:num_unknown, , drop=FALSE]
    ci_values_ci_95 <- ci_values[(num_unknown+1):(num_unknown*2), , drop=FALSE]
    
    colnames(ci_values_ci_5) <- paste(colnames(ci_values_ci_5), "5", sep="_")
    colnames(ci_values_ci_95) <- paste(colnames(ci_values_ci_95), "95", sep="_")
    
    ci_values <- cbind(ci_values_ci_5, ci_values_ci_95)
    
    # Sort column names so that 5% and 95% CIs are next to each other.
    ci_values <- ci_values[ , order(names(ci_values))]
    
    orig_ci_colnames <- colnames(ci_values)    
    ci_values$sequence <- unknown_tips
    ci_values <- ci_values[, c("sequence", orig_ci_colnames)]
    
    write.table(ci_values, file=ci_outfile, sep="\t", quote=FALSE, row.names=FALSE)
  }
  
}

# Add "sequence" as first column of predicted_values.
predicted_values <- data.frame(predicted_values, check.names = FALSE)
predicted_values$sequence <- full_tree$tip.label[unknown_tips_index]
predicted_values <- predicted_values[, c("sequence", colnames(trait_values_ordered))]
# Calculate NSTI per tip and add to output as last column if option set.
if(calc_nsti) {
  predicted_values$metadata_NSTI <- NA
  
  # Calculate NSTIs for tips with previously unknown trait values.
  all_tip_range <- 1:length(full_tree$tip.label)
  
  known_tip_range <- which(! full_tree$tip.label %in% unknown_tips)
  
  predicted_values[, "metadata_NSTI"] <- find_nearest_tips(full_tree,
                                                           target_tips=known_tip_range,
                                                           check_input=check_input_set)$nearest_distance_per_tip[unknown_tips_index]
}

# Write out predicted values.
write.table(predicted_values, file=predict_outfile, row.names=FALSE, quote=FALSE, sep="\t")
