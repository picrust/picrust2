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

# Set random seed if integer specified.
if(seed_setting != "None") {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(as.integer(seed_setting))
}

# Function to get CIs for certain HSP methods.
ci_95_states2values <- function(state_probs, number_of_tips) {
  state_prob_cumsum <- t(apply(state_probs[1:number_of_tips, , drop=FALSE], 1, cumsum))
  ci_5 <- apply(state_prob_cumsum, 1, function(x) { min(which(x >= 0.05)) - 1 })
  ci_95 <- apply(state_prob_cumsum, 1, function(x) { min(which(x >= 0.95)) - 1 } )
  return(c(ci_5, ci_95))
}


# Function to get HSP state probabilities for study (i.e. "unknown" tips only).
get_sorted_prob <- function(in_likelihood, study_tips_i, tree_tips, study_tips) {
  tmp_lik <- in_likelihood[unknown_tips_index, , drop=FALSE]
  
  rownames(tmp_lik) <- tree_tips[study_tips_i]
  
  return(tmp_lik[study_tips, ,drop=FALSE])
}


# Order the trait table to match the tree tip labels. Set all tips without a value to be NA.
unknown_tips <- full_tree$tip.label[which(! full_tree$tip.label %in% rownames(trait_values))]

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

  predicted_values <- mclapply(predict_out, function(x) { x$states[1:num_tip] }, mc.cores = num_cores)
  
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
  
  # Get state with highest probability in each case and subtract 1 to get value.
  predicted_values <- mclapply(hsp_out_models,
                               function(x) { max.col(x$likelihoods[1:num_tip,]) - 1 },
                               mc.cores = num_cores)

  # Get subset of likelihood matrices for previously unknown tips only and output RDS file.
  unknown_tips_index <- which(full_tree$tip.label %in% unknown_tips)
  
  num_unknown <- length(unknown_tips)
  
  hsp_out_models_unknown_lik <- mclapply(names(hsp_out_models), 
                                         function(x) { get_sorted_prob(hsp_out_models[[x]]$likelihoods,
                                                                       study_tips_i=unknown_tips_index, 
                                                                       tree_tips=full_tree$tip.label, 
                                                                       study_tips=unknown_tips)},
                                                                       mc.cores = num_cores)
  names(hsp_out_models_unknown_lik) <- names(hsp_out_models)
  
  saveRDS(object = hsp_out_models_unknown_lik, file = rds_outfile)
  
  # If calc_ci set then figure out what the assigned trait would be at the 95% CI and output resulting matrix.

  if(calc_ci) {
    
    ci_values <- data.frame(mclapply(hsp_out_models_unknown_lik,
                                     function(x) { ci_95_states2values(x) },
                                     mc.cores = num_cores),
                            check.names = FALSE)
    
    colnames(ci_values) <- names(hsp_out_models)
    
    ci_values_ci_5 <- ci_values[1:num_unknown, , drop=FALSE]
    ci_values_ci_95 <- ci_values[(num_unknown+1):(num_unknown*2), , drop=FALSE]
    
    colnames(ci_values_ci_5) <- paste(colnames(ci_values_ci_5), "5", sep="_")
    colnames(ci_values_ci_95) <- paste(colnames(ci_values_ci_95), "95", sep="_")
    
    ci_values <- cbind(ci_values_ci_5, ci_values_ci_95)
    
    # Sort column names so that 5% and 95% CIs are next to each other.
    ci_values <- ci_values[ , order(names(ci_values))]
    
    orig_ci_colnames <- colnames(ci_values)
    
    ci_values$tips <- unknown_tips
    ci_values <- ci_values[, c("tips", orig_ci_colnames)]
    
    write.table(ci_values, file=ci_outfile, sep="\t", quote=FALSE, row.names=FALSE)
  }

}

# Add tips as first column of predicted_values.
predicted_values <- data.frame(predicted_values, check.names = FALSE)
predicted_values$tips <- full_tree$tip.label
predicted_values <- predicted_values[, c("tips", colnames(trait_values_ordered))]

unknown_tip_range <- which(predicted_values$tips %in% unknown_tips)

# Calculate NSTI per tip and add to output as last column if option set.
if(calc_nsti) {
  predicted_values$metadata_NSTI <- NA
  
  # Calculate NSTIs for tips with previously unknown trait values.
  all_tip_range <- 1:length(full_tree$tip.label)
  
  known_tip_range <- which(! predicted_values$tips %in% unknown_tips)

  predicted_values[unknown_tip_range, "metadata_NSTI"] <- find_nearest_tips(
                                                              full_tree,
                                                              target_tips=known_tip_range,
                                                              check_input=check_input_set)$nearest_distance_per_tip[unknown_tip_range]
}

# Subset to previously unknown tips only.
predicted_values <- predicted_values[unknown_tip_range,]

# Write out predicted values.
write.table(predicted_values, file=predict_outfile, row.names=FALSE, quote=FALSE, sep="\t")
