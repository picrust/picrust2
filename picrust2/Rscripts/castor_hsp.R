#!/usr/bin/env Rscript

# Read in quietly to avoid outputting "Loading required package: Rcpp" to stderr.
library(castor, quietly = TRUE)

Args <- commandArgs(TRUE)

# Read in command-line arguments.
full_tree <- read_tree(file=Args[1], check_label_uniqueness = TRUE)
trait_values <- read.delim(Args[2], check.names=FALSE, row.names=1)
hsp_method <- Args[3]
edge_exponent_set <- as.numeric(Args[4])
calc_ci <- as.logical(Args[5])
check_input_set <- as.logical(Args[6])
predict_outfile <- Args[7]
ci_outfile <- Args[8]
seed_setting <- Args[9]

# Set random seed if integer specified.
if(seed_setting != "None") {
    set.seed(as.integer(seed_setting))
}

# Function to get CIs for certain HSP methods.
ci_95_states2values <- function(state_probs) {
  
  if(ncol(state_probs) > 1) {
    state_prob_cumsum <- t(apply(state_probs, 1, cumsum))
  } else {
    state_prob_cumsum <- state_probs
  }
  
  ci_5 <- apply(state_prob_cumsum, 1, function(x) {  as.numeric(colnames(state_probs)[min(which(x >= 0.05))]) })
  ci_95 <- apply(state_prob_cumsum, 1, function(x) {  as.numeric(colnames(state_probs)[min(which(x >= 0.95))]) })
  
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


### Function to wrap hsp_max_parsimony and return probabilities for 
### study sequences only (and for non-zero states only).
mp_study_probs <- function(in_trait, in_tree ,unknown_i, check_input) {

  mp_hsp_out <- hsp_max_parsimony(tree = in_tree,
                                  tip_states = in_trait,
                                  check_input=check_input,
                                  transition_costs = "proportional",
                                  edge_exponent=edge_exponent_set,
                                  weight_by_scenarios = TRUE)
  
  return(get_sorted_prob(mp_hsp_out$likelihoods,
                         study_tips_i=unknown_i, 
                         tree_tips=in_tree$tip.label))
}


### Function to wrap hsp_empirical_probabilities and return probabilities for 
### study sequences only (and for non-zero states only).
emp_prob_study_probs <- function(in_trait, in_tree, unknown_i, check_input) {
  
  emp_prob_hsp_out <- hsp_empirical_probabilities(tree = in_tree,
                                                  tip_states = in_trait,
                                                  check_input=check_input)
  
  return(get_sorted_prob(emp_prob_hsp_out$likelihoods,
                         study_tips_i=unknown_i, 
                         tree_tips=in_tree$tip.label))
}


# First check if any tip labels are present in the trait table but do not match
# because they have quotes around them.
if(length(grep("\"", full_tree$tip.label) > 0) || length(grep("\'", full_tree$tip.label) > 0)) {

    tmp_unknown_tips_index <- which(! full_tree$tip.label %in% rownames(trait_values))
    tmp_unknown_tips <- full_tree$tip.label[tmp_unknown_tips_index]

    unknown_labels_no_quotes <- gsub("\'", "", tmp_unknown_tips)
    unknown_labels_no_quotes <- gsub("\"", "", unknown_labels_no_quotes)

    no_quote_matches = which(unknown_labels_no_quotes %in% rownames(trait_values))

    # Remove quotes from around tips where matches only occur when quotes are removed.
    if(length(no_quote_matches) > 0) {
        indices_to_change <- tmp_unknown_tips_index[no_quote_matches]
        full_tree$tip.label[indices_to_change] <- unknown_labels_no_quotes[no_quote_matches]
    }

    cat(paste("\nWhile running castor_hsp.R: Fixed ", length(no_quote_matches),
        " tip label(s) that had quotations around them that stopped them from matching the function abundance table.",
        sep=""))
}

# Order the trait table to match the tree tip labels. Set all tips without a value to be NA.
unknown_tips_index <- which(! full_tree$tip.label %in% rownames(trait_values))
unknown_tips <- full_tree$tip.label[unknown_tips_index]
num_unknown <- length(unknown_tips)
num_known <- length(full_tree$tip.label) - num_unknown

# Throw error if all tips are unknown.
if(num_unknown == length(full_tree$tip.label)) {
  stop("None of the reference ids within the function abundance table are found within the input tree. This can occur when malformed or mismatched custom reference files are used.")

# Check if only very few tips are "known".
} else if(num_known / nrow(trait_values) < 0.1) {
    stop("Fewer than 10% of reference ids in the function abundance table are in the tree. This could be because the ids are slightly different between the table and the tree.")
}

unknown_df <- as.data.frame(matrix(NA,
                                   nrow=num_unknown,
                                   ncol=ncol(trait_values)))

rownames(unknown_df) = unknown_tips
colnames(unknown_df) = colnames(trait_values)

# Get combined dataframe with known and unknown tips
# (unknown tips have NA as trait values).
trait_values <- rbind(trait_values, unknown_df)

# Remove unknown_df object from memory.
remove(unknown_df)
invisible(gc(verbose = FALSE))

# Order this combined trait table by the order of tips in the tree.
trait_values <- trait_values[full_tree$tip.label, , drop=FALSE]

num_tip <- nrow(trait_values)

if (hsp_method == "pic" || hsp_method == "scp" || hsp_method == "subtree_average") {
  
  if (hsp_method == "pic") {
    predict_out <- lapply(trait_values,
                            hsp_independent_contrasts,
                            tree=full_tree,
                            weighted=TRUE,
                            check_input=check_input_set)
    
  } else if (hsp_method == "scp") {
    
    predict_out <- lapply(trait_values,
                            hsp_squared_change_parsimony,
                            tree=full_tree,
                            weighted=TRUE,
                            check_input=check_input_set)
    
  } else if (hsp_method == "subtree_average") {
    
    predict_out <- lapply(trait_values,
                            hsp_subtree_averaging,
                            tree = full_tree,
                            check_input = check_input_set)
  }
  
  predicted_values <- lapply(predict_out, function(x) { x$states[unknown_tips_index] })
  
  # Remove raw predict_out object from memory.
  remove(predict_out)
  invisible(gc(verbose = FALSE))
  
} else if(hsp_method == "emp_prob" || hsp_method == "mp") {
  
  # Add 1 to all input counts because because traits states need to start at 1.
  trait_values <- trait_values + 1
  
  if (hsp_method == "emp_prob") {
    
    hsp_out_models_unknown_lik <- lapply(trait_values, 
                                         function(x) {
                                           emp_prob_study_probs(in_trait = x,
                                                                in_tree = full_tree,
                                                                unknown_i = unknown_tips_index,
                                                                check_input = check_input_set)})
  } else if (hsp_method == "mp") {
    
    hsp_out_models_unknown_lik <- lapply(trait_values, 
                                         function(x) {
                                           mp_study_probs(in_trait = x,
                                                          in_tree = full_tree,
                                                          unknown_i = unknown_tips_index,
                                                          check_input = check_input_set)})
  }

  # Get state with highest probability in each case.
  predicted_values <- lapply(hsp_out_models_unknown_lik,
                               function(x) { as.numeric(colnames(x)[max.col(x)]) })

  # If calc_ci set then figure out what the assigned trait would be at the 95% CI and output resulting matrix.
  if(calc_ci) {
    
    ci_values <- data.frame(lapply(hsp_out_models_unknown_lik,
                                     function(x) { ci_95_states2values(x) }),
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
predicted_values$sequence <- unknown_tips
predicted_values <- predicted_values[, c("sequence", colnames(trait_values))]

# Check to see if there are any columns that are totally missing.
missing_values_by_column <- colSums(is.na(predicted_values))

if(length(which(missing_values_by_column == nrow(predicted_values))) > 0) {
    stop("\nError - at least one trait in the prediction table was entirely missing values.")
} else if(length(which(missing_values_by_column > 0)) > 0) {
    cat("\nWarning: there are missing values in the output prediction table.")
}

# Write out predicted values.
write.table(predicted_values, file=predict_outfile, row.names=FALSE, quote=FALSE, sep="\t")
