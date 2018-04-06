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

# Function to get CIs for certain HSP methods.
ci_95_states2values <- function(state_probs) {
  state_prob_cumsum <- t(apply(state_probs, 1, cumsum))
  ci_5 <- apply(state_prob_cumsum, 1, function(x) { min(which(x >= 0.05)) - 1 })
  ci_95 <- apply(state_prob_cumsum, 1, function(x) { min(which(x >= 0.95)) - 1 } )
  return(c(ci_5, ci_95))
}


# Function to identify traits that have only 1 state across tips (for pre-processing mk_model input).
identify_single_state <- function(state_vec) {
  state_vec_noNA <- na.exclude(state_vec)
  state_vec_noNA_l <- length(state_vec_noNA)
  
  if(length(which(state_vec_noNA == 1)) == state_vec_noNA_l) {
    return(TRUE)
  } else {
    return(FALSE)
  }
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
  
} else if(hsp_method == "mk_model" | hsp_method == "emp_prob" | hsp_method == "mp") {

  # Add 1 to all input counts because because traits states need to start at 1.
  trait_states_mapped <- trait_values_ordered + 1

  if(hsp_method == "mk_model") {
    
    # mk_model throws an error if there is only 1 known state (e.g. if the trait isn't found in any tips).
    # To get around this can just set all predictions to be this known state for these traits.
    single_state_traits <- which(mclapply(trait_states_mapped, identify_single_state, mc.cores=num_cores) == TRUE)
    
    if(length(single_state_traits) > 0) {
      trait_states_mapped_variable <- trait_states_mapped[-single_state_traits]
    } else {
      trait_states_mapped_variable <- trait_states_mapped
    }
    
    if(length(trait_states_mapped_variable) == 0) {
      stop("Stopping - all input traits have the same state value, which won't work with a Markov model.")
    }
    
    hsp_out_models <- mclapply(trait_states_mapped_variable,
                             hsp_mk_model,
                             tree = full_tree,
                             rate_model = "SUEDE",
                             check_input = check_input_set,
                             mc.cores = num_cores)
    
    # If applicable add the invariable traits into the hsp_mk_model output
    if(length(single_state_traits) > 0) {
      for(single_state_trait in names(single_state_traits)) {
        hsp_out_models[[single_state_trait]] <- list(likelihoods = matrix(rep(1, num_tip), nrow=num_tip, ncol=1), success = TRUE)
      }
    }
    
  } else if (hsp_method == "emp_prob") {

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

  
  # If calc_ci set then figure out what the assigned trait would be at the 95% CI and output resulting matrix.
  # Also output .rds file of list containing all state probabilities.
  if(calc_ci) {
    
    # Get subset of likelihood matrices for previously unknown tips only.
    unknown_tips_index <- which(full_tree$tip.label %in% unknown_tips)
    
    num_unknown <- length(unknown_tips)
    
    hsp_out_models_unknown_lik <- mclapply(names(hsp_out_models), 
                                           function(x) { get_sorted_prob(hsp_out_models[[x]]$likelihoods,
                                            study_tips_i=unknown_tips_index, 
                                            tree_tips=full_tree$tip.label, 
                                            study_tips=unknown_tips)},
                                            mc.cores = num_cores)
    
    names(hsp_out_models_unknown_lik) <- names(hsp_out_models)
    
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
    
    saveRDS(object = hsp_out_models_unknown_lik, file = "/home/gavin/tmp/hsp_out_models_unknown_lik.rds")
    
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
