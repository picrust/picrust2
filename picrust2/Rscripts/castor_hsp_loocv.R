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
tips_to_test <- as.vector(read.delim(Args[3], header=FALSE, stringsAsFactors = FALSE)$V1)
hsp_method <- Args[4]
expect_outfile <- Args[5]
predict_outfile <- Args[6]
metrics_outfile <- Args[7]
num_cores <- as.numeric(Args[8])

trait_values_ordered <- trait_values[full_tree$tip.label, , drop=FALSE]

num_tip <- nrow(trait_values_ordered)

# Restrict tips to test to those found in tree.
tips_to_test <- tips_to_test[which(tips_to_test %in% full_tree$tip.label)]

# Prep predicted out df.
predict_out <- data.frame(matrix(NA, nrow=length(tips_to_test), ncol=ncol(trait_values_ordered)))

rownames(predict_out) <- tips_to_test
colnames(predict_out) <- colnames(trait_values_ordered)

# Prep metric outfile.
metrics_out <- data.frame(matrix(NA, nrow=length(tips_to_test), ncol=3))
rownames(metrics_out) <- tips_to_test
colnames(metrics_out) <- c("rmse", "rho", "nsti")

# Loop through all tips to leave out.
for(tip_name in tips_to_test) {
  
  tip_index <- which(full_tree$tip.label == tip_name)
  
  if (hsp_method == "mp") {

    trait_states_mapped <- as.data.frame(trait_values_ordered, drop=FALSE) + 1

    trait_states_mapped[tip_index,] <- NA
    
    hsp_out_models <- mclapply(trait_states_mapped,
                               hsp_max_parsimony,
                               tree = full_tree,
                               transition_costs = "proportional",
                               weight_by_scenarios = TRUE,
                               mc.cores = num_cores)
  }
  
  predict_out[tip_name,] <- unlist(mclapply(hsp_out_models,
                                            function(x) { max.col(x$likelihoods[tip_index, , drop=FALSE]) - 1 },
                                            mc.cores = num_cores))
  
  # Calculate root mean squared deviation, Spearman's rho, and nearest sequenced taxon index for tip label.
  tip_predict <- as.numeric(predict_out[tip_name,])
  tip_expect <- as.numeric(trait_values_ordered[tip_name,])
  
  rmse_out <- sqrt(mean((tip_predict - tip_expect)^2))
  
  if(length(tip_predict) > 1) {
    rho_out <- as.numeric(cor.test(tip_predict, tip_expect)$estimate)
  } else {
   rho_out <- NA 
  }
  
  nsti_out <- find_nearest_tips(full_tree, target_tips=full_tree$tip.label[-tip_index])$nearest_distance_per_tip[tip_index]
  
  metrics_out[tip_name,] <- c(rmse_out, rho_out, nsti_out)
  
}


# Write output tables.
write.table(data.frame("assembly"=rownames(predict_out), predict_out, check.names=FALSE), file=predict_outfile, row.names=FALSE, quote=FALSE, sep="\t")

expect_out <- trait_values_ordered[tips_to_test, , drop=FALSE]
write.table(data.frame("assembly"=rownames(expect_out), expect_out, check.names=FALSE), file=expect_outfile, row.names=FALSE, quote=FALSE, sep="\t")

write.table(metrics_out, file=metrics_outfile, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")
