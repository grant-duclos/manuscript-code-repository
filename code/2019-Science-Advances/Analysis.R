# load necessary packages
suppressPackageStartupMessages(lapply(c("Biobase", "topicmodels", "MASS", "lsa"), library, character.only = TRUE))

# counts.mat: counts matrix
# type: "gene" or "cells" for gene set or cell cluster identification
# n.states: number of gene sets or cell clusters
# n.starts: number of random starts for LDA
# seed: seed for random number generation used during LDA
# max.q: maximum FDR corrected p-value (< 0.05 is recommended)
# min.coef: minimum regression coefficient (> 1 is recommended)
# min.spec: minimum state specificity (value between 0 and 1)
# min.sim: minimum state similaririty (value between 0 and 1)

run_analysis <- function(
	counts.mat = counts.mat,
	type = c("genes", "cells"),
	n.states = 3,
	n.starts = 5,
	seed = 12345,
	max.q = 0.05,
	min.coef = 1,
	min.spec = 0.1,
	min.sim = 0) {

	# Assign seeds for LDA
	seeds <- seed + 0:(n.starts - 1)

	# Parameters for LDA
	lda.param <- list(estimate.alpha = TRUE, estimate.beta = TRUE,
		verbose = 100, save = 0, keep = 0,
		seed = seeds, nstart = n.starts, best = TRUE)

	# Run LDA
	cat("Begin State Discovery \n")
	if (type == "genes") {
		result.lda <- LDA(t(counts.mat), k=n.states, method="VEM", control=lda.param)
		} else if (type == "cells") {
			result.lda <- LDA(counts.mat, k=n.states, method="VEM", control=lda.param)
			}
	cat("State Discovery Complete \n")

	# Extract States Matrix from LDA output
	states <- posterior(result.lda)$topics

	# Negative Binomial Regression
	result.stats <- data.frame(row.names=colnames(counts.mat))
	result.stats[unlist(lapply(1:ncol(states), paste, c("p", "FDR q", "Coef", "Spec", "Sim")))] <- NA

	# Maximum number of model iterations
	max.it <- 100
	cat("Evaluate State-Feature Association \n")

	if (type == "genes") {
		result.stats <- data.frame(row.names=rownames(counts.mat))
		result.stats[unlist(lapply(1:ncol(states), paste, c("p", "FDR q", "Coef", "Spec", "Sim")))] <- NA
		for (k in 1:ncol(states)) {
			cat("State", k, sep=" ", "\n")
			for (i in 1:nrow(counts.mat)) {
			    if ((i %% 100) == 0) cat("State", k, i, sep=" ", "\r")
			    fit <- try(glm.nb(counts.mat[i,] ~ states[,k], maxit=max.it), silent=TRUE)
			    if (!inherits(fit, "try-error")) {
			        result.stats[i, paste(k, "p")] <- summary(fit)$coefficients[2,4]
			        result.stats[i, paste(k, "Coef")] <- summary(fit)$coefficients[2,1]
			    	}
				}
			}
		} else if (type == "cells") {
			result.stats <- data.frame(row.names=colnames(counts.mat))
			result.stats[unlist(lapply(1:ncol(states), paste, c("p", "FDR q", "Coef", "Spec", "Sim")))] <- NA
			for (k in 1:ncol(states)) {
				cat("State", k, sep=" ", "\n")
				for (i in 1:ncol(counts.mat)) {
		    		if ((i %% 100) == 0) cat("State", k, i, sep=" ", "\r")
		    		fit <- try(glm.nb(counts.mat[,i] ~ states[,k], maxit=max.it), silent=TRUE)
		    		if (!inherits(fit, "try-error")) {
		        		result.stats[i, paste(k, "p")] <- summary(fit)$coefficients[2,4]
		        		result.stats[i, paste(k, "Coef")] <- summary(fit)$coefficients[2,1]
		    		}
				}
			}
		}

	cat("Evaluation Complete \n")

	# FDR correction
	for (prefix in sub(" p$", "", grep(" p$", colnames(result.stats), value=TRUE))) {
	    result.stats[[paste(prefix, "FDR q")]] <- p.adjust(result.stats[[paste(prefix, "p")]], method="fdr")
	}

	# Calculate State-Specificity and add to results matrix
	result.spec <- t(sweep(posterior(result.lda)$terms, 2, apply(posterior(result.lda)$terms, 2, sum), "/"))
	for (k in 1:ncol(states)) {
		result.stats[, paste(k, "Spec")] <- result.spec[,k]
	}
	cat("Calculate State-Specificity \n")

	# Calculate relative expression values: counts divided by total counts per cell
	rel.exp.mat <- sweep(counts.mat, MARGIN=2, STATS=colSums(counts.mat), FUN="/")

	# Calculate State-Similarity and add to results matrix
	if (type == "genes") {
		for (i in 1:nrow(rel.exp.mat)) {
			for (k in 1:ncol(states)) {
				result.stats[i, paste(k, "Sim")] <- cosine(rel.exp.mat[i,], states[,k])
				}
			}
		} else if (type == "cells") {
			for (i in 1:ncol(rel.exp.mat)) {
				for (k in 1:ncol(states)) {
					result.stats[i, paste(k, "Sim")] <- cosine(rel.exp.mat[,i], states[,k])
				}
			}
		}
	cat("Calculate State-Similarity \n")

	# Feature Selection
	result.stats["Assignment"] <- NA

	# Extract FDR q-values and model coefficients from statistical analysis
	for (i in 1:nrow(result.stats)) {
		target.q <- result.stats[i,][grep("FDR",names(result.stats[i,]),value=T)]
		target.coef <- result.stats[i,][grep("Coef",names(result.stats[i,]),value=T)]
		target.spec <- result.stats[i,][grep("Spec",names(result.stats[i,]),value=T)]
		target.sim <- result.stats[i,][grep("Sim",names(result.stats[i,]),value=T)]

		target.sig <- Reduce(intersect, list(
			gsub("[^0-9]", "", names(target.q)[which(target.q < max.q)]),
			gsub("[^0-9]", "", names(target.coef)[which(target.coef > min.coef)]),
			gsub("[^0-9]", "", names(target.spec)[which(target.spec > min.spec)]),
			gsub("[^0-9]", "", names(target.sim)[which(target.sim > min.sim)])
			))

		#
		if (length(target.sig) == 0) {
			result.stats[i, "Assignment"] <- NA
			} else if (length(target.sig) == 1) {
				result.stats[i, "Assignment"] <- target.sig
				} else {
					result.stats[i, "Assignment"] <- gsub("[^0-9]", "", names(result.stats[i,][paste(target.sig, "FDR q")])[which(result.stats[i,][paste(target.sig, "FDR q")] == min(result.stats[i,][paste(target.sig, "FDR q")]) )] )
				}
			}

	cat("State Assignment Complete \n")	
	#
	return(result.stats)
}



# Gene Set analysis
result.gene.sets <- run_analysis(
	counts.mat = counts.mat,
	type = "genes",
	n.states = 19,
	n.starts = 5,
	seed = 12345,
	max.q = 0.00001,
	min.coef = 1,
	min.spec = 0.1,
	min.sim = 0.4)



# Cell Cluster analysis
result.clusters <- run_analysis(
	counts.mat = counts.mat,
	type = "cells",
	n.states = 13,
	n.starts = 5,
	seed = 12345,
	max.q = 0.05,
	min.coef = 1,
	min.spec = 0.1,
	min.sim = 0)

