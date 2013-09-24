PoisMixClusWrapper <- function(y, u, gmin = 1, gmax, conds, lib.size = TRUE, lib.type = "TMM",
	gmin.init.type = "small-em", init.runs = 5, split.init = TRUE, alg.type = "EM", cutoff = 10e-6, iter = 1000,
	fixed.lambda = NA, equal.proportions = FALSE, verbose = FALSE)
{
	all.results <- vector("list", length = gmax - gmin + 1)
	names(all.results) <- paste("g=", seq(gmin,gmax, 1), sep = "")

	## For gmin, run PoisMixClus with regular small-EM initialization
	cat("Running g =", gmin, "...\n")
	all.results[[1]] <- PoisMixClus(y = y, u = u, g = gmin, lib.size = lib.size,
		lib.type = lib.type, conds = conds,
		init.type = gmin.init.type,
		alg.type = alg.type, cutoff = cutoff, iter = iter,
		fixed.lambda = fixed.lambda, equal.proportions = equal.proportions,
		prev.labels = NA, prev.probaPost = NA, init.runs = 5, verbose = verbose)

	## For g > gmin, run PoisMixClus with Panos init using previous results
	index <- 2
	if(gmax > gmin) {
		if(split.init == TRUE) {
			for(K in seq((gmin+1),gmax,1)) {
				cat("Running g =", K, "...\n")
				prev.labels <- all.results[[K-1]]$labels
				prev.probaPost <- all.results[[K-1]]$probaPost
				all.results[[index]] <- PoisMixClus(y = y,
                                        u = u, g = K, lib.size = lib.size,
					lib.type = lib.type, conds = conds,
					init.type = "split.small-em",
					alg.type = alg.type, cutoff = cutoff,
                                        iter = iter,
                                        fixed.lambda = fixed.lambda,
					equal.proportions = equal.proportions,
					prev.labels = prev.labels,
                                        prev.probaPost = prev.probaPost,
					init.runs = 5, verbose = verbose)
				index <- index + 1
			}
		}
		if(split.init == FALSE) {
			for(K in seq((gmin+1),gmax,1)) {
				cat("Running g =", K, "...\n")
				all.results[[index]] <- PoisMixClus(y = y, u = u,g = K, lib.size = lib.size,
					lib.type = lib.type, conds = conds,
					init.type = gmin.init.type,
					alg.type = alg.type, cutoff = cutoff, iter = iter,
					fixed.lambda = fixed.lambda,
					equal.proportions = equal.proportions,
					prev.labels = NA, prev.probaPost = NA, init.runs = 5,
					verbose = verbose)
				index <- index + 1
			}
		}
	}

	logLike.all <- unlist(lapply(all.results, function(x) x$log.like))
        p.logLike.all <- unlist(lapply(all.results, function(x) x$p.log.like))
        entropy.all <- unlist(lapply(all.results, function(x) x$entropy))
        entropyM.all <- unlist(lapply(all.results, function(x) x$entropyM))
        entropyS.all <- unlist(lapply(all.results, function(x) x$entropyS))
        BIC.all <- unlist(lapply(all.results, function(x) x$BIC))
	BIC.choose <- which(BIC.all == max(BIC.all, na.rm = TRUE))
	BIC.select.results <- all.results[[BIC.choose]]
        ICL.all <- unlist(lapply(all.results, function(x) x$ICL))
	ICL.choose <- which(ICL.all == max(ICL.all, na.rm = TRUE))
	ICL.select.results <- all.results[[ICL.choose]]
        SICL.all <- unlist(lapply(all.results, function(x) x$SICL))
	SICL.choose <- which(SICL.all == max(SICL.all, na.rm = TRUE))
	SICL.select.results <- all.results[[SICL.choose]]
	MIL.all <- unlist(lapply(all.results, function(x) x$MIL))
	MIL.choose <- which(MIL.all == max(MIL.all, na.rm = TRUE))
	MIL.select.results <- all.results[[MIL.choose]]
	RESULTS <- list(logLike.all = logLike.all,
                        p.logLike.all = p.logLike.all,
                        entropy.all = entropy.all, entropyM.all = entropyM.all,
                        entropyS.all = entropyS.all, BIC.all = BIC.all,
                 ICL.all = ICL.all, SICL.all = SICL.all, MIL.all = MIL.all,
		BIC.select.results = BIC.select.results,  ICL.select.results = ICL.select.results, SICL.select.results = SICL.select.results, MIL.select.results = MIL.select.results, all.results = all.results)
	class(RESULTS) <- "HTSClusterWrapper"
	return(RESULTS)
}
