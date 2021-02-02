################################################################################
##                   ACTUALLY DOING ABC PARAMETER ESTIMATION                  ## 
################################################################################

## create.abcEst() does ABC parameter estimation, with the following arguments:
## target = the data to fit;
## priors = the priors used to generate the simulation results;
## results = the simulation results, with rows corresponding to 'priors';
## rate = the acceptance rate, i.e., the proportion of runs to accept;
## model = the name of the model: the default is first element of 'priors';
## ss.t = the columns in 'target' to use for ABC;
## ss.p = the columns in 'priors' to use for ABC -
## (the default skips the first column, assuming that it is a model index);
## ss.r = the columns in 'results' to use for ABC

## to use this function, store its output in an object,
## i.e, at the console, type:
## > my.abcEst <- create.abcEst(target = my.target, priors = my.priors,
##      results = my.results, rate = my.accept.rate)
## then look at the result by typing: > summary(my.abcEst)

create.abcEst <- function(target, priors, results, rate, model = priors[1,1],
	ss.t = 1:length(target), ss.p = 2:length(priors[1,]),
	ss.r = 1:length(results[1,])) {
    
    ## store the original function call, so it can be viewed later
    call <- match.call()
	
    ## take the relevant subsets of the provided inputs
    target.ss <- unlist(target)[ss.t]
    results.ss <- results[, ss.r]
	
    ## calculate the standard deviation of each column of simulation results
    sim.sds <- apply(X = results.ss, MARGIN = 2, FUN = sd)
    
    ## don't scale columns with standard deviations of 0
    if (length(sim.sds[sim.sds == 0]) > 0) {
		warning("at least one column of simulation results has a
			standard deviation of 0; not scaling")
		    sim.sds[sim.sds == 0] <- 1
	}
	
    ## scale the simulation results by their standard deviations
    ## (this to correct for different scales between different types of results)
    scaled.results <- sweep(x = results.ss, MARGIN = 2,
    	STATS = sim.sds, FUN = "/")
	
    ## scale the empirical data in the same way as the simulation results
    scaled.data <- target.ss / sim.sds
    		
    ## for each row, calculate the distance (or error) from each simulation
    ## result to the corresponding empirical data point
    errors <- sweep(x = scaled.results, MARGIN = 2, STATS = scaled.data,
        function(x, y) (y - x)^2)
	
    ## for each row, sum the distances between each simulation result and the
    ## corresponding empirical data point (giving the total error)
    summed.errors <- sqrt(apply(X = errors, MARGIN = 1, FUN = sum))

    ## calculate the number of runs to accept given the acceptance rate
    number.to.accept <- ceiling(dim(errors)[1] * rate)
	
    ## calculate the maximum error to accept given the number of runs to accept
    error.to.accept <- sort(summed.errors)[number.to.accept]
	
    ## accept the runs with a total error less than the maximum acceptable error
    accepted <- (summed.errors <= error.to.accept)
    
    ## find the run with the smallest total error
    best.values <- priors[summed.errors == min(summed.errors),]
	
    ## calculate the median values of the accepted parameter values
    med.values <- best.values
    medians <- apply(X = priors[accepted, ss.p], MARGIN = 2, FUN = median)
    med.values[ss.p] <- medians
	
    ## assemble an object describing what's been done, and return it to the user
    outcome <- list(call = call, target = target, priors = priors,
        results = results, rate = rate, model = model, ss.t = ss.t, ss.p = ss.p,
        ss.r = ss.r, errors = summed.errors, accepted = accepted,
        best.values = best.values, med.values = med.values)
    class(outcome) <- c("abcObject", "abcEst")
    return(outcome)
}

################################################################################
##                       PLOT, PRINT & SUMMARY FUNCTIONS                      ## 
################################################################################

## this function computes the main results
## (i.e., the posterior parameter distributions) of an abEst object;
## see this summary by typing > summary(name.of.object) in the console 

summary.abcEst <- function(x, cred.int = 0.95, ...) {
    acc.params <- x$priors[x$accepted, x$ss.p]
    mins <- apply(acc.params, 2, min)
    lows <- apply(acc.params, 2, quantile, ((1 - cred.int) / 2))
    medians <- apply(acc.params, 2, quantile, 0.5)
    means <- apply(acc.params, 2, mean)
    highs <- apply(acc.params, 2, quantile, (1 - ((1 - cred.int) / 2)))
    maxes <- apply(acc.params, 2, max)
    sums <- rbind(mins, lows, medians, means, highs, maxes)
    rownames(sums) <- c("Min.:", paste((((1 - cred.int) / 2) * 100),
        "% Perc.:", sep = ""), "Median:", "Mean:",
        paste(((1 - (1 - cred.int) / 2) * 100), "% Perc.:", sep = ""), "Max.:")
	
    outcome <- list(call = x$call, rate = x$rate, model = x$model,
        sumstats = length(x$target), params = length(acc.params),
        runs = length(x$results[1,]), sums = sums)
    class(outcome) <- "summary.abcEst"
    outcome
}

## this function prints the summary of an abcEst object,
## as produced by the function summary.abcEst()

print.summary.abcEst <- function(x, digits = 3, ...) {
    cat("call:\n")
    show(x$call)

    cat("\ntype of model:\t\t\t'", x$model, "'", sep = "")
    cat("\nacceptance rate:\t\t",  x$rate, sep = "")	
    cat("\n")
    cat("\n# of sumstats:\t\t\t", x$sumstats, sep = "")
    cat("\n# of parameters:\t\t", x$params, sep = "")
    cat("\n# of runs:\t\t\t\t", x$runs, "\n\n", sep = "")
	
    cat("posteriors:\n")
    print(signif(x$sums, digits))
}

## this function plots the main results
## (i.e., the posterior parameter distributions) of an abcEst object
## see this plot by typing > plot(name.of.object) in the console 

plot.abcEst <- function(x, par.names = colnames(x$priors)[x$ss.p],
    cred.int = 0.95, sig = 0.01, ...) {
    params <- x$priors[ , x$ss.p]
    num.params <- length(params[1,])

    prior.medians <- apply(params, 2, median)
    posteriors <- params[x$accepted, ]
	
    scaled.priors <- sweep(x = params, MARGIN = 2, STATS = prior.medians, "/")
    scaled.posteriors <- sweep(x = posteriors, MARGIN = 2,
        STATS = prior.medians, "/")

    names.for.plot <- c()

    for (i in 1:num.params) {
		if (grepl("_", par.names[i])) {
			parts <- strsplit(par.names[i], "_")
			fixed <- paste("italic(", parts[[1]][1], "[", parts[[1]][2], "])",
				sep = "")
			if (length(parts[[1]]) > 2) {
				warning("too many subscripts in par.names; only the first used")
			}
		} else {
			fixed <- paste("italic(\"", par.names[i], "\"[])", sep = "") 
		}
		names.for.plot <- c(names.for.plot, parse(text = fixed))
    }
    sigs <- c()
	
    for (i in 1:num.params) {
		merged <- c(params[,i], posteriors[,i])
		groups <- as.factor(c(rep("PRIOR", dim(params)[1]),
            rep("POST", dim(posteriors)[1])))
		p.value <- leveneTest(merged, groups)$"Pr(>F)"[1]
		corrected <- p.adjust(p.value, method = "holm", n = num.params)
	
		if (corrected < sig) {
	    	sigs <- c(sigs, TRUE)
		} else {
	    	sigs <- c(sigs, FALSE)
		}
    }

    low.posts <- apply(scaled.posteriors, 2, quantile, ((1 - cred.int) / 2))
    top.posts <- apply(scaled.posteriors, 2, quantile, (1 - (1 - cred.int) / 2))

    low.priors <- apply(scaled.priors, 2, quantile, ((1 - cred.int) / 2))
    top.priors <- apply(scaled.priors, 2, quantile, (1 - (1 - cred.int) / 2))

    scaled.prior.medians <- apply(scaled.priors, 2, median)
    scaled.posterior.medians <- apply(scaled.posteriors, 2, median)

    x.cors <- c(0.75, seq(1.75, num.params, by = 1))
    dev.new(width = 0.95 + ((num.params + 1) * 0.35), height = 2)
    par(mai = c(0.6, 0.6, 0.15, 0.15), lwd = 2)
    
    plot(x.cors, scaled.prior.medians, type = "n", axes = FALSE, col = "grey90", 
        xlim = c(0.125 * (0.75 + num.params * 0.35),
        (num.params + 0.75 - (0.125 * (0.75 + num.params * 0.35)))),
        ylim = c(min(low.priors, low.posts) - 0.1 * max(top.priors, top.posts),
        max(top.priors, top.posts) + 0.1 * max(top.priors, top.posts)))
	
	segments(x.cors, low.priors, x.cors, top.priors, col = "grey60",
            lty = 2)
	segments(x.cors - 0.25, low.priors, x.cors + 0.25, low.priors,
            col = "grey60", lty = 2)
	segments(x.cors - 0.25, top.priors, x.cors + 0.25, top.priors,
            col = "grey60", lty = 2)
	points(x.cors, scaled.prior.medians, pch = 21, bg = "grey60",
            col = "grey60")

	segments(x.cors, low.posts, x.cors, top.posts)
	segments(x.cors - 0.25, low.posts, x.cors + 0.25, low.posts)
	segments(x.cors - 0.25, top.posts, x.cors + 0.25, top.posts)
	points(x.cors, scaled.posterior.medians, pch = 21, bg = "white")
		
	box()
	for (i in 1:num.params) {
	    axis(1, at = x.cors[i], labels = names.for.plot[i], cex.axis = 0.9,
                lwd.ticks = 2, mgp = c(3, 0.6, 0), las = 1, tck = -0.05)
	    if (sigs[i]) {
			text(x.cors[i], 1.05 * max(top.priors, top.posts), "*", cex = 2)
	    }
	}
	title(xlab = "parameter", line = 1.5, cex.lab = 0.9)
	title(ylab = "scaled value", line = 1.75, cex.lab = 0.9)
	axis(2, cex.axis = 0.9, lwd.ticks = 2, mgp = c(3, 0.5, 0), las = 1,
            tck = -0.05)
}

################################################################################
##                            CORRELATION FUNCTION                            ## 
################################################################################

create.postCorrs <- function(abcEstObj) {
	call <- match.call()
	num.params <- length(abcEstObj$priors[1, abcEstObj$ss.p])
	correlations <- rcorr(as.matrix(abcEstObj$priors[abcEstObj$accepted,
		abcEstObj$ss.p]), type = "spearman")
	corrected <- apply(correlations$P, 2, p.adjust,
		n = (num.params * num.params) / 2 - num.params, method = "holm")
	outcome <- list(call = call, corrs = correlations$r, p.values = corrected)
	class(outcome) <- "postCorrs"
	outcome 
}

print.postCorrs <- function(x, ...) {
    cat("call:\n")
    show(x$call)
    
    cat("\nSpearman correlations:\n")
    print(round(x$corrs, 2))
    cat("\nCorrected p-values:\n")
    print(round(x$p.values, 2))
}

## This script uses code from the R package 'abc', version 2.0
## http://cran.r-project.org/web/packages/abc/index.html
## Original Copyright (C) 2010 - 2015
## Katalin Csillery, Lemaire Louisiane, Francois Oliver, Michael Blum
## <michael.blum at imag.fr>

## Copyright (C) 2015 Elske van der Vaart, Mark Beaumont, Alice Johnston,
## Richard Sibly <elskevdv at gmail.com>

## This script is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This script is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## The R package 'abc' was originally described in:
## Csillery, K., Francois, O. & Blum, M.G.B.
## "abc: An R package for approximate Bayesian computation (ABC)"
## (2010) Methods in Ecology and Evolution, 3, 475 - 479.

## This script accompanies the following paper:
## van der Vaart, E., Johnston, A.S.A., Beaumont, M.A., & Sibly, R.M.
## "Calibration and evaluation of individual-based models
## using Approximate Bayesian Computation" (2015) Ecological Modelling.