################################################################################
##                    COVERAGE FOR ABC PARAMETER ESTIMATION                   ## 
################################################################################

## create.abcEstCov() calculates coverage for ABC parameter estimation,
## with the following arguments:
## abcEstObj = the abcEst object to calculate coverage for;
## times = the number of coverage calculations to perform

create.abcEstCov <- function(abcEstObj, times) {
	print(paste("Calculating coverage. Total times: ", times, "!", sep = ""))
	
	## store the original function call, so it can be viewed later
	call <- match.call()
	
	## store the current time, so the final time to run can be calculated 
	start.time <- proc.time()
	
	## store the abcEstObj's priors and results in their own objects	
	all.priors <- abcEstObj$priors[, abcEstObj$ss.p]
	all.results <- abcEstObj$results[, abcEstObj$ss.r]
	rownames(all.priors) <- 1:length(all.priors[,1])
	rownames(all.results) <- 1:length(all.results[,1])

	## store the abcEstObj's *accepted* priors in its own object	
	acc.priors <- all.priors[abcEstObj$accepted, ]
	
	## create a matrix to store coverage's p-values in
	p.values <- matrix(nrow = times, ncol = length(all.priors[1,]))
	colnames(p.values) <- colnames(acc.priors)
	
	## calculate the number of runs to accept given the acceptance rate
    number.to.accept <- ceiling((dim(all.results)[1] - 1) * abcEstObj$rate)
	
	## create a matrix to store coverage's accepted runs in
	acc.indexes <- matrix(nrow = number.to.accept, ncol = times)

    ## do ABC parameter estimation for each set of "pseudo-data"
	for (i in 1:times) {
		print(paste("Working on time #", i, sep = ""))
		
		## select the accepted run to serve as "pseudo-data"
		target <- all.results[rownames(all.results) ==
			rownames(acc.priors)[i]]
		
		## select the runs that were *not* selected to serve as pseudo-data
		these.results <- all.results[rownames(all.results)
			!= rownames(acc.priors)[i], ]
		these.priors <-  all.priors[rownames(all.results) !=
			rownames(acc.priors)[i], ]
    	
    	## do ABC parameter estimation
    	sim.sds <-  apply(X = these.results, MARGIN = 2, FUN = sd)
	    scaled.results <- these.results
    	scaled.results[, sim.sds != 0] <- sweep(x = these.results[, sim.sds != 0],
        	MARGIN = 2, STATS = sim.sds[sim.sds != 0], FUN = "/")
    	scaled.data <- target
    	scaled.data[sim.sds != 0] <- target[sim.sds != 0] / sim.sds[sim.sds != 0]
    	errors <- sweep(x = scaled.results, MARGIN = 2, STATS = scaled.data,
        	function(x, y) (y - x)^2)
		summed.errors <- sqrt(apply(X = errors, MARGIN = 1, FUN = sum))	
    	error.to.accept <- sort(summed.errors)[number.to.accept]
    	accepted <- (summed.errors <= error.to.accept)
    	these.acc.priors <- these.priors[accepted, ]
    	
    	## for each parameter value, calculate the proportion of accepted values
    	## smaller than the true value
    	for (j in 1:length(acc.priors[1,])) {
    		true.value <- acc.priors[i, j]
    		p.value <- (1 + sum(these.acc.priors[, j] < true.value)) /
    			(2 + number.to.accept)
    		p.values[i, j] <- p.value
    	}
    	
    	## store the accepted runs for this set of "pseudo-data"
    	acc.indexes[, i] <- rownames(these.acc.priors)
	}
	
	## calculate the total time it took to calculate coverage
	total.time <- proc.time() - start.time
	print(paste("Done calculating coverage!"))	

	## assemble an object describing what's been done, and return it to the user
	outcome <- list(call = call, time = total.time,
		p.values = p.values, acc.indexes = acc.indexes)
	class(outcome) <- c("abcCovProp", "abcObject")
	outcome
}

## plot.abcEstCov() plots the p.values for all parameters as calculated by
## coverage - i.e., the proportion of accepted parameter values which was
## smaller than the true value.

## note that this is not quite a generally applicable function -
## it requires exactly 14 parameters

################################################################################
##                    COVERAGE FOR ABC PARAMETER ESTIMATION                   ## 
################################################################################

plot.abcEstCov <- function(x, sig = 0.01, ...) {
	par.names <- colnames(x$p.values)
	names.for.plot <- c()

	for (i in 1:14) {
		if (grepl("_", par.names[i])) {
			parts <- strsplit(par.names[i], "_")
			fixed <- paste("italic(", parts[[1]][1], "[", parts[[1]][2],
				"])", sep = "")
		} else {
			fixed <- paste("italic(\"", par.names[i], "\"[])", sep = "") 
		}
		names.for.plot <- c(names.for.plot, parse(text = fixed))
	}

    quartz(width = 8.8, height = 5.85)
    test <- layout(matrix(c(rep(c(1, rep(3, 5), rep(4, 5), rep(5, 5), rep(6, 5),
        rep(7, 5)), 5), rep(c(1, rep(8, 5), rep(9, 5), rep(10, 5), rep(11, 5),
        rep(12, 5)), 5), rep(c(1, rep(13, 5), rep(14, 5), rep(15, 5), rep(16, 5),
        rep(0, 5)), 5), c(0, rep(2, 25))), nrow = 16, ncol = 26, byrow = TRUE))
	
    par(mai = c(0, 0, 0, 0))

	plot(x.cors <- c(0, 1), y.cors <- c(0, 1), type = "n", ann = FALSE,
		axes = FALSE)
	text(0.5, 0.5, "number", srt = 90, cex = 1.8)
	
	plot(x.cors <- c(0, 1), y.cors <- c(0, 1), type = "n", ann = FALSE,
		axes = FALSE)
	text(0.5, 0.4, expression(paste(italic(p), " values")), cex = 1.8)

	par(mai = c(0.15, 0.2, 0.3, 0.1), lwd = 2, cex = 0.8)

	for (i in 1:14) {
		this.hist <- hist(x$p.values[, i], breaks = seq(0, 1, by = 0.05),
			ylim = c(0, 30), xlim = c(0, 1), axes = FALSE, ann = FALSE, col = "grey")
		axis(1, mgp = c(3, 0.4, 0), at = c(this.hist$mids[1], this.hist$mids[20]),
			lwd.ticks = 2, cex.axis = 1.1, tcl = -0.25)
		axis(2, mgp = c(3, 0.4, 0), at = c(0, 10, 20, 30), lwd.ticks = 2,
			cex.axis = 1.1, tcl = -0.25) 
		title(main = names.for.plot[i], line = 0.9, cex.main = 1.5)
		box(lwd = 2)
		p.value <- suppressWarnings(ks.test(x$p.values[,i], punif,
			alternative = "two.sided")$p.value)
		corrected <- p.adjust(p.value, method = "holm", n = 14)
		if (corrected < sig) {
			text(this.hist$mids[1], 29, "*", cex = 2)
		}
	}
}

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

## This script implements the algorithm initially described in:
## Prangle, D., Blum, M.G.B, Popovic, G. & Sisson, S.A.
## "Diagnostic tools for approximate Bayesian computation
## using the coverage property" (2010) 
## Australian & New Zealand Journal of Statistics, 56, 309 - 329.

## This script accompanies the following paper:
## van der Vaart, E., Johnston, A.S.A., Beaumont, M.A., & Sibly, R.M.
## "Calibration and evaluation of individual-based models
## using Approximate Bayesian Computation" (2015) Ecological Modelling.