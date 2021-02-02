################################################################################
##                   ACTUALLY DOING ABC PARAMETER ESTIMATION                  ## 
################################################################################

## create.abcSel() does ABC model selection, with the following arguments:
## target = the data to fit;
## indexes = a list specifying which rows of results belong to which model;
## results = the simulation results;
## rate = the acceptance rate, i.e., the proportion of runs accepted
## ss.t = the columns in 'target' to use for ABC;
## ss.r = the columns in 'results' to use for ABC

## to use this function, store its output in an object,
## i.e, at the console, type:
## > my.abcSel <- create.abcSel(target = my.target, indexes = my.indexes,
##      results = my.results, rate = my.accept.rate)
## then look at the result by typing: > summary(my.abcEst)

create.abcSel <- function(target, indexes, results, rate,
	ss.t = 1:length(target), ss.r = 1:length(results[1,])) {
    ## store the original function call, so it can be viewed later
    call <- match.call()
    
    ## take the relevant subsets of all provided inputs
    target.ss <- unlist(target)[ss.t]
    results.ss <- results[ , ss.r]
	
    ## calculate the standard deviation of each column of simulation results
    sim.sds <- apply(X = results.ss, MARGIN = 2, FUN = sd)
    
    ## TODO
	## throw a warning if this occurs
    sim.sds[sim.sds == 0] <- 1
	
    ## scale the simulation results by their standard deviations
    ## (this to correct for different scales between different types of results)
    scaled.results <- sweep(x = results.ss, MARGIN = 2, STATS = sim.sds, FUN = "/")
	
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
    
    ## for each model, calculate the number of runs that were accepted
    ## unique(indexes) gives the number of different models
    models.accepted <- c()
    for (i in 1:length(unique(indexes))) {
		models.accepted <- c(models.accepted,
            length(indexes[accepted][indexes[accepted] == unique(indexes)[i]]))
    }
    
    ## assemble an object describing what's been done, and return it to the user
    outcome <- list(call = call, target = target, indexes = indexes,
        results = results, rate = rate, ss.t = ss.t, ss.r = ss.r,
        errors = summed.errors, accepted = accepted,
        models.accepted = models.accepted)
	class(outcome) <- c("abcObject", "abcSel")
    return(outcome)
}

################################################################################
##                         PRINT & SUMMARY FUNCTIONS                          ## 
################################################################################

## this function computes the main results
## (i.e., the Bayes factors) of an abcSel object;
## see this summary by typing > summary(name.of.object) in the console 

summary.abcSel <- function(x, ...) {
	accepted <- x$models.accepted
	factors <- matrix(ncol = length(accepted), nrow = length(accepted))
	colnames(factors) <- paste(unique(x$indexes))
	rownames(factors) <- paste(unique(x$indexes))
	for (i in 1:length(accepted)) {
		for (j in 1:length(accepted)) {
			factors[i, j] <- round(accepted[i] / accepted[j], 2)
		}
	}
	
	outcome <- list(call = x$call, num.sumstats = length(x$target),
	num.models = length(unique(x$indexes)),
	num.runs = length(x$results[,1]), rate = x$rate, factors = factors)
	
	class(outcome) <- "summary.abcSel"
	outcome
}

## this function prints the summary of an abcSel object,
## as produced by the function summary.abcSel()

print.summary.abcSel <- function(x, ...) {
	cat("call:\n")
	show(x$call)

	cat("\nacceptance rate:\t",  x$rate, sep = "")	
	cat("\n")
	cat("\n# of models:\t\t\t", x$num.models, sep = "")
	cat("\n# of sumstats:\t\t",  x$num.sumstats, sep = "")
	cat("\n# of runs:\t\t\t", x$num.runs, sep = "")
	
	cat("\n\n")
	cat("Bayes factors:\n")
	print(x$factors)
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