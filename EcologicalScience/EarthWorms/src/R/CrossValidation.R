################################################################################
##                CROSS-VALIDATION FOR ABC PARAMETER ESTIMATION               ## 
################################################################################

## create.abcEstCross() does cross-validation for ABC parameter estimation,
## with the following arguments:
## priors = the priors used to generate the simulation results;
## results = the simulation results, with rows corresponding to 'priors';
## rate = the acceptance rate, i.e., the proportion of runs to accept;
## model = the name of the model: the default is first element of 'priors';
## times = the number of cross-validations to perform;
## ss.p = the columns in 'priors' to use for ABC -
## (the default skips the first column, assuming that it is a model index);
## ss.r = the columns in 'results' to use for ABC

create.abcEstCross <- function(priors, results, rate, model = priors[1,1],
	times, ss.p = 2:length(priors[1,]), ss.r = 1:length(results[1,])) {
	print(paste("Doing cross-validation for parameter estimation. Total times: ",
		times, "!", sep = ""))

    ## store the original function call, so it can be viewed later
    call <- match.call()

    ## take the relevant subset of the provided results
    results.ss <- results[, ss.r]

	## randomly select the rows to use as "pseudo-data" for
	## the cross-validation
    c.indexes <- sample.int(length(priors[,1]), times)

	## store the runs used as "pseudo-data"
    c.targets <- results.ss[c.indexes,]
    
    ## store the remaining runs *not* used as pseudo-data    
    c.results <- results.ss[-c.indexes,]
    c.priors <- priors[-c.indexes,]
    
    ## calculate the number of runs to accept given the acceptance rate
    number.to.accept <- ceiling(dim(c.results)[1] * rate)
	
	## create a matrix to store ABC's estimates in
    estimates <- priors[c.indexes,]

    ## calculate the standard deviation of each column of simulation results
    c.sds <-  apply(c.results, 2, sd)
    scaled.c.results <- c.results
    scaled.c.results[, c.sds != 0] <- sweep(x = c.results[ , c.sds != 0],
        MARGIN = 2, STATS = c.sds[c.sds != 0], FUN = "/")

    ## do ABC parameter estimation for each set of "pseudo-data"	
    for (i in 1:times) {
    	print(paste("Working on time #", i, sep = ""))
		scaled.data <- c.targets[i,] / c.sds
        errors <- sweep(x = scaled.c.results, MARGIN = 2, STATS = scaled.data,
        	function(x, y) (y - x)^2)
        summed.errors <- sqrt(apply(X = errors, MARGIN = 1, FUN = sum))
		error.to.accept <- sort(summed.errors)[number.to.accept]
		accepted <- (summed.errors <= error.to.accept)
		
		## store the median of the accepted priors as ABC's posterior estimate
		estimates[i, ss.p] <- apply(X = c.priors[accepted, ss.p], MARGIN = 2,
            FUN = median)
    }
	
	## store the "pseudo-data" in its own object	
    true.values <- priors[c.indexes,]
    rownames(true.values) <- paste("row", c.indexes)
    
   	print(paste("Done doing cross-validation for parameter estimation!"))	
	
	## assemble an object describing what's been done, and return it to the user
	outcome <- c(list(call = call, rate = rate, selected = c.indexes,
		true.values = true.values, est.values = estimates, 
		ss.p = ss.p, ss.r = ss.r))

    class(outcome) <- c("abcObject", "abcEstCross")
    return(outcome)
}

## plot.abcEstCross() plots the true versus estimated parameter values
## for this paper's abcEstCross object 

## note that this is not a generally applicable function - it requires
## exactly 14 parameters and makes data-specific assumptions about digits

plot.abcEstCross <- function(x, sig = 0.01, ...) {
    cross.true <- x$true.values[ , x$ss.p]
    cross.est <- x$est.values[, x$ss.p]

    par.names <- c()
    
    for (i in 1:14) {
		if (grepl("_", colnames(cross.true)[i])) {
	    	parts <- strsplit(colnames(cross.true)[i], "_")
	    	fixed <- paste("bolditalic(", parts[[1]][1], "[", parts[[1]][2],
	    		"])", sep = "")
		} else {
	    	fixed <- paste("bolditalic(\"", colnames(cross.true)[i], "\"[])",
	    		sep = "") 
		}
		par.names <- c(par.names, parse(text = fixed))
    }
	
    mins <- c()
    maxes <- c()
	
    for (i in 1:14) {
		mins <- c(mins, min(cross.true[,i], cross.est[,i]))
		maxes <- c(maxes, max(cross.true[,i], cross.est[,i]))
    }
		
    min.lims <- mins - 0.05 * maxes
    max.lims <- 1.05 * maxes
		
    my.digits <- c(0, 3, 1, 1, 1, 2, 2, 3, 3, 2, 3, 3, 3, 3)	

    quartz(width = 8.8, height = 5.85)
    test <- layout(matrix(c(rep(c(1, rep(3, 5), rep(4, 5), rep(5, 5), rep(6, 5),
        rep(7, 5)), 5), rep(c(1, rep(8, 5), rep(9, 5), rep(10, 5), rep(11, 5),
        rep(12, 5)), 5), rep(c(1, rep(13, 5), rep(14, 5), rep(15, 5), rep(16, 5),
        rep(0, 5)), 5), c(0, rep(2, 25))), nrow = 16, ncol = 26, byrow = TRUE))
	
    par(mai = c(0, 0, 0, 0))

    plot(x <- c(0, 1), y <- c(0, 1), type = "n", ann = FALSE, axes = FALSE)
    text(0.5, 0.5, "estimated parameter value", srt = 90, cex = 1.8)
	
    plot(x <- c(0, 1), y <- c(0, 1), type = "n", ann = FALSE, axes = FALSE)
    text(0.5, 0.4, "true parameter value", cex = 1.8)

    par(mai = c(0.15, 0.2, 0.3, 0.1), lwd = 2, cex = 0.8)

    for (i in 1:14) {
		plot(cross.true[, i], cross.est[, i], ann = FALSE, axes = FALSE,
            ylim = c(min.lims[i], max.lims[i]),
            xlim = c(min.lims[i], max.lims[i]))
		axis(1, mgp = c(3, 0.4, 0), lwd.ticks = 2, at = c(mins[i], maxes[i]),
            labels = c(round(mins[i], my.digits[i]),
            round(maxes[i], my.digits[i])), cex.axis = 1.1, tcl = -0.25) 
		axis(2, mgp = c(3, 0.4, 0), lwd.ticks = 2, at = c(mins[i], maxes[i]),
            labels = c(round(mins[i], my.digits[i]),
            round(maxes[i], my.digits[i])), cex.axis = 1.1, tcl = -0.25) 
		title(main = par.names[i], line = 0.9, cex.main = 1.5)
		box()
		corr <- cor.test(cross.true[,i], cross.est[,i])
		p.value <- p.adjust(corr$p.value, n = 14, method = "holm")
		if (p.value > sig) {
	    	text(mins[i] + 0.22 * (max.lims[i] - min.lims[i]), max.lims[i],
                substitute(expression(paste(italic("r " ), "= ", V1)),
                list(V1 = format(corr$estimate, digits = 1, nsmall = 2)))[[2]],
                offset = 0.25, pos = 1, cex = 1.2)
		} else {
	    	text(mins[i] + 0.22 * (max.lims[i] - min.lims[i]), max.lims[i],
	    		substitute(expression(paste(italic("r " ), "= ", V1, "*")),
	    		list(V1 = format(corr$estimate, digits = 1, nsmall = 2)))[[2]],
	    		offset = 0.25, pos = 1, cex = 1.2)
		}
    }
}

################################################################################
##                   CROSS-VALIDATION FOR ABC MODEL SELECTION                 ## 
################################################################################

## create.abcSelCross() does cross-validation for ABC parameter model selection,
## indexes = the model indexes used to generate the simulation results;
## results = the simulation results, with rows corresponding to 'indexes';
## rate = the acceptance rate, i.e., the proportion of runs to accept;
## times = the number of cross-validations to perform;
## ss.r = the columns in 'results' to use for ABC

create.abcSelCross <- function(indexes, results, rate, times,
	ss.r = 1:length(results[1,])) {
	
	## store the original function call, so it can be viewed later
    call <- match.call()
    
    ## take the relevant subset of the provided results
    results.ss <- results[, ss.r]
    row.names(results.ss) <- c(1:length(results.ss[,1]))
	
	## find the names & numbers of unique models
	unique.models <- unique(indexes)
	number.of.models <- length(unique.models)
	
	print(paste("Doing cross-validation for model selection. Total times: ",
		times, ", for ", number.of.models, " models!", sep = ""))
	
	## select 'times' runs of each model as "pseudo-data"
	c.indexes <- c()
	for (i in 1:number.of.models) {
		these.indexes <- sample(as.numeric(rownames(
			results.ss[indexes == unique.models[i], ])), times)
		c.indexes <- c(c.indexes, these.indexes)
	}
	
	## store the runs used as "pseudo-data"
    c.targets <- results.ss[c.indexes,]
    
    ## store the remaining runs *not* used as pseudo-data    
    c.results <- results.ss[-c.indexes,]
    c.indexes <- indexes[-c.indexes]
    
    ## calculate the number of runs to accept given the acceptance rate
    number.to.accept <- ceiling(dim(c.results)[1] * rate)

	## create a matrix to store ABC's model classifications in
    classifications <- matrix(nrow = number.of.models, ncol = number.of.models)
    colnames(classifications) <- unique.models
    rownames(classifications) <- unique.models
    classifications[ , ] <- 0

    ## calculate the standard deviation of each column of simulation results
	c.sds <-  apply(c.results, 2, sd)
	c.sds[c.sds == 0] <- 1
    scaled.c.results <- sweep(x = c.results, MARGIN = 2, STATS = c.sds, FUN = "/")

    ## do ABC parameter estimation for each set of "pseudo-data"
	for (i in 1:number.of.models) {
		for (j in 1:times) {
			print(paste("Working on time #", j, " for model #", i, sep = ""))
			scaled.data <- c.targets[((i - 1) * times) + j, ] / c.sds
        	errors <- sweep(x = scaled.c.results, MARGIN = 2, STATS = scaled.data,
        	function(x, y) (y - x)^2)
        	summed.errors <- sqrt(apply(X = errors, MARGIN = 1, FUN = sum))
			error.to.accept <- sort(summed.errors)[number.to.accept]
			accepted <- (summed.errors <= error.to.accept)
			
			models.accepted <- c()
    		for (k in 1:number.of.models) {
				models.accepted <- c(models.accepted,
            	length(c.indexes[accepted][c.indexes[accepted] ==
            		unique.models[k]]))
    		}
    		
			## store the model accepted most often as ABC's classification
    		classifications[i, c(1:number.of.models)[models.accepted ==
    			max(models.accepted)]] <- classifications[i, c(1:number.of.models)
    			[models.accepted == max(models.accepted)]] + 1
    	}
    }
    
    print(paste("Done doing cross-validation for model selection!"))	
    
   	## assemble an object describing what's been done, and return it to the user
	outcome <- c(list(call = call, rate = rate, selected = c.indexes,
		classifications = classifications, ss.r = ss.r))
	class(outcome) <- c("abcObject", "abcSelCross")
	outcome
}

## plot.abcSelCross() plots the number of times each model was classified as
## each other model for this paper's abcSelCross object 

## note that this is not a generally applicable function - it assumes exactly
## two models and labels them as appropriate for this paper

plot.abcSelCross <- function(x, ..) {	
	quartz(width = 4.5, height = 3)
	layout(matrix(c(rep.int(1, 2), rep.int(2, 1)), 1, 3))
	par(mai = c(0.7, 0.7, 0.2, 0.1), lwd=2)
	mp <- barplot(t(x$classifications), axes = FALSE, ann = FALSE,
		axisnames = FALSE, col = c("white", "black" ))
	title(xlab = "model version", line = 3.25, cex.lab = 1.8)
	title(ylab = "number of classifications", line = 2.75, cex.lab = 1.8)
	axis(2, at = seq(0, sum(x$classifications[1,]), by = sum(x$classifications[1,]) / 4), lwd = 2, cex.axis = 1.8,
		mgp = c(3, 0.6, 0))
	axis(1, at = mp, lab =c ("Full", "Simple"), lwd = 2, cex.axis = 1.8,
		mgp = c(3, 1.2, 0), tick = FALSE)
	par(mai=c(0.2, 0.1, 0.2, 0.1), lwd = 2)
	plot(x <- c(0, 1), y <- c(0, 1), type = "n", ann = FALSE, axes = FALSE,
		ylim = c(0, 1), xlim = c(0, 1))
	legend(0.06, 0.9, legend = c("Full", "Simple"), col = c(""), bty = "n",
		fill = c("white", "black"), cex = 1.8)
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