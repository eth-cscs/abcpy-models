################################################################################
##                          DOING POSTERIOR CHECKING                          ## 
################################################################################

## This script contains many of the R commands necessary to perform the
## posterior checks described in the following paper:

## van der Vaart, E., Johnston, A.S.A., Beaumont, M.A., & Sibly, R.M.
## "Calibration and evaluation of individual-based models
## using Approximate Bayesian Computation" (2015) Ecological Modelling

## To actually plot any results, use the script "Commands.R". See end of file
## for references concerning the NetLogo models and empirical data used. 

################################################################################
##                    GETTING READY FOR POSTERIOR CHECKING                    ## 
################################################################################

## load in the full models for the four experiments described in the paper;
## f.path is a path variable that should be set in 'Commands.R'

#### Gunadi et al. 2002 - Full
g2002.1.netlogo <- "g2002.1.netlogo"
NLStart(nl.path, gui = FALSE, nl.obj = g2002.1.netlogo)

g2002.1.model <- paste(f.path, "/src/models/g2002_1.nlogo", sep = "")
NLLoadModel(g2002.1.model, nl.obj = g2002.1.netlogo)

#### Gunadi & Edwards 2003 - Full
ge2003.1.netlogo <- "ge2003.1.netlogo"
NLStart(nl.path, gui = FALSE, nl.obj = ge2003.1.netlogo)

ge2003.1.model <- paste(f.path, "/src/models/ge2003_1.nlogo", sep = "")
NLLoadModel(ge2003.1.model, nl.obj = ge2003.1.netlogo)

#### Reinecke & Viljoen 1990, Variable Condition - Full
rv1990v.1.netlogo <- "rv1990v.1.netlogo"
NLStart(nl.path, gui = FALSE, nl.obj = rv1990v.1.netlogo)

rv1990v.1.model <- paste(f.path, "/src/models/rv1990v_1.nlogo", sep = "")
NLLoadModel(rv1990v.1.model, nl.obj = rv1990v.1.netlogo)

#### Reinecke & Viljoen 1990, Constant Condition - Full
rv1990c.1.netlogo <- "rv1990c.1.netlogo"
NLStart(nl.path, gui = FALSE, nl.obj = rv1990c.1.netlogo)

rv1990c.1.model <- paste(f.path, "/src/models/rv1990c_1.nlogo", sep = "")
NLLoadModel(rv1990c.1.model, nl.obj = rv1990c.1.netlogo)

## loading in the simple models is only necessary for reproducing the
## paper's Supplementary Information - just set 'FALSE' to 'TRUE'

if (FALSE) {
	#### Gunadi et al. 2002 - Simple
	g2002.0.netlogo <- "g2002.0.netlogo"
	NLStart(nl.path, gui = FALSE, nl.obj = g2002.0.netlogo)

	g2002.0.model <- paste(f.path, "/src/models/g2002_0.nlogo", sep = "")
	NLLoadModel(g2002.0.model, nl.obj = g2002.0.netlogo)

	#### Gunadi & Edwards 2003 - Simple
	ge2003.0.netlogo <- "ge2003.0.netlogo"
	NLStart(nl.path, gui = FALSE, nl.obj = ge2003.0.netlogo)

	ge2003.0.model <- paste(f.path, "/src/models/ge2003_0.nlogo", sep = "")
	NLLoadModel(ge2003.0.model, nl.obj = ge2003.0.netlogo)

	#### Reinecke & Viljoen 1990, Variable Condition - Simple
	rv1990v.0.netlogo <- "rv1990v.0.netlogo"
	NLStart(nl.path, gui = FALSE, nl.obj = rv1990v.0.netlogo)

	rv1990v.0.model <- paste(f.path, "/src/models/rv1990v_0.nlogo", sep = "")
	NLLoadModel(rv1990v.0.model, nl.obj = rv1990v.0.netlogo)

	#### Reinecke & Viljoen 1990, Constant Condition - Simple
	rv1990c.0.netlogo <- "rv1990c.0.netlogo"
	NLStart(nl.path, gui = FALSE, nl.obj = rv1990c.0.netlogo)

	rv1990c.0.model <- paste(f.path, "/src/models/rv1990c_0.nlogo", sep = "")
	NLLoadModel(rv1990c.0.model, nl.obj = rv1990c.0.netlogo)
}

################################################################################
##                        PLOTTING THE POSTERIOR CHECK                        ## 
################################################################################

plot.postCheck <- function(abcEstObj, draws = 5, rerun = FALSE) {	
	indexes <- sort(sample.int(length(abcEstObj$priors[abcEstObj$accepted, 1]),
		draws, replace = TRUE))
		
	all.data <- as.vector(unlist(abcEstObj$target))
	
	results.lit <- data.frame(matrix(nrow = draws,
		ncol = length(abcEstObj$results[1,])))
	rsqs.lit <- data.frame(matrix(nrow = draws, ncol = 6))
	
	results.best <- data.frame(matrix(nrow = draws,
		ncol = length(abcEstObj$results[1,])))
	rsqs.best <- data.frame(matrix(nrow = draws, ncol = 6))
	
	results.ppc <- data.frame(matrix(nrow = draws,
		ncol = length(abcEstObj$results[1,])))
	rsqs.ppc <- data.frame(matrix(nrow = draws, ncol = 6))
			
	for (i in 1:draws) {
		results.lit[i,] <- do.run(make.default(abcEstObj$priors[1,1]))
		
		if (rerun) {
			results.best[i,] <- do.run(as.data.frame(t(abcEstObj$best)))
			results.ppc[i,] <-
				do.run(as.data.frame(t(abcEstObj$priors[abcEstObj$accepted,]
				[indexes[i],])))
		} else {
			results.best[i,] <- abcEstObj$results[abcEstObj$errors ==
				min(abcEstObj$errors),]
			results.ppc[i,] <-
				abcEstObj$results[abcEstObj$accepted,][indexes[i],]
		}
		for (j in 1:6) {
			ss.tot <- sum((all.data[all.cols[[j]]] -
				mean(all.data[all.cols[[j]]]))^2)
			rsqs.lit[i, j] <- 1 - (sum((all.data[all.cols[[j]]] -
				results.lit[i, all.cols[[j]]])^2) / ss.tot)
			rsqs.best[i, j] <- 1 - (sum((all.data[all.cols[[j]]] -
				results.best[i, all.cols[[j]]])^2) / ss.tot)
			rsqs.ppc[i, j] <- 1 - (sum((all.data[all.cols[[j]]] -
				results.ppc[i, all.cols[[j]]])^2) / ss.tot)
		}
	}

	my.colors <- c("black", "dimgrey")
	my.x <- list(c(0:26)*7, c(0:16)*10, c(0:16)*10, c(0:60)*7, c(0:18)*10,
		c(0:18)*10)
	my.xmins <- list(-1, -1, -1, -2.5, -1,-1)
	my.xmaxes <- list(183, 183, 183, 432.5, 183, 183)
	my.ymaxes <- list(0.75, 0.75, 0.75, 0.75, 60, 60)
	my.ylabs <- list("mass (grams)", "mass (grams)", "mass (grams)",
		"mass (grams)", "number of cocoons", "number of cocoons")
	my.xats <- list(seq(0, 180, by = 60), seq(0, 180, by = 60),
		seq(0, 180, by = 60), seq(0, 420, by = 140), seq(0, 180, by = 60),
		seq(0, 180, by = 60))
	my.yats <- list(c(0, 0.25, 0.5, 0.75), c(0, 0.25, 0.5, 0.75),
		c(0, 0.25, 0.5, 0.75), c(0, 0.25, 0.5, 0.75), c(0, 20, 40, 60),
		c(0, 20, 40, 60))
	my.legends <- list("A", "B", "C", "D", "E", "F")
	my.rsq.locs <- list("bottomright", "topleft", "bottomright",
		"bottomright", "topleft", "bottomright")
	my.arrow.x <- list(c(0), c(0, 10, 60, 140), seq(0, 160, 20),
		c(0, 161, 315), c(0, 10, 60, 140), seq(0, 160, 20))
	my.arrow.y <- list(c(-0.1977273, -0.1568182),
		c(-0.1977273, -0.1568182), c(-0.1977273, -0.1568182),
			c(-0.1977273, -0.1568182), c(-15.81818, -12.54545),
			c(-15.81818, -12.54545))

	quartz(width = 9, height = 6)
	par(mai = c(0.7, 0.6, 0.15, 0.1), lwd=2, mfrow = c(2,3))
	
	for (i in 1:6) {
		plot(x <- my.x[[i]], all.data[all.cols[[i]]], type = "n", pch = 20,
			xlim = c(my.xmins[[i]], my.xmaxes[[i]]), axes = FALSE, ann = FALSE,
			ylim = c(0, my.ymaxes[[i]]))
		title(xlab = "time (days)", line = 4, cex.lab = 1.8)
		title(ylab = my.ylabs[[i]], line = 2.75, cex.lab = 1.8)
		axis(2, at = my.yats[[i]], lwd = 2, cex.axis = 1.8, mgp = c(3, 0.8, 0))
		axis(1, at = my.xats[[i]], lwd = 2, cex.axis = 1.8, mgp = c(3, 1.3, 0))
		box(lwd = 2)
		legend("topright", legend = my.legends[[i]], text.font = 2, cex = 1.6,
			bty = "n", inset = c(-0.015, -0.0375))
		legend(my.rsq.locs[[i]], legend = make.rsqs.texts(c(colMeans(rsqs.lit[i]),
			colMeans(rsqs.best[i]))), lty = 1, lwd = 4, col = my.colors,
			bty = "n", cex = 1.5, seg.len = 0.4, x.intersp = 0.6,
			y.intersp = 1.2, inset = c(0.02, 0))
		for (j in 1:draws) {
			lines(my.x[[i]], results.ppc[j, all.cols[[i]]], lwd = 4, 
				col = rgb(0.75, 0.75, 0.75, 0.2), type = "l", lty = 1)
		}
		for (j in 1:length(my.arrow.x[[i]])) {
			arrows(my.arrow.x[[i]][j], my.arrow.y[[i]][1], my.arrow.x[[i]][j],
				my.arrow.y[[i]][2], xpd = TRUE, length = 0.05)
		}
		if (i == 2 || i == 5) {
			arrows(100, my.arrow.y[[i]][2], 100, my.arrow.y[[i]][1], xpd = TRUE,
				length = 0.05)
		}
		lines(my.x[[i]], colMeans(results.lit[all.cols[[i]]]), cex = 1.5,
			type = "l", col = "black", lwd = 4)
		lines(my.x[[i]], colMeans(results.best[all.cols[[i]]]), cex = 1.5,
			type = "l", col = "dimgrey", lwd = 4)
		lines(my.x[[i]], all.data[all.cols[[i]]], cex = 1.2, type = "o",
			col = "black", bg = "white", pch = 21)
	}
}

make.rsqs.texts <- function(values) {
	texts = vector('expression', length(values))
	texts[1] <- substitute(expression(paste(MYVALUE1, " ", bar(italic(R)^2))),
		list(MYVALUE1 = round(values[1], 2)))[2]
	texts[2] <- substitute(expression(paste(MYVALUE1, " ", bar(italic(R)^2))),
		list(MYVALUE1 = round(values[2], 2)))[2]
	texts
}

################################################################################
##                         RUNNING THE NETLOGO MODELS                         ## 
################################################################################

do.run <- function(this_row) {
	result <- c(do.run.g2002(this_row), do.run.ge2003(this_row),
		do.run.rv1990v(this_row), c(do.run.rv1990c(this_row)))
	result
}

do.run.g2002 <- function(this_row) {
	instance <- make.instance(this_row[1,1], "g2002")
	set.parameters(this_row, instance)
	NLCommand("setup-interface", nl.obj = instance)
	NLCommand("setup", nl.obj = instance)
	birth <- NLReport("mean-mass", nl.obj = instance)
	mass <- unlist(NLDoReport(iterations = 26, command = "go-7",
		reporter = "mean-mass", nl.obj = instance))
	result <- c(birth, mass)
	result
}

do.run.ge2003 <- function(this_row) {
	instance <- make.instance(this_row[1,1], "ge2003")
	set.parameters(this_row, instance)
	NLCommand("setup-interface", nl.obj = instance)
	NLCommand("setup", nl.obj = instance)
	birth <- NLReport("mean-mass", nl.obj = instance)
	mass <- unlist(NLDoReport(iterations = 60, command = "go-7",
		reporter = "mean-mass", nl.obj = instance))
	result <- c(birth, mass)
	result
}

do.run.rv1990v <- function(this_row) {
	instance <- make.instance(this_row[1,1], "rv1990v")
	set.parameters(this_row, instance)
	NLCommand("setup-interface", nl.obj = instance)
	NLCommand("setup", nl.obj = instance)
	NLCommand("go-25", nl.obj = instance)
	start <- c(NLReport("sum-hatchlings", nl.obj = instance),
		NLReport("mean-mass", nl.obj = instance))
	both <- unlist(NLDoReport(iterations=16, command="go-10",
		reporter = c("sum-hatchlings", "mean-mass"), nl.obj = instance))
	cocoons <- unlist(NLDoReport(iterations=2, command="go-10",
		reporter = c("sum-hatchlings"), nl.obj = instance))
	result <- c(start, both, cocoons)
	result
}

do.run.rv1990c <- function(this_row) {
	instance <- make.instance(this_row[1,1], "rv1990c")
	set.parameters(this_row, instance)
	NLCommand("setup-interface", nl.obj = instance)
	NLCommand("setup", nl.obj = instance)
	NLCommand("go-25", nl.obj = instance)
	start <- c(NLReport("sum-hatchlings", nl.obj = instance),
		NLReport("mean-mass", nl.obj = instance))
	both <- unlist(NLDoReport(iterations=16, command="go-10",
		reporter = c("sum-hatchlings", "mean-mass"), nl.obj = instance))
	cocoons <- unlist(NLDoReport(iterations=2, command="go-10",
		reporter = c("sum-hatchlings"), nl.obj = instance))
	result <- c(start, both, cocoons)
	result
}

################################################################################
##                         RUNNING THE NETLOGO MODELS                         ## 
################################################################################

set.parameters <- function(values, instance) {
	NLCommand(paste("set B_0 ",  round(values[2], 6), sep = ""),
		nl.obj = instance)
	NLCommand(paste("set activation_energy ", round(values[3], 6), sep = ""),
		nl.obj = instance)
	NLCommand(paste("set energy_tissue ", round(values[4], 6), sep = ""),
		nl.obj = instance)
	subtraction <- 1
	if (values[1] == 1) {
		NLCommand(paste("set energy_food ", round(values[5], 6), sep = ""),
			nl.obj = instance)
		NLCommand(paste("set half_saturation_coef ", round(values[7], 6),
			sep = ""), nl.obj = instance)
		NLCommand(paste("set speed ", round(values[15], 3), sep = ""),
			nl.obj = instance)
		subtraction <- 0
	}
	NLCommand(paste("set energy_synthesis ", round(values[6 - subtraction], 6),
		sep = ""), nl.obj = instance)
	if (values[1] == 0) {
		subtraction <- 2
	}
	NLCommand(paste("set max_ingestion_rate ", round(values[8 - subtraction], 
		6), åsep = ""), nl.obj = instance)
	NLCommand(paste("set mass_birth ", round(values[9 - subtraction], 6),
		sep = ""), nl.obj = instance)
	NLCommand(paste("set mass_cocoon ", round(values[10 - subtraction], 6),
		sep = ""), nl.obj = instance)

	NLCommand(paste("set mass_maximum ", round(values[11 - subtraction], 6),
		sep = ""), nl.obj = instance)	
	NLCommand(paste("set mass_sexual_maturity ",
		round(values[12 - subtraction], 6), sep = ""), nl.obj = instance)
	NLCommand(paste("set growth_constant ", round(values[13 - subtraction], 6),
		sep = ""), nl.obj = instance)
	NLCommand(paste("set max_reproduction_rate ",
		round(values[14 - subtraction], 6), sep = ""), nl.obj = instance)
}

make.instance <- function(model, experiment) {
	instance <- NULL
	
	if (model == 0) {
		if (experiment == "g2002") { instance = g2002.0.netlogo }
		if (experiment == "ge2003") { instance = ge2003.0.netlogo }
		if (experiment == "rv1990v") { instance = rv1990v.0.netlogo }
		if (experiment == "rv1990c") { instance = rv1990c.0.netlogo }
	} else {
		if (experiment == "g2002") { instance = g2002.1.netlogo }
		if (experiment == "ge2003") { instance = ge2003.1.netlogo }
		if (experiment == "rv1990v") { instance = rv1990v.1.netlogo }
		if (experiment == "rv1990c") { instance = rv1990c.1.netlogo }
	}
	instance
}

make.default <- function(model) {
	values <- data.frame(
		M = model,
		B_0 = 967,
		E = 0.25,
		E_c = 7,
		E_f = 10.6,
		E_s = 3.6,
		h = 3.5,
		IGm = 0.7,
		M_b = 0.011,
		M_c = 0.015,
		M_m = 0.5,
		M_p = 0.25,
		r_B = 0.177,
		r_m = 0.182,
		s = 0.0004
	)
	if (model == 0) {
		values$IGm[1] <- 0.15
		values <- values[, !colnames(values) %in% c("E_f", "h", "s")]
	}
	values
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

## This script accompanies the following paper:
## van der Vaart, E., Johnston, A.S.A., Beaumont, M.A., & Sibly, R.M.
## "Calibration and evaluation of individual-based models
## using Approximate Bayesian Computation" (2015) Ecological Modelling.

## ...and the earthworm IBM it runs was originally described here:
## Johnston, A.S.A., Hodson, M.E., Thorbek, P., Alvarez, T. & Sibly, R.M.
## "An energy budget agent-based model of earthworm populations and its
## application to study the effects of pesticides"
## (2015) Ecological Modelling, 280, 5 - 17.

## ...and finally, the empirical data it plots was taken from:
## Gunadi, B., Blount, C. & Edwards, C.A.
## "The growth and fecundity of Eisenia fetida (Savigny) in cattle solids
## pre-composted for different periods" (2002) Pedobiologia, 46, 15 - 23.

## Gunadi, B. & Edwards, C.A.
## "The effects of multiple applications of different organic wastes on the
## growth, fecundity and survival of Eisenia fetida (Savigny) (Lumbricidae)"
## (2003) Pediobiologia, 47, 321 - 329.

## Reinecke, A.J. & Viljoen, S.A.
## "The influence of feeding patterns on growth and reproduction of the
## vermicomposting earthworm Eisenia fetida (Oligochaeta)" (1990)
## Biology and Fertility of Soils, 10, 184 - 187.