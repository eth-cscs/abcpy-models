## This script contains the full set of R commands necessary to generate
## all of the results presented the following paper:

## van der Vaart, E., Johnston, A.S.A., Beaumont, M.A., & Sibly, R.M.
## "Calibration and evaluation of individual-based models
## using Approximate Bayesian Computation" (2015) Ecological Modelling

################################################################################
##                    GETTING READY FOR PARAMETER ESTIMATION                  ## 
################################################################################

## set an absolute path to the home directory of this set of scripts; the
## recommendation is to avoid shortcuts like '~' as this can introduce problems
f.path <- "/Users/Elske/Gitbucket/FigWorms"

## load all simulation results for the full model - 1e6 rows by 15 and 160
## columns; not advised on machines with less than 16 GB of RAM
if (TRUE) {
	priors.full <- readRDS(paste(f.path, "/results/priors_full.rds", sep = ""))
	results.full.1 <- readRDS(paste(f.path, "/results/results_full_1.rds",
		sep = ""))
	results.full.2 <- readRDS(paste(f.path, "/results/results_full_2.rds",
		sep = ""))
	results.full <- rbind(results.full.1, results.full.2)
}

## on machines with less than 16 GB of RAM, load only the first 100,000 rows
## of the full model; supplied as text files for readability
if (FALSE) {
	priors.full <- read.table(paste(f.path, "/results/priors_full_1e5.txt",
		sep = ""), header = TRUE)
	results.full <- read.table(paste(f.path, "/results/results_full_1e5.txt",
		sep = ""), header = TRUE)
}

## load all empirical data; see 'PosteriorChecking.R' for references
all.data <- read.table(paste(f.path, "/results/all_data.txt", sep = ""))

## specify what experiments the columns of the simulation results and
## empirical data refer to; see end of file for references
cols.g02 <- 1:27
cols.ge03 <- 28:88
cols.rv90v.m <- seq(90, 122, 2)
cols.rv90v.c <- c(seq(89, 123, 2), 124)
cols.rv90c.m <- seq(126, 158, 2)
cols.rv90c.c <- c(seq(125, 159, 2), 160)
all.cols <- list(cols.g02, cols.rv90v.m, cols.rv90c.m, cols.ge03, cols.rv90v.c,
	cols.rv90c.c)

## load the packages & R scripts necessary to do ABC parameter estimation
library(car)
library(Hmisc)
source(paste(f.path, "/src/R/ABCObject.R", sep = ""))
source(paste(f.path, "/src/R/ParameterEstimation.R", sep = ""))
source(paste(f.path, "/src/R/CrossValidation.R", sep = ""))
source(paste(f.path, "/src/R/Coverage.R", sep = ""))

################################################################################
##                          DOING PARAMETER ESTIMATION                        ## 
################################################################################

## run 'Getting Ready for Parameter Estimation' before running this code section

## create an abcEst object for the 'full model' to do ABC parameter estimation
full.abcEst <- create.abcEst(target = all.data, priors = priors.full,
	results = results.full, rate = 0.0001, model = "Full")
	
## look at the full.abcEst object
print(full.abcEst)

## show a summary of the full.abcEst object (i.e., its posteriors)
summary(full.abcEst)

## plot Figure 3
plot(full.abcEst)

################################################################################
##                       CHECKING PARAMETER ESTIMATION                        ## 
################################################################################

## run 'Doing Parameter Estimation' before running this code section

## create an abcEstCross object to do cross-validation for ABC parameter
## estimation; for time reasons, this example only uses 5 runs as "pseudo-data"
full.abcEstCross.mini <- create.abcEstCross(priors = priors.full,
	results = results.full, rate = 0.0001, times = 5)

## plot the mini cross-validation for ABC parameter estimation
plot(full.abcEstCross.mini)

## read in and plot the ABC parameter estimation cross-validation in Figure 5
full.abcEstCross <- readRDS(paste(f.path, "/results/cross_est.rds", sep = ""))
plot(full.abcEstCross)

## calculate and print the correlations in ABC's posteriors,
## the basis of Figure 6
corr.full.abcEst <- create.postCorrs(full.abcEst)
print(corr.full.abcEst)

## create an abcEstCov object to calculate coverage for ABC parameter
## estimation; for time reasons, this example only uses 5 runs as "pseudo-data"
full.abcEstCov.mini <- create.abcEstCov(abcEstObj = full.abcEst, times = 5)

## plot the mini coverage calculation
plot.abcEstCov(full.abcEstCov.mini)

## read in and plot the coverage calculation in Figure 7
full.abcEstCov <- readRDS(paste(f.path, "/results/cov_est.rds", sep = ""))
plot.abcEstCov(full.abcEstCov)

################################################################################
##                      GETTING READY FOR MODEL SELECTION                     ## 
################################################################################

## run 'Getting Ready for Parameter Estimation' before running this code section

## loading all simulation results for the full model - 1e6 rows by 12 and 160
## columns; not advised on machines with less than 16 GB of RAM
if (TRUE) {
	priors.simple <- readRDS(paste(f.path, "/results/priors_simple.rds",
		sep = ""))
	results.simple.1 <- readRDS(paste(f.path, "/results/results_simple_1.rds",
		sep = ""))
	results.simple.2 <- readRDS(paste(f.path, "/results/results_simple_2.rds",
		sep = ""))
	results.simple <- rbind(results.simple.1, results.simple.2)
}

## on machines with less than 16 GB of RAM, load only the first 100,000 rows
## of the simple model; supplied as text files for readability
if (FALSE) {
	priors.simple <- read.table(paste(f.path, "/results/priors_simple_1e5.txt",
		sep = ""), header = TRUE)
	results.simple <- read.table(paste(f.path,
		"/results/results_simple_1e5.txt", sep = ""), header = TRUE)
}

## create a vector of all model indexes
all.indexes <- c(rep("Full", 1e6), rep("Simple", 1e6))

## combine all model results into a single matrix
all.results <- rbind(results.full, results.simple)

## load the R script for ABC model selection
source(paste(f.path, "/src/R/ModelSelection.R", sep = ""))

################################################################################
##                            DOING MODEL SELECTION                           ## 
################################################################################

## run 'Getting Ready for Model Selection' before running this code section

## create an abcSel object, which does ABC model selection
all.abcSel <- create.abcSel(target = all.data, indexes = all.indexes,
	results = all.results, rate = 0.0001)
	
## look at the all.abcSel object
print(all.abcSel)

## show a summary of the all.abcSel object (i.e., it's Bayes factors; Table 2)
summary(all.abcSel)

################################################################################
##                          CHECKING MODEL SELECTION                          ## 
################################################################################

## run 'Doing Model Selection' before running this code section

## create an abcSelCross object to do cross-validation for ABC model
## selection; for time reasons, this example only uses 3 runs as "pseudo-data"
all.abcSelCross.mini <- create.abcSelCross(indexes = all.indexes,
	results = all.results, rate = 0.0001, times = 3)

## plot the mini cross-validation for ABC model selection
plot(all.abcSelCross.mini)

## read in and plot the ABC model selection cross-validation in Figure 8
all.abcSelCross <- readRDS(paste(f.path, "/results/cross_sel.rds", sep = ""))
plot(all.abcSelCross)

################################################################################
##                          DOING POSTERIOR CHECKING                          ## 
################################################################################

## run 'Doing Parameter Estimation' before running this code section

## load in the necessary RNetLogo package
library(RNetLogo)
	
## set the path to the NetLogo application
nl.path <- "/Applications/NetLogo"

## load the R script for doing posterior checking with the earthworm IBM
source(paste(f.path, "/src/R/PosteriorChecking.R", sep = ""))

## do a posterior check, recreating Figure 4
## note that this function draws randomly from ABC's accepted runs, and then
## re-runs the earthworm IBM; as there's randomness in the IBM, the results
## may differ slightly from the paper; it may also take a while to run
plot.postCheck(abcEstObj = full.abcEst, draws = 100, rerun = TRUE)

################################################################################
##                         THE SUPPLEMENTARY MATERIAL                         ## 
################################################################################

## run 'Doing Parameter Estimation', 'Getting Ready for Model Selection' and
## 'Doing Posterior Checking' before running this code section

## to re-create Figure S1, the posterior distributions of the 'simple model',
## create an abEst object and plot it
simple.abcEst <- create.abcEst(target = all.data, priors = priors.simple,
	results = results.simple, rate = 0.0001, model = "Simple")
plot(simple.abcEst)

## to re-create Figure S2, the posterior check for the 'simple model', open
## 'PosteriorChecking.R' and load the simple models by running the code in the
## 'if (FALSE)' block, then use:
plot.postCheck(abcEstObj = simple.abcEst, draws = 100, rerun = TRUE)

## to re-create Figure S3, the posterior distributions of the 'full model'
## with a wider tolerance, create an abcEst object and plot it
full.abcEst.001 <- create.abcEst(target = all.data, priors = priors.full,
	results = results.full, rate = 0.001, model = "Full")
plot(full.abcEst.001)

## to re-create Table S2, ABC model selection with a wider tolerance,
## create an abcSel object and summarize it
all.abcSel.001 <- create.abcSel(target = all.data, indexes = all.indexes,
	results = all.results, rate = 0.001)
summary(all.abcSel.001)

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