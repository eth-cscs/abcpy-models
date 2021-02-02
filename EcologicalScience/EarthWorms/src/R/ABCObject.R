################################################################################
##                        PRINTING ABC OBJECTS PRETTILY                       ## 
################################################################################

print.abcObject <- function(x, as.is = FALSE, num.items = 5,...) {
    if (num.items == 0) {
		num.items <- 1
    }
    for (i in 1:length(x)) {
		cat("$", names(x)[i], "\n", sep = "")
		if (!as.is & !is.call(x[[i]]) & !is.function(x[[i]]) &
            !("abcObject" %in% class(x[[i]])) &
            !("proc_time" %in% class(x[[i]])) & length(x[[i]]) > num.items) {
	    	cat("total: ", calcElements(x[[i]]), " elements\n", sep = "")
		}
		if (as.is | is.call(x[[i]]) | is.function(x[[i]]) |
            "abcObject" %in% class(x[[i]]) | "proc_time" %in% class(x[[i]])) {
	    	print(x[[i]])
		} else {
	    	if (is.array(x[[i]]) | is.data.frame(x[[i]])) {
				print(x[[i]][1:min(length(x[[i]][,1]), num.items),])
	    	} else {
				print(x[[i]][1:min(length(x[[i]]), num.items)])
	    	}
		}		
		if (!("list" %in% class(x[[i]])) & i != length(x)) {
	    	cat("\n")
		}
    }
}

calcElements <- function(x) {
    if (is.list(x)) {
		do.call(sum, lapply(x, calcElements))
    } else {
		length(x)
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

## This script accompanies the following paper:
## van der Vaart, E., Johnston, A.S.A., Beaumont, M.A., & Sibly, R.M.
## "Calibration and evaluation of individual-based models
## using Approximate Bayesian Computation" (2015) Ecological Modelling.