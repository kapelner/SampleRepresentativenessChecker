#' Illustrates how equivalent the sample (p-hat) is assessed to be relative to the truth (pi)
#' based on the simultaneous multinomial proportion confidence interval.
#'  
#' @param pi 		The true, gold-standard probabilities which defines the population characteristics
#' 					on d cells and is thus of length d. The resulting illustration will plot these
#' 					as blue dots.
#' @param delta 	The practical equivalence tolerance level. This can be a scalar (implying
#' 					a d-length vector with the same value) or a vector of length d. The resulting 
#' 					illustration will plot this tolerance around pi as black bars.
#' @param x 		The experimental sample consisting of a count vector of length d. The default is
#' 					\code{NULL} as this argument is optional. If specified, the sample proportions will
#' 					be plotted in the resulting illustration as green dots.
#' @param ci 		The confidence interval for the experiment. This data can be grabbed from the
#' 					returned list in \code{runRepresentativenessTest} or it can be generated independently
#' 					by using the \code{multinomialCI} function. The default is \code{NULL} as this argument 
#' 					is optional. If specified, the interval will be plotted in the resulting illustration as 
#' 					red line segments. The representativeness test passes if all these intervals are completed
#' 					contained within the black intervals.
#' @param space 	Spacing between bars in the resulting barplot illustration. Default is \code{0.33}.
#' @param ... 		Additional arguments to be passed to the \code{barplot} function. 
#' 
#' @author Adam Kapelner
#' @examples
#'  \dontrun{
#'		gender = c(0.5, 0.5)
#'		ethnicity = c(.13, .70, .13, .04)
#'		age = c(.13, .41, .30, .16)
#'		region = c(.22, .18, .37, .23)
#'		# install.packages("tensor")
#'		library(tensor)
#'		#vectorize the tensor multiplication
#'		pi = c(tensor(tensor(tensor(gender, ethnicity), age), region))
#'		#pretend to collect data
#'		x = rmultinom(1, 1000, pi)
#'		delta = 0.03
#'		res = runRepresentativenessTest(pi, x, delta = delta)
#'		illustrateEquivalence(pi, delta, x, res$ci)
#' 	}
#' @export
illustrateEquivalence = function(pi, delta, x = NULL, ci = NULL, space = 0.33, ...){
	d = length(pi)
	pi_min = pmax(pi - delta, 0)
	pi_max = pmin(pi + delta, 1)
	
	if (is.null(ci)){
		y_max = max(pi_max)
	} else {
		y_max = max(pi_max, ci[, 2])
	}
	barplot(pi_max, 
			col = "black", 
			ylim = c(0, y_max), 
			xlab = "cells", ylab = "probabilities",
			space = space, border = NA,  ...)
	b = barplot(pi_min, col = "white", add = TRUE,  space = space, border = NA, ...)
	
	if (!is.null(ci)){
		for (i in 1 : d){
			segments(b[i], ci[i, 1], b[i], ci[i, 2], col = "red", lwd = 2)
		}
	}  
	if (!is.null(x)){
		points(x = b, y = x / sum(x), col = "green", pch = 16) #y = p-hat
	}
	points(x = b, y = pi, col = "blue", pch = 16)
}



