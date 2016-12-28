#' Runs the Equivalence Hypothesis Test with the null of non-representativeness for a fixed sample size.
#' This makes use of the package \code{MultinomialCI} which is based on the method of Sison & Glaz (1995).
#' 
#' @param pi 			The true, gold-standard probabilities which defines the population characteristics
#' 						on d cells and is thus of length d.
#' 						We are testing against this to ascertain if the sample is representative.
#' @param x 			The experimental sample consisting of a count vector of length d.
#' @param delta 		The practical equivalence tolerance level. This can be a scalar (implying
#' 						a d-length vector with the same value) or a vector of length d.
#' @param alpha 		The minimum Type I error rate (test size). Default is \code{5\%}.
#' @param illustrate 	Should the d cells be illustrated as a bar plot showing pi +- delta with
#' 						the p-hat sample and its confidence interval? Default is code{TRUE}.
#' @return 				A list which contains function arguments as well as the following keys: 
#' 							(a) reject_null - a boolean which indicates
#' 							if the non-representativeness null was rejected
#' 							(b) ci - this is the simultaneous confidence interval for the
#' 							multinomial sample proportion vector as a d x 2 matrix 
#' 							(c) pval - the significance level for the test
#' 							(d) conventional_null_chisq - the chi-squared value for the 
#' 							conventional goodness-of-fit test
#' 							(e) conventional_null_pval - the significance level based on
#' 							the conventional goodness-of-fit test
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
#'		res = runRepresentativenessTest(pi, x, delta = 0.03)
#'		res$pval
#' 	}
#' @export
runRepresentativenessTest = function(
	pi,
	x,
	delta,
	alpha = 0.05,
	illustrate = TRUE
){
	
	if (length(x) != length(pi)){
		stop("dimension of x != dimension of pi")
	}
	
	ci = multinomialCI(x, alpha = 2 * alpha, verbose = FALSE)
	ci_a = ci[, 1]
	ci_b = ci[, 2]
	pi_a = pi - delta
	pi_b = pi + delta
	
	if (all(ci_a > pi_a) && all(ci_b < pi_b)){
		reject_null = TRUE
	} else {
		reject_null = FALSE
	}
	
	##find p-val using Newton-Raphson
	pval_function = function(x, xsamp){
		#x comes in as a vector, so we have to vectorize everything in here
		res = array(NA, length(x))
		for (i in 1 : length(x)){
			ci_r = multinomialCI(xsamp, alpha = 2 * x[i], verbose = FALSE)
			ci_r_a = ci_r[, 1]
			ci_r_b = ci_r[, 2]
			res[i] = min(min(ci_r_a - pi_a), min(pi_b - ci_r_b))
			# if (res[i] > -0.01){
			#   illustrateEquivalence(pi, delta, xsamp / sum(xsamp), ci_r, main = paste("pval = ", 2 * x[i]))
			# }
		}
		# print(cbind(x, res))
		res
	}
	pval = uniroot.all(pval_function, c(0, 0.4999999999), xsamp = x)
	if (length(pval) == 0){
		pval = 0.5
	}
	
	##do GOF test as well
	n = sum(x)
	conventional_null_chisq = sum((x - n * pi)^2 / n * pi)
	conventional_null_pval = 1 - pchisq(conventional_null_chisq, length(x) - 1)
	
	if (illustrate){
		illustrateEquivalence(pi, delta, x, ci)
	}
	
	list(
		reject_null = reject_null,
		ci = cbind(ci_a, ci_b),
		pval = pval,
		conventional_null_chisq = conventional_null_chisq,
		conventional_null_pval = conventional_null_pval,
		pi = pi,
		x = x,
		delta = delta,
		alpha = alpha
	)
}
