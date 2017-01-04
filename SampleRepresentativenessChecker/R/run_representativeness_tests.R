#' Runs a single Equivalence Hypothesis Test with the null of non-representativeness for a fixed sample size.
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
#'		pi_gender = c(0.5, 0.5)
#'		pi_ethnicity = c(.13, .70, .13, .04)
#'		pi_age = c(.13, .41, .30, .16)
#'		pi_region = c(.22, .18, .37, .23)
#'		# install.packages("tensor")
#'		library(tensor)
#'		#vectorize via tensor multiplication
#'		pi = c(tensor(tensor(tensor(pi_gender, pi_ethnicity), pi_age), pi_region))
#'		#pretend to collect data
#' 		n = 1000
#'		x = rmultinom(1, n, pi)
#'		res = runSingleRepresentativenessTest(pi, x, delta = 0.03, alpha = 0.05)
#'		res$pval
#' 	}
#' @export
runSingleRepresentativenessTest = function(
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
#			 if (res[i] > -0.01){
#			   illustrateEquivalence(pi, delta, xsamp / sum(xsamp), ci_r, main = paste("pval = ", 2 * x[i]))
#			 }
		}
#		 print(cbind(x, res))
		res
	}
	pval = uniroot.all(pval_function, c(1e-6, 0.5 - 1e-6), xsamp = x, n = 100)
	if (length(pval) == 0){
		#now test really small alpha
		if (pval_function(1e-6, x) > 0){
			cat("pval < 1e-6 but no more accuracy available.\n")
			pval = 1e-6
		} else {
			pval = 0.5
		}
	}
	#uniroot.all may return multiple roots due to instabilities
	pval = min(pval)
	
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


#' Runs a set of L Equivalence Hypothesis Tests with the null of non-representativeness for a fixed sample size.
#' This makes use of the \code{runSingleRepresentativenessTest} function.
#' 
#' @param pi 			A list of length L where the entries are the true, gold-standard probability 
#' 						vectors which define the population characteristics.
#' 						We are testing against these to ascertain if the sample is representative.
#' @param x 			A list of length L where the entries are the experimental samples consisting of a count 
#' 						vector. The length of each vector must be equal to the corresponding vector length in \code{pi}.
#' @param delta 		The practical equivalence tolerance level. This can be a scalar (implying
#' 						a d-length vector with the same value) or a vector of length d.
#' @param alpha 		The minimum family-wise Type I error rate (test size). Default is \code{5\%}.
#' @param illustrate 	Should each of the L tests be illustrated? For each test, each of the d cells can be 
#' 						illustrated as a bar plot showing pi +- delta with the p-hat sample and its 
#' 						confidence interval. Default is code{TRUE}.
#' @return 				A list which contains function arguments as well as the following keys: 
#' 							(a) reject_null - a boolean which indicates
#' 							if the non-representativeness null was rejected
#' 							(b) ci - this is the simultaneous confidence interval for the
#' 							multinomial sample proportion vector as a d x 2 matrix 
#' 							(c) pval - the significance level for the family of tests
#' 							(d) pvals - the significance levels for each test using the 
#' 							Bonferroni-Holm correction
#' 
#' @author Adam Kapelner
#' @examples
#'  \dontrun{
#'		pi_gender = c(0.5, 0.5)
#'		pi_ethnicity = c(.13, .70, .13, .04)
#'		pi_age = c(.13, .41, .30, .16)
#'		pi_region = c(.22, .18, .37, .23)
#'		#pretend to collect data
#'		n = 1000
#'		x_gender = rmultinom(1, n, pi_gender)
#'		x_ethnicity = rmultinom(1, n, pi_ethnicity)
#'		x_age = rmultinom(1, n, pi_age)
#'		x_region = rmultinom(1, n, pi_region)
#'		res = runMultipleRepresentativenessTests(
#' 				list(pi_gender, pi_ethnicity, pi_age, pi_region),
#' 				list(x_gender, x_ethnicity, x_age, x_region),  
#' 				delta = 0.03,
#' 				alpha = 0.05)
#'		res$pval
#' 		#you can also see this by running all tests and then doing 
#' 		#the p-val correction at the end
#'		res1 = runSingleRepresentativenessTest(pi_gender, x_gender, 
#' 				delta = 0.03, alpha = 0.05)
#'		res2 = runSingleRepresentativenessTest(pi_ethnicity, x_ethnicity, 
#' 				delta = 0.03, alpha = 0.05)
#'		res3 = runSingleRepresentativenessTest(pi_age, x_age, 
#' 				delta = 0.03, alpha = 0.05)
#'		res4 = runSingleRepresentativenessTest(pi_region, x_region, 
#' 				delta = 0.03, alpha = 0.05)
#'		p.adjust(p = c(res1$pval, res2$pval, res3$pval, res4$pval), method = "holm")
#' 	}
#' @export
runMultipleRepresentativenessTests = function(
		pi,
		x,
		delta,
		alpha = 0.05,
		illustrate = TRUE
){
	if (class(pi) != "list"){
		stop("Argument \"pi\" must be a list of vectors with length p")
	}
	if (class(x) != "list"){
		stop("Argument \"x\" must be a list of vectors with length p")
	}
	p = length(pi)
	if (length(x) != p){
		stop("The arguments \"pi\" and \"x\" must be of the same length")
	}
	#run all tests individually and collect data from each run
	pvals = array(NA, p)
	cis = list()
	for (j in 1 : p){
		res = runSingleRepresentativenessTest(pi[[j]], x[[j]], delta, alpha, illustrate)
		pvals[j] = res$pval
		cis[[j]] = res$ci
	}
	#adjust the p-vals using the Bonferroni-Holm method i.e. the most powerful method
	#which makes no assumptions on the tests whatsoever
	pvals_adj = p.adjust(p = pvals, method = "holm")
	#the overall p-val for all rejctions can be thought of as the maximum value here
	pval = max(pvals_adj)
	list(
		reject_null = pval < alpha, 
		pval = pval, 
		pvals = pvals, 
		pvals_adj = pvals_adj, 
		cis = cis,
		pi = pi,
		x = x,
		delta = delta,
		alpha = alpha
	)
}