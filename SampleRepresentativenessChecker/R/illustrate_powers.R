#' Runs simulations of sample size and delta to investigate the power of the equivalence test
#' 
#' @param pi 				The true, gold-standard probabilities which defines the population characteristics.
#' 							We are simulating under this scenario to see how often we can detect it (i.e., the test power).
#' @param alpha 			The minimum Type I error rate (test size). Default is \code{5\%}. 
#' @param ns 				Which sample sizes should be simulated? The default is \code{round(10^(seq(1.67, 3.33, by = 0.33)))}
#' 							corresponding to \code{47, 100, 214, 457, 977, 2089}.
#' @param deltas 			Which deltas (the practical equivalence range) should be simulated? The default
#' 							are the scalars \code{0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.10}.
#' @param Nsim 				How many simulations per n and delta combination? The higher the value, the
#' 							greater resolution. Default is \code{100} for speed.
#' @param plot 				Should we plot power by n and delta upon completion? Default is \code{TRUE}.
#' @param verbose			Should we print out messages indicating progress during a potentionally
#' 							long simulation? Default is \code{FALSE}.
#' @return 					Simulated powers in a matrix of dimension number of 
#' 							\code{ns} x number of \code{deltas}. 
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
#'		#vectorize the tensor multiplication
#'		pi = c(tensor(tensor(tensor(pi_gender, pi_ethnicity), pi_age), pi_region))
#'		illustratePowerSingleTest(pi = pi, Nsim = 500, verbose = TRUE)
#' 	}
#' @export
illustratePowerSingleTest = function(
	pi,
	alpha = 0.05,
	ns = round(10^(seq(1.67, 3.33, by = 0.33))),
	deltas = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.10),
	Nsim = 100,
	plot = TRUE,
	verbose = FALSE
){
	
	if (sum(pi) != 1){
		warning("sum of the pi vector was not 1. Normalizing...")
		pi = pi / sum(pi)
	}
	powers = matrix(NA, length(ns), length(deltas))
	rownames(powers) = ns
	colnames(powers) = deltas
	
	for (delta_i in 1 : length(deltas)){
		delta = deltas[delta_i]
		#get bounds now since they don't change
		pi_a = pi - delta
		pi_b = pi + delta
		
		if (verbose){cat("delta = ", delta, "\n")}
		for (n_i in 1 : length(ns)){
			n = ns[n_i]
			if (verbose){cat(" n = ", n, "\n")}
			
			res = array(0, Nsim)
			for (i in 1 : Nsim){
				x = rmultinom(1, size = n, prob = pi)
				#the TOST procedure gives us 2 * alpha here
				ci = multinomialCI(x, alpha = 2 * alpha, verbose = FALSE)
				
				ci_a = ci[, 1]
				ci_b = ci[, 2]
				
				if (verbose && plot && i %% 100 == 0){
					illustrateEquivalence(pi, delta, x, ci)
				}
				if (all(ci_a > pi_a) && all(ci_b < pi_b)){
					res[i] = 1
				}
			} 
			powers[n_i, delta_i] = mean(res)
		}
	}
	if (plot){
		plot(0, 0, 
				xlim = c(0, max(ns)), 
				ylim = c(0, 1), 
				xlab = "Sample Size",
				ylab = "Power",
				type = "n")
		for (delta_i in 1 : length(deltas)){
			lines(ns, powers[, delta_i])
		}
	}
	powers
}

#' Runs simulations of sample size and delta to investigate the power of multiple equivalence tests run simultaneously
#' 
#' @param pi 				A list of length L where the entries are the true, gold-standard probability 
#' 							vectors which define the population characteristics. We are simulating under 
#' 							these scenarios to see how often we can detect them while preserving the
#' 							family-wise error rate (i.e., the test power).
#' @param alpha 			The minimum family-wise Type I error rate (test size). Default is \code{5\%}. 
#' @param ns 				Which sample sizes should be simulated? The default is \code{round(10^(seq(1.67, 3.33, by = 0.33)))}
#' 							corresponding to \code{47, 100, 214, 457, 977, 2089}.
#' @param deltas 			Which deltas (the practical equivalence range) should be simulated? The default
#' 							are the scalars \code{0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.10}.
#' @param Nsim 				How many simulations per n and delta combination? The higher the value, the
#' 							greater resolution. Default is \code{100} for speed.
#' @param plot 				Should we plot power by n and delta upon completion? Default is \code{TRUE}.
#' @param verbose			Should we print out messages indicating progress during a potentionally
#' 							long simulation? Default is \code{FALSE}.
#' @return 					Simulated powers in a matrix of dimension number of 
#' 							\code{ns} x number of \code{deltas}. 
#' 
#' @author Adam Kapelner
#' @examples
#'  \dontrun{
#'		pi_gender = c(0.5, 0.5)
#'		pi_ethnicity = c(.13, .70, .13, .04)
#'		pi_age = c(.13, .41, .30, .16)
#'		pi_region = c(.22, .18, .37, .23)
#'		illustratePowerMultipleTests(
#' 			pi = list(pi_gender, pi_ethnicity, pi_age, pi_region), 
#' 			Nsim = 500, verbose = TRUE)
#' 	}
#' @export
illustratePowerMultipleTests = function(
		pi,
		alpha = 0.05,
		ns = round(10^(seq(1.67, 3.33, by = 0.33))),
		deltas = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.10),
		Nsim = 100,
		plot = TRUE,
		verbose = FALSE
){

	powers = matrix(NA, length(ns), length(deltas))
	rownames(powers) = ns
	colnames(powers) = deltas
	
	for (j in 1 : length(pi)){
		if (sum(pi[[j]]) != 1){
			warning("sum of the pi vector was not 1. Normalizing...")
			pi[[j]] = pi[[j]] / sum(pi[[j]])
		}	
	}
	
	for (delta_i in 1 : length(deltas)){
		delta = deltas[delta_i]
		
		if (verbose){cat("delta = ", delta, "\n")}
		for (n_i in 1 : length(ns)){
			n = ns[n_i]
			if (verbose){cat(" n = ", n, "\n")}
			
			res = array(0, Nsim)
			for (i in 1 : Nsim){
				x = list()
				for (j in 1 : length(pi)){
					x[[j]] = rmultinom(1, size = n, prob = pi[[j]])
				}
				
				if (verbose && plot && i %% 100 == 0){
					for (j in 1 : length(pi)){
						ci = multinomialCI(x[[j]], alpha = 2 * alpha, verbose = FALSE)
						illustrateEquivalence(pi[[j]], delta, x[[j]], ci)
					}
				}
				res[i] = ifelse(runMultipleRepresentativenessTests(pi, x, delta, alpha, illustrate = FALSE)$reject_null, 1, 0)
			} 
			powers[n_i, delta_i] = mean(res)
		}
	}
	if (plot){
		plot(0, 0, 
				xlim = c(0, max(ns)), 
				ylim = c(0, 1), 
				xlab = "Sample Size",
				ylab = "Power",
				type = "n")
		for (delta_i in 1 : length(deltas)){
			lines(ns, powers[, delta_i])
		}
	}
	powers
}
