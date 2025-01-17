% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_representativeness_tests.R
\name{runSingleRepresentativenessTest}
\alias{runSingleRepresentativenessTest}
\title{Runs a single Equivalence Hypothesis Test with the null of non-representativeness for a fixed sample size.
This makes use of the package \code{MultinomialCI} which is based on the method of Sison & Glaz (1995).}
\usage{
runSingleRepresentativenessTest(pi, x, delta, alpha = 0.05,
  illustrate = TRUE)
}
\arguments{
\item{pi}{The true, gold-standard probabilities which defines the population characteristics
on d cells and is thus of length d.
We are testing against this to ascertain if the sample is representative.}

\item{x}{The experimental sample consisting of a count vector of length d.}

\item{delta}{The practical equivalence tolerance level. This can be a scalar (implying
a d-length vector with the same value) or a vector of length d.}

\item{alpha}{The minimum Type I error rate (test size). Default is \code{5\%}.}

\item{illustrate}{Should the d cells be illustrated as a bar plot showing pi +- delta with
the p-hat sample and its confidence interval? Default is code{TRUE}.}
}
\value{
A list which contains function arguments as well as the following keys: 
							(a) reject_null - a boolean which indicates
							if the non-representativeness null was rejected
							(b) ci - this is the simultaneous confidence interval for the
							multinomial sample proportion vector as a d x 2 matrix 
							(c) pval - the significance level for the test
							(d) conventional_null_chisq - the chi-squared value for the 
							conventional goodness-of-fit test
							(e) conventional_null_pval - the significance level based on
							the conventional goodness-of-fit test
}
\description{
Runs a single Equivalence Hypothesis Test with the null of non-representativeness for a fixed sample size.
This makes use of the package \code{MultinomialCI} which is based on the method of Sison & Glaz (1995).
}
\examples{
 \dontrun{
	pi_gender = c(0.5, 0.5)
	pi_ethnicity = c(.13, .70, .13, .04)
	pi_age = c(.13, .41, .30, .16)
	pi_region = c(.22, .18, .37, .23)
	# install.packages("tensor")
	library(tensor)
	#vectorize via tensor multiplication
	pi = c(tensor(tensor(tensor(pi_gender, pi_ethnicity), pi_age), pi_region))
	#pretend to collect data
	n = 1000
	x = rmultinom(1, n, pi)
	res = runSingleRepresentativenessTest(pi, x, delta = 0.03, alpha = 0.05)
	res$pval
	}
}
\author{
Adam Kapelner
}

