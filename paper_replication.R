# Section 4
library(SampleRepresentativenessChecker)

#declare the different pi vectors (Fig 2)
pi_gender = c(0.5, 0.5)
pi_ethnicity = c(.13, .70, .13, .04)
pi_age = c(.13, .41, .30, .16)
pi_region = c(.22, .18, .37, .23)

# Section 4.1
#install.packages("tensor")
library(tensor)
#assume independent to find joint probabilities
#vectorize the tensor product
pi = c(tensor(tensor(tensor(pi_gender, pi_ethnicity), pi_age), pi_region))
#Fig 3
illustratePowerSingleTest(pi = pi, Nsim = 500, verbose = TRUE)
#run test on a random sample at the agreed-upon sample size
n = 1500
x = rmultinom(1, n, pi)
res = runSingleRepresentativenessTest(pi, x, delta = 0.03, alpha = 0.05)
res$reject_null
res$pval
#Fig 4
illustrateEquivalence(pi, delta = 0.03, x = x, ci = res$ci)

# Section 4.2

#Fig 5
illustratePowerMultipleTests(
		pi = list(pi_gender, pi_ethnicity, pi_age, pi_region), 
		Nsim = 500, verbose = TRUE)

#run test on a random sample at the agreed-upon sample size
n = 2000
x_gender = rmultinom(1, n, pi_gender)
x_ethnicity = rmultinom(1, n, pi_ethnicity)
x_age = rmultinom(1, n, pi_age)
x_region = rmultinom(1, n, pi_region)
res = runMultipleRepresentativenessTests(
		list(pi_gender, pi_ethnicity, pi_age, pi_region),
		list(x_gender, x_ethnicity, x_age, x_region),  
		delta = 0.03,
		alpha = 0.05)
res$reject_null
res$pval

#Fig 6
par(mfrow = c(2, 2))
illustrateEquivalence(pi_gender, delta = 0.05, x = x_gender, ci = res$ci)
illustrateEquivalence(pi_ethnicity, delta = 0.05, x = x_ethnicity, ci = res$ci)
illustrateEquivalence(pi_age, delta = 0.05, x = x_age, ci = res$ci)
illustrateEquivalence(pi_region, delta = 0.05, x = x_region, ci = res$ci)
