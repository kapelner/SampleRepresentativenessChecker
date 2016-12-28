gender = c(0.5, 0.5)
ethnicity = c(.13, .70, .13, .04)
age = c(.13, .41, .30, .16)
region = c(.22, .18, .37, .23)

# install.packages("tensor")
library(tensor)
#vectorize the tensor multiplication
#pi = c(tensor(tensor(tensor(gender, ethnicity), age), region))
#illustratePower(pi = pi, Nsim = 500, verbose = TRUE)