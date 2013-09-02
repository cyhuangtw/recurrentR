##################
# Initialization #
##################

library(recurrentR)
library(RPPGen)
library(microbenchmark)

set.seed(1000)

# Lambda function
lambda <- function(x) exp(-x/10)
T_0 <- rpois(1, 40)
gen_zi <- function() runif(1, 0.5, 1.5)
n <- 150
beta <- c(1, -1)
X <- cbind(sin(1:n), sample(c(0, 1), n, TRUE))

y <- rpois(n, T_0)
y <- as.numeric(ifelse(y < T_0, y, T_0))
t <- sapply(seq_along(y), function(i) {
	# browser()
	# z <- gen_z() * exp(as.vector(X[i,] %*% beta))
	zi <- gen_zi() * exp((X[i,] %*% beta)[[1]])
	lambda_i <- function(t) zi * lambda(t)
	retval <- gen_inhomo_poisson(lambda_i, y[i] - 1, lambda_i(0))
	if (is.null(retval)) return(vector("numeric", 0))
	return(retval)
	# return(as.numeric(ceiling(retval)))
})

# round up for some demicals, so two event could be same in simulation
# for DEBUG only
time.digit <- 2  # 0: integer, 1: with 1 decimals, ...
t <- sapply(t, function(t) {
	round(t + 5 * 10^(-time.digit-1), digits=time.digit)
})

# not use W yet, passed with an empty data frame
obj <- new("recurrent-data", X, y, t, data.frame(), T_0)

################
# Q.hat, R.hat #
################

# Q.hat() and Q.hat.c() do not behave the same!
q_hat_gen <- recurrentR:::Q.hat(obj)
q_hat_gen.c <- recurrentR:::Q.hat.c(obj)
time.steps <- q_hat_gen.c$x[q_hat_gen.c$x >= 0.95 & q_hat_gen.c$x <= 1]
q_hat_gen.c$y[q_hat_gen.c$x %in% time.steps]
q_hat_gen(time.steps)

microbenchmark(
	q_hat_gen <- recurrentR:::Q.hat(obj),
	q_hat_gen.c <- recurrentR:::Q.hat.c(obj),
	times=5L
)
# Unit: milliseconds
#                        expr       min         lq     median         uq        max neval
#     recurrentR:::Q.hat(obj) 537.58193 546.787917 565.205981 594.809339 625.133629     5
#   recurrentR:::Q.hat.c(obj)   3.61891   4.379248   4.733636   5.215207   5.570719     5

r_hat_gen <- recurrent:::R.hat(obj)
r_hat_gen.c <- recurrent:::R.hat.c(obj)
r_hat_gen.c$[r_hat_gen.c$x %in% time.steps]
r_hat_gen(time.steps)

microbenchmark(
	recurrentR:::R.hat(obj),
	recurrentR:::R.hat.c(obj),
	times=5L
)
# Unit: milliseconds
#                       expr        min         lq     median         uq        max neval
#    recurrentR:::R.hat(obj) 755.695248 776.170171 776.399472 821.658846 852.265457     5
#  recurrentR:::R.hat.c(obj)   3.302091   3.689455   3.982129   4.056021   4.807083     5
