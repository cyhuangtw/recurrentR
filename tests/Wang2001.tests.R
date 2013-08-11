# implemtent by the method based on "Acceptance Rejection"
gen_inhomo_poisson <- function(lambda, T_0) {
  result <- c()
  # give bound of lambda(t)
  # 0 < lambda(t) <= lambda_u
  lambda_u <- optimize(lambda, c(0, T_0), maximum = TRUE)$objective
  t <- 0
  while(t <= T_0) {
    u <- runif(1)   # u1: Uniform(0,1)
    t <- t - log(u) / lambda_u
    if (t > T_0) break
    u2 <- runif(1)  # u2: Uniform(0,1)
    # the chance p to keep generated t
    # if lambda(t) == lambda_u, p = 1 >= u2 (Uniform(0,1))
    # t will be definitely kept
    if (u2 <= lambda(t)/lambda_u) result <- append(result, t)
  }
  result
}

T_0 <- 30   # experiment end time
lambda.fun <- function(t) exp(sin(t) - 1)
lambda.fun <- function(t) exp(-t)
lambda.const <- function(t) 10

## find example for multiple events occurred in [t, t+1]
find_ex_seed <- function(i) {
  set.seed(i)
  gen_inhomo_poisson(lambda.fun, T_0)
}
# sapply(1:10, find_ex_seed)


## Validate F.hat ##########
library(plyr, quietly=TRUE)
tau <- T_0  # [0, tau] = time interval of interest
x.N <- 500  # number of patients
x.seq <- rlply(x.N, gen_inhomo_poisson(lambda=lambda.fun, T_0=T_0))
x.merged <- sort(unlist(x.seq))  # merge all events and sort

table(round(x.merged + 0.5, digits=0))

time.digit <- 1  # 0: integer, 1: with 1 decimals, ...
x.merged <- round(x.merged + 5 * 10^(-time.digit-1), digits=time.digit)
d <- table(x.merged)
s <- as.numeric(names(d))
N <- cumsum(d)

d[s > 6]
N[s > 6]
1 - d[s > 6] / N[s > 6]
prod(1 - d[s > 6] / N[s > 6])

F.hat <- function(t_seq) {
  sapply(t_seq, function(t) prod(1 - d[s > t] / N[s > t]))
}

Lambda.fun <- function(t) {
  sapply(t, function(t) integrate(lambda.fun, 0, t)$value)
}

F.fun <- function(t) {
  sapply(t, function(t) Lambda.fun(t) / Lambda.fun(T_0))
}

curve(F.fun, 0, tau, col = 2)
curve(F.hat, 0, tau, add=TRUE)

## Generate Censor Time and Validate Lambda.hat ##########
time.digit <- 1
x.seq <- sapply(x.seq, function(obs) round(obs + 5 * 10^(-time.digit-1), digits=time.digit))

# y.unobs ~ Uniform(0, 1.5 * T_0)
# for observation i, events berfore y.unobs[[i]] will be dropped, 
# if y > T_0 then no censor time y[[i]] == NA

# real censor time that cannot be observed in real world
y.unobs <- runif(x.N, 0.5 * T_0, 1.5 * T_0)   
censor_x <- function(x, y) x <- x[x < y]
x.seq.censored <- mapply(censor_x, x.seq, y.unobs)

y.obs <- ifelse(y.unobs > T_0, NA, y.unobs)  # observed y

x.merged <- sort(unlist(x.seq.censored))
d <- table(x.merged)
s <- as.numeric(names(d))
N <- cumsum(d)
curve(F.fun, 0, tau, col = 2)
curve(F.hat, 0, tau, add=TRUE)

y <- sapply(y.obs, function(y) min(y, tau, na.rm=TRUE))
m <- mapply(function(obs, y) length(obs[obs < y]), x.seq.censored, y)

Lambda.hat.T_0 <- sum(m / F.hat(y)) / x.N
# compared with theoretical value
Lambda.fun(T_0)

Lambda.hat.fun <- function(t_seq) {
  sapply(t_seq, function(t) F.hat(t) * Lambda.hat.T_0)
}

curve(Lambda.fun, 0, tau, col = 2)
curve(Lambda.hat.fun, 0, tau, add=TRUE)
