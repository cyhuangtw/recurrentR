# implementation of generating function from the slides:
# http://rpubs.com/wush978/7008
# theorem reference: 
# http://filebox.vt.edu/users/pasupath/papers/nonhompoisson_streams.pdf
# the following code is adapted from Wush Wu

library(recurrentR)

# === Verify the Inhomogeneous Poisson Random Process Generator ===
T_0 <- 20   # experiment end time
T_seq.num <- 10   # total time steps for testing 
# generate a sequence of time steps of length T_seq.num
# [0, ..., t_i, ... , T_0] with 0 excluded
T_seq <- seq(T_0/T_seq.num, T_0, length.out=T_seq.num)

# define the rate function
rate.fun <- function(t) exp(sin(t) - 1)
# rate.fun <- function(t) exp(-t)


require(plyr)
sample.N <- 1000
# use bulit-in function
# sample.seq <- replicate(sample.N, gen_inhomo_poisson(lambda=rate.fun, T_0=T_0))
# use wrapper by plyr for progress report
sample.seq <- rlply(sample.N, 
                    gen_inhomo_poisson(lambda=rate.fun, T_0=T_0),
                    .progress='text')


count_until_t <- function(sample_seq, t) {
  sapply(sample_seq, function(x) sum(x <= t))
}

count_along_t <- function(sample, t_seq){
  sapply(t_seq, function(t) sum(sample <= t))
}

df.count <- ldply(sample.seq, function(s) count_along_t(s, T_seq), .progress="text")
colnames(df.count) <- paste("T=", T_seq, sep="")


# === quantile-quantile plot visual verification ===

# generate the sequence of cumulative rate 
# CRF(t_i) = $\int_0^t_i\lambda(t)dt$ for t_i in T_seq
# from property of inhomogeneous Poisson process,
# N[0,t_i] ~ Poisson(cum_rate.seq[[i]])
cum_rate.seq <- sapply(T_seq, function(end_t) integrate(rate.fun, 0, end_t)$value)

# use system bulti-in qqplot() function
qt.N <- 501  # number of quantiles
i <- 5 # i-th time step
plot_qq <- function(qt.N, i){
  qqplot(qpois(seq(0, 1, length.out=qt.N), lambda=cum_rate.seq[[i]]),
         df.count[[i]],
         xlab='Therotical', ylab='Generated')
  qqline(df.count[[i]], distribution=function(p) qpois(p, lambda=cum_rate.seq[[i]]))
}
plot_qq(qt.N, i)

# use ggplot2 for scatter plot
require(ggplot2)
plot_qq.ggplot <- function(qt.N, i) {
  g <- qplot(qpois(seq(0, 1, length.out=qt.N),  cum_rate.seq[[i]]), 
             quantile(df.count[[i]], seq(0, 1, length.out=qt.N)), 
             geom='point') 
  g + expand_limits(x=0, y=0) + labs(x="Therotical", y="Generated")
}
plot_qq.ggplot(qt.N, i)+theme_bw()


# TODO
# Kolmogorov-Smirnov Testsks.test distance between emperical and therotical CDF