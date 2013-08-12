# implemtent by the method based on "Acceptance Rejection"
gen_inhomo_poisson <- function(lambda, T_0, lambda_u = NULL) {
  result <- c()
  ## give bound of lambda(t), `lambda_u` should satisfy the condition:
  # 0 < lambda(t) <= lambda_u
  # `lambda_u` can be given and speed up the generator
  if (missing(lambda_u)) 
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