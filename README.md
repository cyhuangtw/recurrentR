recurrentR: non-parametric recurrent data analysis in R
=======================================================

Brief Introduction
==================

The package *recurrentR* implements three semi-parametric model of
recurrent data analysis in the following papers:

-   (M-C., Qin, and Chiang 2001)
-   (C.-Y. Huang and Wang 2004)
-   (C.-Y. Y. Huang, Qin, and Wang 2010)

TODO: Describe some background of non-parametric analysis of recurrent
event data here.

This technical note unifies the mathematical notations in these three
paper and describes the implemented mathematical formulas closely.

Data
====

Without loss of generality, we assume that there is a dataset of
recurrent event data containing $n$ instances. Each instance, $i$,
includes the following fields:

-   Censor time, $y_i$. TODO: explanation of censor time
-   Censor type, $D_i$. TODO:
-   Time period of observation: $[0, T_0]$.
-   Recurrent event time, $t_{i,1}$, $t_{i,2}$, ..., $t_{i,m_i}$. These
    are the realization of poisson process $N_i(.)$ TODO:
-   $q$-dim vector of time independent covariates, $W_i$. We assume that
    $W_i \in \mathbb{R}^{q \times 1}$. For simplicity, we denote
    $W \in \mathbb{R}^{q \times n}$ as the corresponding matrix. TODO:
-   $p$-dim time-dependent covariate process, $X_i(t)$. We assme that
    $X_i(t) \in \mathbb{R}^{p \times 1}$. TODO:

The field $X_i(t)$ is required for model of (C.-Y. Y. Huang, Qin, and
Wang 2010). The user could omit this field for model (M-C., Qin, and
Chiang 2001) and (C.-Y. Huang and Wang 2004).

In *recurrentR*, the data is stored in a S4-class object:
`recurrent-data`. The following is the structure of `recurrent-data`
with 100 instances, named `obj`:

    Formal class 'recurrent-data' [package "recurrentR"] with 6 slots
      ..@ W  : num [1:100, 1:2] 1 1 1 1 1 1 1 1 1 1 ...
      ..@ y  : num [1:100] 9.1729 8.8428 10 10 0.0597 ...
      ..@ t  :List of 100
      .. .. [list output truncated]
      ..@ X  :List of 100
      .. .. [list output truncated]
      ..@ T_0: num 10
      ..@ D  : logi [1:100] TRUE TRUE FALSE FALSE TRUE TRUE ...

The name of the slot is consistent to the variable name described above.
For example, for instance 1:

-   The censor time $y_1$ is `obj@y[1]`.
-   The censor type $D_1$ is `obj@D[1]`. `FALSE` stands for informative
    censoring(TODO: verify!!!).
-   The recurrent events $t_{1,1}, t_{1, 2}, ..., t_{1, m_1}$ is the
    numeric vector `obj@t[[1]]`.
-   The $T_0$ is `obj@T_0`
-   The $W_1$ is `obj@W[1,]`. And the $W \in \mathbb{R}^{q \times n}$ is
    `t(obj@W)`
-   The $X_1(t)$ is the function `obj@X[[1]]`

The user could create the object with the following function:

~~~~ {.r}
str(create_recurrent_data)
obj <- create_recurrent_data()
obj
~~~~

Usage
=====

(M-C., Qin, and Chiang 2001) {#section}
----------------------------

For each instance $i$, the occurrence of recurrent event follows a
inhomogenous poisson process with the following intensity:

$$\lambda_i(t) = \lambda_0(t) z_i exp(W_i \gamma)$$

where:

-   $z_i$ is a nonnegative-valued latent variable such that
    $E(z_i | W_i) = E(z_i)$.
-   The baseline intensity function $lambda_0(t)$ is a probability
    function:
    -   $\lambda_0(t) \neq 0$
    -   $\Lambda_0(T_0) = \int_0^{T_0} \lambda_0(u) du = 1$

-   $\gamma$ is a $\mathbb{R}^{1 \times q}$ vector.

In *recurrentR*:

~~~~ {.r}
library(recurrentR)
Wang2001(obj)
~~~~

(C.-Y. Huang and Wang 2004) {#section-1}
---------------------------

The intensity is the same:

$$\lambda_i(t) = \lambda_0(t) z_i exp(W_i \gamma)$$

where:

-   $z_i$ is a nonnegative-valued latent variable such that
    $E(z_i | W_i) = E(z_i)$.
-   The baseline intensity function $lambda_0(t)$ is a probability
    function:
    -   $\lambda_0(t) \neq 0$
    -   $\Lambda_0(T_0) = \int_0^{T_0} \lambda_0(u) du = 1$

-   $\gamma$ is a $\mathbb{R}^{1 \times q}$ vector.

Moreover, the hazard function of the censor time is modeled as

$$h_i(t) = h_0(t) z_i exp(W_i \alpha)$$

where:

-   $\alpha$ is a $\mathbb{R}^{1 \times q}$ vector.

Conditional on $(W_i, z_i)$, $N_i(.)$ and $y_i$ are independent.

~~~~ {.r}
library(recurrentR)
Huang2004(obj)
~~~~

(C.-Y. Y. Huang, Qin, and Wang 2010) {#section-2}
------------------------------------

The intensity is:

$$\lambda_i(t) = \lambda_0(t) z_i exp(X_i(t) \beta + \gamma W_i)$$

where:

-   $z_i$ is a nonnegative-valued latent variable such that
    $E(z_i | W_i) = E(z_i)$.
-   The baseline intensity function $lambda_0(t)$ is a probability
    function:
    -   $\lambda_0(t) \neq 0$
    -   $\Lambda_0(T_0) = \int_0^{T_0} \lambda_0(u) du = 1$

-   $\gamma$ is a $\mathbb{R}^{1 \times q}$ vector.
-   $\beta$ is a $\mathbb(R)^{1 \times p}$ vector.

Conditional on $(W_i, z_i, X_i)$, $N_i(.)$ and $y_i$ are independent.

~~~~ {.r}
library(recurrentR)
Huang2010(obj)
~~~~

Implementation Details
======================

(M-C., Qin, and Chiang 2001) {#section-3}
----------------------------

The inference are all included in the output of `Wang2001`.

~~~~ {.r}
library(recurrentR)
result <- Wang2001(obj)
~~~~

Recall that $$\lambda_i(t) = \lambda_0(t) z_i exp(W_i \gamma)$$ and
$$\Lambda_0(t) = \int_0^t \lambda_0(u) du$$. The nonparametric maximal
likelihood estimator $\hat{\Lambda}_0(t)$ is:

$$\hat{\Lambda}_0(t) = \prod_{s_{(l)} > t}(1 - \frac{d_{(l)}}{R_{(l)}})$$

where:

-   $s_{(l)}$ is the ordered and distinct values of event times
    ${t_{ij}}$.
-   $d_{(l)}$ is the number of events occurred at $s_{(l)}$.
-   $R_{(l)}$ is the total number of events ${t_{i,j}}$ which satisfies
    $t_{ij} \leq s_{(l)} \leq y_i$.

The user can obtain $\hat{\Lambda}_0(t)$:

~~~~ {.r}
str(result$Lambda.hat)
result$Lambda.hat(rexp(10))
~~~~

The $\hat{\gamma}$ is estimated by solving the following equation:

$$\frac{1}{n} \sum_{i=1}^n w_i \bar{W}_i^T ( \frac{m_i}{\hat{\Lambda}_0(y_i)} - exp(\bar{W}_i \bar{\gamma}) = 0 \in \mathbb{R}^{1 \times (q+1)},$$

where:

-   $\bar{W}_i = (1, W_i)$
-   $\bar{\gamma} = (\mu_Z, \gamma^T)^T$ and $E(z_i) = \mu_Z \forall i$.

(M-C., Qin, and Chiang 2001) provides the best $w_i$ to estimate
$\gamma$, but it involves $\gamma$ which might produce potential
instability. In *recurrentR*, we let $w_i = 1$. If the instance has
$\hat{\Lambda}_0(y_i) = 0$, we will let
$\frac{m_i}{\hat{\Lambda}_0(y_i)} = 0$ as usual convention.

Let
$$V(\bar{\gamma}) = \frac{1}{n} \sum_{i=1}^n {\bar{W}_i^T ( \frac{m_i}{\hat{\Lambda}_0(y_i)} - exp(\bar{W}_i \bar{\gamma})},$$

Then
$$\frac{dV}{d\bar{\gamma}}(\bar{\gamma}) = \frac{-1}{n} \sum_{i=1}^n{\bar{W}_i^T \bar{W}_i exp(\bar{W}_i \bar{\gamma})}$$

The *recurrentR* solves $V(\bar{\gamma}) = 0$ with the function
`nleqslv` from the package *nleqslv*. The $\hat{\bar{\gamma}}$ could be
accessed by:

~~~~ {.r}
result$gamma.bar.hat
~~~~

The $\hat{\gamma}$ is:

~~~~ {.r}
result$gamma.hat
~~~~

*recurrentR* provides both bootstrap and asymptotic variance estimator.
The user could choose one by passing the argument `method` in function
`Huang2001`.

~~~~ {.r}
result <- Wang2001(obj, method = "asymptotic")
result <- Wang2001(obj, method = "bootstrap")
result$gamma.bar.hat.var
result$gamma.hat.var
str(result$Lambda.hat.var)
~~~~

The default method is (TODO)

To calculate the asymptotic variance of $\hat{\gamma}$ and
$\hat{\Lambda}_0(t)$, we need the following formulas given
$\hat{\bar{\gamma}}$ and $\hat{\Lambda}_0(t)$:

-   $\hat{Q}(u) = \frac{1}{n} \sum_{i=1}^{n} { \sum_{j=1}^{m_i} {I(t_{i,j} \leq u)} }$
-   $\hat{R}(u) = \frac{1}{n} \sum_{i-1}^{n} { \sum_{j=1}^{m_i} {I(t_{i,j} \leq u \leq y_i)}}$
-   $\hat{b}_i(t) = \sum_{j=1}^{m_i}{ \int_{t}^{T_0} {\frac{I(t_{i,j} \leq u \leq y_i) d\hat{Q}(u)}{\hat{R}(u)^2}} - \frac{I(t \leq t_{i,j})}{\hat{R}(t_{i,j})}}$
    -   $\int_{t}^{T_0} {\frac{I(t_{i,j} \leq u \leq y_i) d\hat{Q}(u)}{\hat{R}(u)^2}} = \sum_{l} {\frac{I(t_{i,j} \leq s_{(l)} \leq y_i) d_{(l)} I(t \leq s_{(l)}) }{n \hat{R}(s_{(l)})^2}}$

-   $\hat{c}_i(t) = - \sum_{j=1}^n {\frac{m_j b_i(y_j)}{n \hat{\Lambda_0(y_j)}}} + \frac{m_i}{\hat{\Lambda}_0(y_i)} - \hat{\mu}_Z$
-   $\hat{d}_i(t) = \hat{\Lambda}_0(t) (\hat{c}_i + \hat{\mu}_Z \hat{b}_i(t) )$
-   $\hat{e}_i =\sum_{j=1}^n{ \frac{\bar{W}_j^T m_j b_i(y_j)}{n \hat{\Lambda_0}(y_j)}} + \bar{W}_i^T(\frac{m_i}{\hat{\Lambda}_0(y_i)} - exp(\bar{W}_i \hat{\bar{\gamma}}))$
    -   Let $\bar{\psi} = E[- \frac{d e_i}{d \bar{\gamma}}]$, then
        $\hat{\bar{\psi}} = \frac{1}{n} \sum_{i=1}^n{ \bar{W}_i^T \bar{W}_i exp(\bar{W}_i \hat{\bar{\gamma}})}$

-   $\hat{\bar{f}}_i(\hat{\bar{\gamma}}) = \hat{\bar{\psi}}^{-1} \hat{e}_i$
    -   Let
        $\hat{f}_i(\hat{\gamma}) = \hat{\bar{f}}_i(\hat{\bar{\gamma}})$
        without the first entry.

According to (M-C., Qin, and Chiang 2001), the asymptotic variacne of
$\hat{\gamma}$ is $\frac{1}{n} \sum_{i=1}^n{\hat{f}_i(\hat{\gamma})}$.
The asymptotic variance of $\hat{\Lambda}_0(t)$ is
$\hat{Lambda}_0(t)^2 * \frac{\sum_{i}{b_i(t)^2}}{n^3}$

(C.-Y. Huang and Wang 2004) {#section-4}
---------------------------

The estimator and asymptotic variance related to $\Lambda_0$ and
$\gamma$ are the same as the one in (M-C., Qin, and Chiang 2001). To
obtain the estimator of $\alpha$ and $H_0(t) = \int_0^t h(u) du$, we
need the estimator of random effect $z_i$ first:

$$\hat{z}_i = \frac{m_i}{\hat{\Lambda}_0(y_i) exp(W_i \hat{\gamma)}}.$$

Let
$$U(\alpha) = \frac{1}{n} \sum_{i=1}^n {D_i W_i^T \frac{\sum_{j=1}^n{W_j^T \hat{z}_j exp(W_j \hat{\alpha}) I(y_j \geq y_i)}}{\sum_{j=1}^n{\hat{z}_j exp(W_j \hat{\alpha})I(y_j \geq y_i)}} },$$

Then $\hat{\alpha}$ is the one satisfies $U(\hat{\alpha}) = 0$.

Moreover, Let
$$\Gamma(\alpha) = \frac{dU}{d\alpha}(\alpha) = \frac{1}{n} \sum_{i=1}^n{D_i(-\frac{\sum_{j=1}^n{W_j^2 \hat{z}_j exp( W_j \alpha ) I(y_j \geq y_i)}}{\sum_{j=1}^n{\hat{z}_j exp( W_j \alpha ) I(y_j \geq y_i) }} + \frac{(\sum_{j=1}^n{W_j \hat{z}_j exp( W_j \alpha ) I(y_j \geq y_i) })^2}{(\sum_{j=1}^n{\hat{z}_j exp( W_j \alpha ) I(y_j \geq y_i)})^2})},$$

Then we can solve $\hat{\alpha}$ with `nleqslv` again. Note that $a^2$
is the convention of $a^T a$ if $a$ is a vector.

With $\hat{\alpha}$, the $\hat{H}_0(t)$ will be:

$$\hat{H}_0(t) = \sum_{i=1}^n{D_i I(y_i \leq t) \frac{1}{\sum_{j=1}^n{\hat{z}_j exp(W_j \alpha) I(y_j \geq y_i)}}}.$$

~~~~ {.r}
result <- Huang2004(obj)
result$alpha.hat
str(result$H0.hat)
~~~~

To evaluate the asymptotic variance, we need:

-   $\psi_{3i}(t, \alpha) = \frac{1}{n}\sum_{j=1}^n{\frac{m_j}{\hat{Lambda}_0(y_j)} exp(W_j(\alpha - \gamma)) I(y_j \geq t) (W_j \hat{f}_i(\alpha) + b_i(y_j))} + \frac{m_i}{\hat{Lambda}_0(y_i)} exp(W_i(\alpha - \gamma)) I(y_i \geq t) - \frac{1}{n}\sum_{j=1}^n{\hat{z}_j exp(W_j \alpha) I(y_j \geq t)}$
-   $\psi_{4i}(t, \alpha) = \frac{1}{n}\sum_{j=1}^n{W_j \frac{m_j}{\hat{Lambda}_0(y_j)} exp(W_j(\alpha - \gamma)) I(y_j \geq t) (W_j \hat{f}_i(\alpha) + b_i(y_j))} + W_i \frac{m_i}{\hat{Lambda}_0(y_i)} exp(W_i(\alpha - \gamma)) I(y_i \geq t) - \frac{1}{n}\sum_{j=1}^n{W_j \hat{z}_j exp(W_j \alpha) I(y_j \geq t)}$
-   $\psi_i(\alpha) = W_i D_i - n^{-1}\sum_{j=1}^{n}{W_j D_j} + \sum_{j=1}^n{D_j \psi_{3i}(y_j, \alpha) \frac{\sum_{k=1}^n{W_k \hat{z}_k exp(W_k \alpha) I(y_k \geq y_j)}}{(\sum_{k=1}^n{\hat{z}_k exp(W_k \alpha) I(y_k \geq y_j)})^2}} - \sum_{j=1}^n{D_j \frac{\psi_{4i}(y_j, \alpha)}{\sum_{k=1}^n{\hat{z}_k exp(W_k \alpha) I(y_k \geq y_j)}}} + \frac{1}{n} \sum_{j=1}^{n}{D_j \frac{\sum_{k=1}^n{W_k \hat{z}_k exp(W_k \alpha) I(y_k \geq y_j)}}{\sum_{k=1}^n{\hat{z}_k exp(W_k \alpha) I(y_k \geq y_j)}}} - D_i \frac{\sum_{k=1}^n{\hat{z}_k exp(W_k \alpha) I(y_k \geq y_i)}}{\sum_{k=1}^n{\hat{z}_k exp(W_k \alpha) I(y_k \geq y_i)}}$
-   $\psi^*(\alpha) = \frac{1}{n} \sum_{i=1}^n \psi_i(\alpha)$
-   $\hat{\Sigma}(\alpha) = n^{-1} \sum_{i=1}^n{(\psi_i(\alpha) - \psi^*(\alpha))(\psi_i(\alpha) - \psi^*(\alpha))^T}$

According to (C.-Y. Huang and Wang 2004), the estimator of asymptotic
variance of $\alpha$ will be:

$$\frac{1}{n} \Gamma(\hat{\alpha})^{-1} \hat{\Sigma}(\hat{\alpha}) \Gamma(\hat{\alpha})^{-1}.$$

For the asymptotic variance of $\hat{H}_0(t)$, we need

-   $\phi_i(t) = n \sum_{j=1}^n{D_j I(y_j \leq t) \frac{\psi_{3i}(y_j, \hat{\alpha})}{(\sum_{k=1}^n{\hat{z}_k exp(W_k \alpha) I(y_k \geq y_j))^2}}} - \sum_{j=1}^n{D_j I(y_j \leq t) \frac{1}{\sum_{k=1}^n{ \hat{z}_k exp( W_k \alpha) I(y_k \geq y_j) }}} + D_i I(y_i \leq t) \frac{n}{\sum_{k=1}^n{\hat{z}_k exp(W_k \alpha) I(y_k \geq y_i)}} - \frac{d \hat{H}_0}{d \alpha}(t, \hat{\alpha}) \Gamma(\hat{\alpha})^{-1}\psi_i(\hat{\alpha})$
    -   $\frac{d\hat{H}_0}{d\alpha}(t, \alpha) = -\sum_{i=1}^n{D_i I(y_i \leq t) \frac{\sum_{j=1}^n{W_j \hat{z}_j exp(W_j \alpha) I(y_j \geq y_i)}}{(\sum_{j=1}^n{\hat{z}_j exp(W_j \alpha) I(y_j \geq y_i)})^2}}$

Then the estimator of asymptotic variance of $\hat{H}_0(t)$ is the
sample variance of $\frac{1}{n} \phi_i(t)$.

Huang2010
---------

Recall that the intensity is:

$$\lambda_i(t) = z_i \lambda_0(t) exp(X_i(t) \beta + W_i \gamma)$$

The estimator of $\hat{\beta}$ does not involve $W_i$ and $\gamma$.

The derivative of logged pairwise pseudolikelihood is:

$$g_{i,j}(\beta) = \sum_{k=1}^{m_i}{ \sum_{l=1}^{m_j}{ I(t_{i,k} \leq y_{i,j}) I(t_{j,l} \leq y_{i,j}) \rho_{i,j}(t_{i,k}, t_{j,l}) \frac{- exp(\rho_{i,j}(t_{i,k}, t_{j,l}) \beta)}{1 + exp(\rho_{i,j}(t_{i,k}, t_{j,l}) \beta)} } },$$

where

-   $\rho_{i,j}(u, v) = X_i(v) + X_j(u) - X_i(u) - X_j(v)$

Let
$$S(\beta) = \frac{1}{\left(\begin{array}{c} n \\ 2 \end{array}\right)} \sum_{i < j}{g_{i,j}(\beta)},$$

Then
$$\frac{dS}{d\beta}(\beta) = \frac{1}{\left(\begin{array}{c} n \\ 2 \end{array}\right)} \sum_{i<j} {\frac{dg_{i,j}}{d\beta}(\beta)} ,$$

where:

-   $\frac{dg_{i,j}}{d\beta}(\beta) = \sum_{k=1}^{m_i}{ \sum_{l=1}^{m_j}{ I(t_{i,k} \leq y_{i,j}) I(t_{j,l} \leq y_{i,j}) \frac{- \rho_{i,j}(t_{i,k}, t_{j,l})^2 exp(\rho_{i,j}(t_{i,k}, t_{j,l}) \beta)}{(1 + exp(\rho_{i,j}(t_{i,k}, t_{j,l}) \beta))^2}  } }$

The $\hat{\beta}$ is the one satisfies $S(\beta) = 0$.

To obtain the asymptotic variance, we need:

-   $\hat{V_1} = \frac{4}{n}\sum_{i=1}^n{ \frac{1}{\left(\begin{array}{c}n-1 \\ 2 \end{array}\right)} \sum_{i < j < k}{g_{i,j}(\hat{\beta}) g_{i,k}(\hat{\beta})} }$
-   $\hat{V_2} = \frac{-1}{\left(\begin{array}{c}n \\ 2 \end{array}\right)} \sum_{i < k}{\frac{dg_{i,k}}{d\beta}(\hat{\beta})}$

Recall that in (M-C., Qin, and Chiang 2001), the $\hat{\Lambda}_0(t)$ is
based on:

$$\hat{\Lambda}_0(t) = \prod_{s_{(l)} > t}(1 - \frac{d_{(l)}}{R_{(l)}})$$

where:

-   $s_{(l)}$ is the ordered and distinct values of event times
    ${t_{ij}}$.
-   $d_{(l)}$ is the number of events occurred at $s_{(l)}$.
-   $R_{(l)}$ is the total number of events ${t_{i,j}}$ which satisfies
    $t_{ij} \leq s_{(l)} \leq y_i$.

To correct the effect of time-dependent covariates $X(t)$ and $\beta$,
we let

$$\hat{\Lambda}_0(t, \beta) = \prod_{s_{(l)} > t}(1 - \frac{d_{(l)}(\beta)}{R_{(l)}(\beta)}),$$

where:

-   $s_{(l)}$ is the ordered and distinct values of event times
    ${t_{ij}}$.
-   $d_{(l)}(\beta) = \frac{1}{n} \sum_{i=1}^n { \sum_{j=1}^{m_i} { I(t_{i,j} == s_{(l)}) exp(-X_i(t_{i,j]} \beta)) } }$
-   $R_{(l)}(\beta) = \frac{1}{n} \sum_{i=1}^n { \sum_{j=1}^{m_i} { I(t_{i,j} \leq s_{(l)} \leq y_i) exp( -X_i(t_{i,j}) \beta ) } }$

Note that $d_{(l)}(0) = d_{(l)}$ and $R_{(l)}(0) = R_{(l)}$.

The asymptotic variance of $\hat{\Lambda}_0(t, \beta)$ will be

$$4 \Lambda_0(t)^2 E\left[\kappa_{1,2} (t, \beta) \kappa_{1, 3}( t, \beta)\right].$$

The estimator of $\hat{\gamma}$ is the root of the following equation:

$$\frac{1}{n} \sum_{i=1}^n{\bar{W}_i \left[\frac{m_i}{\sum_{s_{(l)} \leq y_i}{exp(X_i(s_{(l)}) \hat{\beta})(\hat{\Lambda}_0(s_{(l)}) - \hat{\Lambda}_0(s_{(l-1)})})} - exp(\bar{W}_i \bar{\gamma})\right]},$$

where:

-   $\bar{W}_i = (1, W_i)$
-   $\bar{\gamma} = (\mu_Z, \gamma^T)$

The asymptotic variance of $\hat{\gamma}$ will be:

$$E(\frac{d\xi}{d\gamma})^{-1} \Sigma E(\frac{d\xi}{d\gamma})^{-1}$$.

To obtain the definition of $\kappa$, we need:

-   $\hat{Q}(u) = \frac{1}{n} \sum_{i=1}^n{ \sum_{j=1}^{m_i} { I(t_{i,j} \leq u) exp(- X_i(t_{i,j}) \hat{\beta}) } }$
-   $\hat{R}(u) = \frac{1}{n} \sum_{i=1}^n{ \sum_{j=1}^{m_i} { I(t_{i,j} \leq u \leq y_i) exp(- X_i(t_{i,j}) \hat{\beta}) } }$
-   $\hat{V}_{\tilde{Q}}(u) = E \frac{dQ}{d\beta} = \frac{1}{n} \sum_{i=1}^n{ \sum_{j=1}^{m_i} { - X_i(t_{i,j}) I(t_{i,j} \leq u) exp(- X_i(t_{i,j}) \hat{\beta}) } }$
-   $\hat{V}_{\tilde{R}}(u) = E \frac{dR}{d\beta} = \frac{1}{n} \sum_{i=1}^n{ \sum_{j=1}^{m_i} { - X_i(t_{i,j}) I(t_{i,j} \leq u \leq y_i) exp(- X_i(t_{i,j}) \hat{\beta}) } }$
-   $\hat{\phi}_{i,j}(t) = \left(\int_t^{T_0} { \frac{d\hat{V}_{\tilde{Q}}(u)}{\hat{R}(u)} - \frac{\hat{V}_{\tilde{R}}(u) d\hat{Q}(u)}{\hat{R}(u)^2} }\right)\hat{V}_2^{-1} g_{i,j}(\hat{\beta}) = \left(\sum_{t \leq s_{(l)}}{\frac{\hat{V}_{\tilde{Q}}(s_{(l)}) - \hat{V}_{\tilde{Q}}(s_{(l-1)})}{\hat{R}(s_{(l)})}} - \sum_{t \leq s_{(l)}}{\frac{\hat{V}_{\tilde{R}}(s_{(l)})(\hat{Q}(s_{(l)}) - \hat{Q}(s_{(l-1)}))}{\hat{R}(s_{(l)})^2}} \right)\hat{V}_2^{-1} g_{i,j}(\hat{\beta})$
-   $\hat{\psi}_i(t) = \sum_{j=1}^{m_i}{I( t < t_{i,j}) \frac{1}{\hat{R}(t_{i,j})}} - \int_t^{T_0} {I(t_{i,j} \leq u \leq y_i) exp(X_i(t_{i,j}) \hat{\beta}) \frac{d \hat{Q}(u)}{\hat{R}(u)^2} } = \sum_{j=1}^{m_i}{I( t < t_{i,j}) \frac{1}{\hat{R}(t_{i,j})}} - \left( \sum_{s_{(l)} > t} { I(t_{i,j} \leq s_{(l)} \leq y_i) exp(-X_i(t_{i,j}) \hat{\beta}) \frac{(\hat{Q}(s_{(l)}) - \hat{Q}(s_{(l-1)}))}{\hat{R}(s_{(l)})^2} } \right)$

Then

$$\hat{\kappa}_{i,j}(t) = \hat{\phi}_{i,j}(t) + \frac{\hat{\psi}_i(t) + \hat{\psi}_j(t)}{2}$$

The estimator of asymptotic variance of $\hat{\Lambda}_0(t)$ will be

$$4 \hat{\Lambda_0}(t)^2 \left(\frac{2}{n(n-1)(n-2)} \sum_{i=1}^n{ \sum_{j \neq i, k \neq i, j < k}{ \hat{\kappa}_{i,j}(t)\hat{\kappa}_{i,k}(t) } }\right)$$

The asymptotic variance of $\hat{\gamma}$ involves:

-   $\hat{\xi}_{i,j} = \frac{1}{n} \sum_{k=1}^n{ \frac{- \bar{W}_k m_k}{\left[\sum_{s_{(l)} \leq y_k} { exp(X_k(s_{(l)}) \hat{\beta}) (\hat{\Lambda}_0(s_{(l)}) - \hat{\Lambda}_0(s_{(l-1)})) }\right]^2} \left[\sum_{s_{(l)} \leq y_k} { X_k(s_{(l)}) \hat{V}_2^{-1} g_{i,j}(\hat{\beta}) exp(X_k(s_{(l)}) \hat{\beta}) ( \hat{\Lambda}_0(s_{(l)}) - \hat{\Lambda}_0(s_{(l-1)})) + exp(X_k(s_{(l)}) \hat{\beta}) (\hat{\kappa}_{i,j}(s_{(l)}) \hat{\Lambda}_0(s_{(l)}) - \hat{\kappa}_{i,j}(s_{(l-1)}) \hat{\Lambda}_0(s_{(l-1)})) } \right]  } + \frac{1}{2} \bar{W}_i \left[ \frac{m_i}{\sum_{s_{(l)} \leq y_i} { exp(X_i(s_{(l)}) \hat{\beta}) (\hat{\Lambda}_0(s_{(l)}) - \hat{\Lambda}_0(s_{(l-1)})) }} - exp(\bar{W}_i \hat{\bar{\gamma}}) \right] + \frac{1}{2} \bar{W}_j \left[ \frac{m_j}{\sum_{s_{(l)} \leq y_j} { exp(X_j(s_{(l)}) \hat{\beta}) (\hat{\Lambda}_0(s_{(l)}) - \hat{\Lambda}_0(s_{(l-1)})) }} - exp(\bar{W}_j \hat{\bar{\gamma}}) \right]$
-   $- \frac{d\xi_{i,j}}{d\gamma} = \frac{1}{2} \bar{W}_i^2 exp(\bar{W}_i \hat{\bar{\gamma}}) + \frac{1}{2} \bar{W}_j^2 exp(\bar{W}_j \hat{\bar{\gamma}})$

Therefore, the estimator of asymptotic variance of $\hat{\bar{\gamma}}$
is:

$$TODO$$

Reference
=========

Huang, C.-Y. Y., J. Qin, and M.-C. C. Wang. 2010. “Semiparametric
analysis for recurrent event data with time-dependent covariates and
informative censoring.” *Biometrics* 66 (1) (mar 12): 39–49.
doi:10.1111/j.1541-0420.2009.01266.x.
<http://dx.doi.org/10.1111/j.1541-0420.2009.01266.x>.

Huang, Chiung-Yu, and Mei-Cheng Wang. 2004. “Joint Modeling and
Estimation for Recurrent Event Processes and Failure Time Data.”
*Journal of the American Statistical Association* 99: 1153–1165.
<http://EconPapers.repec.org/RePEc:bes:jnlasa:v:99:y:2004:p:1153-1165>.

M-C., Wang, J. Qin, and C.-T. Chiang. 2001. “Analyzing Recurrent Event
Data With Informative Censoring.” *Journal of the American Statistical
Association* 96: 1057–1065.
<http://EconPapers.repec.org/RePEc:bes:jnlasa:v:96:y:2001:m:september:p:1057-1065>.
