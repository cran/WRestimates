################################################################################
# Sample size calculation
## Author: Autumn O Donnell
## Date: 02/08/2023
## R version 4.1.2 (2021-11-01) - Bird Hippie
################################################################################
################################################################################
# The following functions are derived from Ron Xiaolong Yu and Jitendra Ganju
# 2022 paper "Sample size formula for a win ratio endpoint"
# DOI: 10.1002/sim.9297

# Let no.W denote the number of wins, no.L the number of losses, and no.T the
# number of ties
# Further, let N denote the total sample size
# k denote the proportion of patients allocated to one group (so (1 - k) is the
# proportion allocated to the other group)
# let p.tie = no.T/(k(k - 1)N^2) represent the proportion of ties
# Then, k(1 - k)N^2 = no.W + no.L + no.T where k(1 - k)N^2 is the product of the
# sample size per group.

# The approximate variance of the log of the win ratio, Var(ln(WR)) or
# var.ln.WR, is (1/N)*sigma.sqr; where
# sigma.sqr = (4(1 + p.tie)/3k(1 - k)(1 - p.tie))

# The required sample size N is approximately equal to
# (sigma.sqr * (Z.1minusalpha + Z.1minusbeta)^2)/(ln^2(WR.true)) where; WR.true
# is the assumed or true win ration and ln(WR.true) is its natural log-transform
# alpha and beta refer, respectively, to type I and type II error rates, with Z
# denoting the quantile value from the standard normal distribution

################################################################################
################################################################################
# Functions

## Sigma^2 function
wr.sigma.sqr <- function(k, p.tie){
  if (p.tie >=0 & p.tie <= 1){
    sigma.sqr = (4*(1 + p.tie))/(3*k*(1 - k)*(1 - p.tie))

    return(list(sigma.sqr = sigma.sqr,
                k = k,
                p.tie = p.tie))
  } else {
    return(stop("p.tie is not in range"))
  }
}

## Sample Size function
wr.ss <- function(alpha = 0.025,
                        beta = 0.1,
                        WR.true = 1,
                        k,
                        p.tie,
                        sigma.sqr){
  if(missing(sigma.sqr)){
    sigma.sqr <- wr.sigma.sqr(k, p.tie)$sigma.sqr
  }
  if ((alpha < 0 | alpha > 1) | beta < 0 | beta > 1){
    return(stop("Alpha or beta value out of range."))
    } else {
      if (!is.null(k) & !is.null(p.tie) & !is.null(sigma.sqr)){
        check <- wr.sigma.sqr(k, p.tie)$sigma.sqr
        if (check != sigma.sqr){
          return(warning("sigma.sqr is incorrect given the k and p.tie values provided."))
        }
      }
      Z.1minusalpha <- qnorm(1 - alpha)
      Z.1minusbeta <- qnorm(1 - beta)

      N = (sigma.sqr*(Z.1minusalpha + Z.1minusbeta)^2)/((log(WR.true))^2)

      if (!is.null(k) & !is.null(p.tie)){
        return(list(N = ceiling(N),
                    alpha = alpha,
                    beta = beta,
                    WR.true = WR.true,
                    sigma.sqr = sigma.sqr,
                    k = k,
                    p.tie = p.tie))
        } else {
          return(list(N = ceiling(N),
                      alpha = alpha,
                      beta = beta,
                      WR.true = WR.true,
                      sigma.sqr = sigma.sqr))
        }
    }
  }

# Power function
wr.power <- function(N,
                     alpha = 0.025,
                     WR.true = 1,
                     sigma.sqr,
                     k,
                     p.tie) {
  if (WR.true < 1){
    return(stop("WR.true is not in range."))
  }

  if (alpha < 0 | alpha > 1){
    return(stop("Alpha value out of range."))
  } else {
    if (missing(sigma.sqr)){
      sigma.sqr <- wr.sigma.sqr(k, p.tie)$sigma.sqr
      }
    Z.1minusalpha <- qnorm(1 - alpha)

    power = 1 - (pnorm((Z.1minusalpha - (log(WR.true)*(sqrt(N/sigma.sqr))))))

    if (!is.null(k) & !is.null(p.tie)){
      return(list(power = power,
                  N = N,
                  alpha = alpha,
                  WR.true = WR.true,
                  sigma.sqr = sigma.sqr,
                  k = k,
                  p.tie = p.tie))
      } else {
        return(list(power = power,
                    N = N,
                    alpha = alpha,
                    WR.true = WR.true,
                    sigma.sqr = sigma.sqr))
      }
  }
}

## Variance function
wr.var <- function(N, sigma.sqr, k, p.tie){

  if (!is.null(k) & !is.null(p.tie) & !is.null(sigma.sqr)){
    check <- wr.sigma.sqr(k, p.tie)$sigma.sqr
    if (check != sigma.sqr){
      return(warning("sigma.sqr is incorrect given the k and p.tie values provided."))
    }
  }
  if (missing(sigma.sqr)){
    if (missing(k) | missing(p.tie)){
      return(stop("Cannot calculate sigma.sqr"))
    }
    else {
      sigma.sqr <- wr.sigma.sqr(k, p.tie)$sigma.sqr
    }
  }

  var.ln.WR = sigma.sqr/N

  if (!is.null(k) & !is.null(p.tie)){
    return(list(var.ln.WR = var.ln.WR,
                N = N,
                sigma.sqr = sigma.sqr,
                k = k,
                p.tie = p.tie))
  }
  else {
    return(list(var.ln.WR = var.ln.WR,
                N = N,
                sigma.sqr = sigma.sqr))
  }
}

## CI function
wr.ci <- function(WR = 1,
                        Z = 1.96,
                        var.ln.WR,
                        N,
                        sigma.sqr,
                        k,
                        p.tie){
  if (missing(var.ln.WR)){
    if (missing(sigma.sqr)){
      if (missing(k) | missing(p.tie)){
        return(stop("Cannot calculate sigma.sqr"))
      }
      else {
        sigma.sqr <- wr.sigma.sqr(k, p.tie)$sigma.sqr
      }
    }
    var.ln.WR <- wr.var(N, sigma.sqr, k, p.tie)$var.ln.WR
  }

  ci = c(exp(log(WR) - Z*(sqrt(var.ln.WR))),
    exp(log(WR) + Z*(sqrt(var.ln.WR))))
  names(ci) <- c("lower.CI", "upper.CI")

  if (!is.null(N) & !is.null(sigma.sqr) & !is.null(k) & !is.null(p.tie)){
    return(list(ci = ci,
                WR = WR,
                Z = Z,
                var.ln.WR = var.ln.WR,
                N = N,
                sigma.sqr = sigma.sqr,
                k = k,
                p.tie = p.tie))
  }
  else {
    return(list(ci = ci,
                WR = WR,
                Z = Z,
                var.ln.WR = var.ln.WR))
  }
}

################################################################################
################################################################################
## Example 1 (with functions)
### 1:1 allocation, one-sided alpha = 2.5%, power = 90% (beta = 10%), a small
### proportion of ties p.tie = 0.1, and 50% more wins on treatment than control

## Calculate Sample Size
# wr.ss(WR.true = 1.5, k = 0.5, p.tie = 0.1)

### Calculate the Power (N = 117)
# wr.power(N = 417, WR.true = 1.5, k = 0.5, p.tie = 0.1)

### Calculation 95% CI (N = 117)
# wr.ci(N = 417, WR = 1.5, k = 0.5, p.tie = 0.1)


################################################################################
## Example 1 (without functions)
### 1:1 allocation, one-sided alpha = 2.5%, power = 90% (beta = 10%), a small
### proportion of ties p.tie = 0.1, and 50% more wins on treatment than control

### Set-up
# k <- 0.5
# alpha <- 0.025
# beta <- 0.1
# Z.1minusalpha <- qnorm(1 - alpha)
# Z.1minusbeta <- qnorm(1 - beta)
# p.tie <- 0.1
# WR.true <- 1.5
# WR <- 1.5

### Calculate sigma^2
# sigma.sqr = (4*(1 + p.tie))/(3*k*(1 - k)*(1 - p.tie))

### Calculate N
# N = (sigma.sqr*(Z.1minusalpha + Z.1minusbeta)^2)/(log^2(WR.true))

### Calculation Var(ln(WR))
# var.ln.WR = sigma.sqr/N

### Calculate the Power
# Phi <- pnorm(1.96)
# Z.alpha <- qnorm(alpha, lower.tail = FALSE)
# Power = 1 - (Phi*(Z.alpha - (log2(WR.true)*sqrt(N/sigma.sqr))))

### Calculation 95% CI
# CI.95 = c(exp(log2(WR) - 1.96*(sqrt(var.ln.WR))),
#           exp(log2(WR) + 1.96*(sqrt(var.ln.WR))))
# names(CI.95) <- c("lower.CI", "upper.CI")
