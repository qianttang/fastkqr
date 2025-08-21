## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(fastkqr)
library(MASS)
data(GAGurine)
x <- as.matrix(GAGurine$Age)
y <- GAGurine$GAG

## ------------------------------------------------------------------------
lambda <- 10^(seq(1, -4, length.out=10))
fit <- kqr(x, y, lambda=lambda, tau=0.1, is_exact=TRUE)

## ------------------------------------------------------------------------
cv.fit <- cv.kqr(x, y, lambda=lambda, tau=0.1)

## ------------------------------------------------------------------------
coef <- coef(fit, s = c(0.02, 0.03))
predict(fit, x, tail(x), s = fit$lambda[2:3])

## ------------------------------------------------------------------------
l2 <- 1e-4
tau <- c(0.1, 0.3, 0.5)
l1_list <- 10^seq(-8, 2, length.out=10)
fit1 <- nckqr(x ,y, lambda1 = l1_list, lambda2 = l2,  tau = tau)

## ------------------------------------------------------------------------
l2_list <- 10^(seq(1, -4, length.out=10))
cv.fit1 <- cv.nckqr(x, y, lambda1=10, lambda2=l2_list, tau=tau)

## ------------------------------------------------------------------------
coef <- coef(fit1, s2=1e-4, s1 = l1_list[2:3])
predict(fit1, x, tail(x), s1=l1_list[1:3], s2=l2)

