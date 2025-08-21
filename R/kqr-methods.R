#' Extract model coefficients from a `kqr` object.
#'
#' Computes the coefficients at the requested value(s) for `lambda` from a
#' [kqr()] object.
#'
#' `s` is the new vector of `lambda` values at which predictions are requested.
#' If `s` is not in the lambda sequence used for fitting the model, the `coef`
#' function will use linear interpolation to make predictions. The new values
#' are interpolated using a fraction of coefficients from both left and right
#' `lambda` indices.
#' @importFrom methods rbind2
#' @param object Fitted [kqr()] object.
#' @param s Value(s) of the penalty parameter `lambda` at which
#'  coefficients are required. Default is the entire sequence.
#' @param ... Not used.
#' @seealso [kqr()] and [predict.kqr()].
#'
#' @return The coefficients at the requested values for `lambda`.
#'
#' @method coef kqr
#' @export
#' @examples
#' library(MASS)
#' data(GAGurine)
#' x <- as.matrix(GAGurine$Age)
#' y <- GAGurine$GAG
#' lambda <- 10^(seq(1, -4, length.out=10))
#' fit <- kqr(x, y, lambda=lambda, tau=0.1)
#' coef(fit)

coef.kqr <- function(object, s = NULL, ...) {
  # rlang::check_dots_empty()
  b0 <- matrix(object$alpha[1,], nrow = 1)
  rownames(b0) <- "(Intercept)"
  alpha <- rbind2(b0, object$alpha[-1,,drop=FALSE])
  if (!is.null(s)) {
    vnames <- dimnames(alpha)[[1]]
    dimnames(alpha) <- list(NULL, NULL)
    lambda <- object$lambda
    lamlist <- lambda.interp(lambda, s)
    ls <- length(s)
    if (ls == 1) {
      alpha = alpha[, lamlist$left, drop = FALSE] * lamlist$frac +
        alpha[, lamlist$right, drop = FALSE] * (1 - lamlist$frac)
    } else {
      alpha = alpha[, lamlist$left, drop = FALSE] %*%
        Matrix::Diagonal(ls, lamlist$frac) +
        alpha[, lamlist$right, drop = FALSE] %*%
        Matrix::Diagonal(ls, 1 - lamlist$frac)
    }
    if (is.null(names(s))) names(s) <- paste0("s", seq_along(s))
    dimnames(alpha) <- list(vnames, names(s))
  }
  return(alpha)
}



#' Predict the fitted values for a \code{kqr} object.
#'
#' @importFrom methods cbind2
#' @importFrom stats coef 
#' @param object A fitted \code{kqr} object.
#' @param x The predictor matrix, i.e., the \code{x} matrix used when fitting the \code{kqr} object.
#' @param newx A matrix of new values for \code{x} at which predictions are to be made. Note
#' that \code{newx} must be of a matrix form, predict function does not accept a vector or other
#' formats of \code{newx}.
#' @param s Value(s) of the penalty parameter `lambda` at which
#'   predictions are required. Default is the entire sequence used to create the
#'   model.
#' @param ... Not used.
#'
#' @details
#' The result is \eqn{\beta_0 + K_i' \alpha} where \eqn{\beta_0} and \eqn{\alpha} are from the
#' \code{kqr} object and \eqn{K_i} is the ith row of the kernel matrix.
#'
#' @return
#' Returns the fitted values.
#' @keywords classification kernel
#' @useDynLib fastkqr, .registration=TRUE
#' @method predict kqr
#' @export
#' @examples
#' library(MASS)
#' data(GAGurine)
#' x <- as.matrix(GAGurine$Age)
#' y <- GAGurine$GAG
#' lambda <- 10^(seq(1, -4, length.out=30))
#' fit <- kqr(x, y, lambda=lambda, tau=0.1, is_exact=TRUE)
#' predict(fit, x, tail(x))

predict.kqr <- function(object, x, newx=NULL, s=NULL,...) {
  sigma <- object$sigma
  alpha <- coef(object, s)
  newK <- kernelMat(newx, x, sigma=sigma)
  fit <- as.matrix(cbind2(1, newK) %*% alpha)
  fit
}

