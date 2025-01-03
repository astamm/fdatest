#' Two population Interval Wise Testing procedure
#'
#' The function implements the Interval Wise Testing procedure for testing mean
#' differences between two functional populations. Functional data are tested
#' locally and unadjusted and adjusted p-value functions are provided. The
#' unadjusted p-value function controls the point-wise error rate. The adjusted
#' p-value function controls the interval-wise error rate.
#'
#' @param data1 First population's data. Either pointwise evaluations of the
#'   functional data set on a uniform grid, or an \code{\link[fda]{fd}} object.
#'   If pointwise evaluations are provided, it should be a matrix of dimensions
#'   `c(n1, J)`, with `J` evaluations on columns and `n1` units on rows.
#' @param data2 Second population's data. Either pointwise evaluations of the
#'   functional data set on a uniform grid, or an \code{\link[fda]{fd}} object.
#'   If pointwise evaluations are provided, it should be a matrix of dimensions
#'   `c(n2, J)`, with `J` evaluations on columns and `n2` units on rows.
#' @param mu Functional mean difference under the null hypothesis. Three
#'   possibilities are available for \code{mu}: 
#'   
#'   - a constant (in this case, a constant function is used);
#'   - a \code{J}-dimensional vector containing the evaluations on the same grid
#'   which \code{data} are evaluated;
#'   - a \code{fd} object from the package \code{fda} containing one function.
#'   
#'   Defaults to `mu = 0`.
#' @inheritParams IWT1
#' @param paired Flag indicating whether a paired test has to be performed.
#'   Defaults to `FALSE`.
#' @param alternative A character string specifying the alternative hypothesis.
#'   Must be one of `"two.sided"` (default), `"greater"` or `"less"`.
#' @param verbose Logical: if \code{FALSE}, reduces the amount of output.
#'   Default is \code{TRUE}.
#'
#' @return An object of class \code{\link{IWT2}}, which is a list containing at
#'   least the following components:
#' - `test`: String vector indicating the type of test performed. In this case
#' equal to `"2pop"`.
#' - `mu`: Evaluation on a grid of the functional mean difference under the null
#' hypothesis (as entered by the user).
#' - `unadjusted_pval`: Evaluation on a grid of the unadjusted p-value function.
#' - `pval_matrix`: Matrix of dimensions `c(p, p)` of the p-values of the
#' interval-wise tests. The element `(i, j)` of matrix `pval.matrix` contains
#' the p-value of the test contains the p-value of the test of interval indexed
#' by `(j,j+1,...,j+(p-i))`.
#' - `adjusted_pval`: Evaluation on a grid of the adjusted p-value function.
#' - `data.eval`: Evaluation on a grid of the functional data.
#' - `ord_labels`: Vector of labels indicating the group membership of
#' `data.eval`.
#'
#' @seealso See also \code{\link{plot.fdatest2}} and \code{\link{IWTimage}} for
#'   plotting the results.
#'
#' @references 
#' A. Pini and S. Vantini (2017). The Interval Testing Procedure: Inference
#' for Functional Data Controlling the Family Wise Error Rate on Intervals.
#' *Biometrics*, 73(3): 835–845.
#' 
#' A. Pini and S. Vantini (2017). Interval-wise testing for functional data.
#' *Journal of Nonparametric Statistics*, 29(2), 407-424.
#'
#' @export
#' @examples
#' # Performing the IWT for two populations
#' IWT.result <- IWT2(NASAtemp$paris, NASAtemp$milan, B = 10L)
#'
#' # Plotting the results of the IWT
#' plot(
#'   IWT.result, 
#'   xrange = c(0, 12), 
#'   main = 'IWT results for testing mean differences'
#' )
#'
#' # Plotting the p-value heatmap
#' IWTimage(IWT.result, abscissa_range = c(0, 12))
#'
#' # Selecting the significant components at 5% level
#' which(IWT.result$adjusted_pval < 0.05)
IWT2 <- function(data1, data2, 
                 mu = 0, 
                 B = 1000L, 
                 dx = NULL, 
                 recycle = TRUE, 
                 paired = FALSE, 
                 alternative = "two.sided", 
                 verbose = TRUE) {
  alternative <- rlang::arg_match(alternative, values = AVAILABLE_ALTERNATIVES())

  # data preprocessing
  inputs <- twosamples2coeffs(data1, data2, mu, dx = dx)
  coeff1 <- inputs$coeff1
  coeff2 <- inputs$coeff2
  mu.eval <- inputs$mu

  n1 <- dim(coeff1)[1]
  n2 <- dim(coeff2)[1]
  p <- dim(coeff1)[2]
  n <- n1 + n2
  etichetta_ord <- c(rep(1, n1), rep(2, n2))
  coeff1 <- coeff1 - matrix(data = mu.eval, nrow = n1, ncol = p)

  #splines coefficients:
  eval <- coeff <- rbind(coeff1, coeff2)

  data.eval <- eval
  data.eval[1:n1, ] <- data.eval[1:n1, ] + matrix(
    data = mu.eval, nrow = n1, ncol = p
  )

  if (verbose)
    cli::cli_h1("Point-wise tests")
  
  #univariate permutations
  meandiff <- colMeans(coeff[1:n1, , drop = FALSE], na.rm = TRUE) - 
    colMeans(coeff[(n1 + 1):n, , drop = FALSE], na.rm = TRUE)
  sign.diff <- sign(meandiff)
  sign.diff[which(sign.diff == -1)] <- 0
  T0 <- switch(
    alternative,
    two.sided = (meandiff)^2,
    greater   = (meandiff * sign.diff)^2,
    less      = (meandiff * (sign.diff - 1))^2
  )

  T_coeff <- matrix(ncol = p, nrow = B)
  for (perm in 1:B) {
    if (paired) {
      if.perm <- stats::rbinom(n1, 1, 0.5)
      coeff_perm <- coeff
      for (couple in 1:n1) {
        if (if.perm[couple] == 1) {
          coeff_perm[c(couple, n1 + couple), ] <- coeff[c(n1 + couple, couple), ]
        }
      }
    } else {
      permutazioni <- sample(n)
      coeff_perm <- coeff[permutazioni, ]
    }
    meandiff <- colMeans(coeff_perm[1:n1, , drop = FALSE], na.rm = TRUE) - 
      colMeans(coeff_perm[(n1 + 1):n, , drop = FALSE], na.rm = TRUE)
    sign.diff <- sign(meandiff)
    sign.diff[which(sign.diff == -1)] <- 0
    T_coeff[perm, ] <- switch(
      alternative,
      two.sided = (meandiff)^2,
      greater   = (meandiff * sign.diff)^2,
      less      = (meandiff * (sign.diff - 1))^2
    )
  }
  
  pval <- numeric(p)
  for (i in 1:p) {
    pval[i] <- sum(T_coeff[, i] >= T0[i]) / B
  }

  #combination
  if (verbose)
    cli::cli_h1("Interval-wise tests")

  #asymmetric combination matrix:
  matrice_pval_asymm <- matrix(nrow = p, ncol = p)
  matrice_pval_asymm[p,] <- pval[1:p]
  T0_2x <- c(T0, T0)
  T_coeff_2x <- cbind(T_coeff, T_coeff)

  maxrow <- 1

  if (recycle) {
    for (i in (p - 1):maxrow) { # rows
      for (j in 1:p) { # columns
        inf <- j
        sup <- (p - i) + j
        T0_temp <- sum(T0_2x[inf:sup])
        T_temp <- rowSums(T_coeff_2x[, inf:sup])
        pval_temp <- sum(T_temp >= T0_temp) / B
        matrice_pval_asymm[i, j] <- pval_temp
      }
      
      if (verbose)
        cli::cli_h1("Creating the p-value matrix: end of row {p - i + 1} out of {p}")
    }
  } else { # without recycling
    for (i in (p - 1):maxrow) { # rows
      for (j in 1:i) { # columns
        inf <- j
        sup <- (p - i) + j
        T0_temp <- sum(T0_2x[inf:sup])
        T_temp <- rowSums(T_coeff_2x[, inf:sup])
        pval_temp <- sum(T_temp >= T0_temp) / B
        matrice_pval_asymm[i, j] <- pval_temp
      }
      
      if (verbose)
        cli::cli_h1("Creating the p-value matrix: end of row {p - i + 1} out of {p}")
    }
  }

  corrected.pval.matrix <- pval_correct(matrice_pval_asymm)
  corrected.pval <- corrected.pval.matrix[1, ]

  if (verbose)
    cli::cli_h1("Interval-Wise Testing completed")
  
  out <- list(
    test = '2pop',
    mu = mu.eval,
    adjusted_pval = corrected.pval,
    unadjusted_pval = pval,
    pval_matrix = matrice_pval_asymm,
    data.eval = data.eval,
    ord_labels = etichetta_ord
  )
  class(out) <- 'fdatest2'
  out
}
