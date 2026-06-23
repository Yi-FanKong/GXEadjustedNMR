#' Create nonlinear MR summary
#'
#' @param data Analysis data.frame.
#' @param outcome_col Outcome column.
#' @param exposure_col Exposure column.
#' @param prs_col PRS column.
#' @param covar_cols Covariate columns.
#' @param family "gaussian" or "binomial".
#' @param strata_method "ranked" or "residual".
#' @param q Number of strata.
#'
#' @return Output from SUMnlmr::create_nlmr_summary.
#' @export
create_nlmr_summary <- function(
    data,
    outcome_col = "outcome",
    exposure_col = "exposure_adjusted",
    prs_col = "PRS",
    covar_cols = NULL,
    family = c("gaussian", "binomial"),
    strata_method = c("ranked", "residual"),
    q = 10
) {
  family <- match.arg(family)
  strata_method <- match.arg(strata_method)

  data <- data.frame(data, check.names = FALSE)
  data <- stats::na.omit(data)

  y <- data[[outcome_col]]
  x <- data[[exposure_col]]
  g <- data[[prs_col]]

  covar <- NULL
  if (!is.null(covar_cols)) {
    covar <- as.matrix(data[, covar_cols, drop = FALSE])
  }

  SUMnlmr::create_nlmr_summary(
    y = y,
    x = x,
    g = g,
    covar = covar,
    family = family,
    strata_method = strata_method,
    q = q
  )
}
