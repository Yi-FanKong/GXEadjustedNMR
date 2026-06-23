#' Plot LACE estimates across strata
#'
#' @param lace_data data.frame with beta, lci, uci, strata, Method.
#' @param exposure Exposure name.
#' @param outcome Outcome name.
#' @param family "gaussian" or "binomial".
#'
#' @return ggplot object.
#' @export
plot_lace <- function(
    lace_data,
    exposure = "Exposure",
    outcome = "Outcome",
    family = c("gaussian", "binomial")
) {
  family <- match.arg(family)

  plot_data <- data.frame(lace_data)

  if (family == "binomial") {
    ylab <- paste0("OR for ", outcome)
    hline <- 1
  } else {
    ylab <- paste0("Beta for ", outcome)
    hline <- 0
  }

  ggplot2::ggplot(plot_data, ggplot2::aes(x = strata, y = beta, color = Method)) +
    ggplot2::geom_point(
      position = ggplot2::position_dodge(width = 0.5),
      size = 3
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lci, ymax = uci),
      width = 0,
      position = ggplot2::position_dodge(width = 0.5),
      linewidth = 1
    ) +
    ggplot2::geom_hline(yintercept = hline, linetype = "solid") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = paste0("Strata of ", exposure),
      y = ylab,
      color = "Method"
    )
}


#' Plot genetic association heterogeneity
#'
#' @param heter_data data.frame with beta, lci, uci, strata, Method.
#' @param exposure Exposure name.
#'
#' @return ggplot object.
#' @export
plot_heterogeneity <- function(
    heter_data,
    exposure = "Exposure"
) {
  plot_data <- data.frame(heter_data)

  ggplot2::ggplot(plot_data, ggplot2::aes(x = strata, y = beta, color = Method)) +
    ggplot2::geom_point(
      position = ggplot2::position_dodge(width = 0.5),
      size = 3
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lci, ymax = uci),
      width = 0,
      position = ggplot2::position_dodge(width = 0.5),
      linewidth = 1
    ) +
    # ggplot2::theme_minimal() +
    ggplot2::labs(
      x = paste0("Strata of ", exposure),
      y = paste0("Genetic association for ", exposure),
      color = "Method"
    )
}
