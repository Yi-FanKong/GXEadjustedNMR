#' Compute GxE adjustment score
#'
#' Computes:
#' \deqn{
#'   A_i = \sum_{(j,k) \in S} G_{ij} E_{ik} \hat\beta^{GxE}_{jk}
#' }
#'
#' @param geno Genotype data.frame. Must contain ID column and SNP dosage columns.
#' @param env Environmental data.frame. Must contain ID column and environmental variables.
#' @param beta_wide data.frame with one row per SNP and one column per environmental factor.
#'        The first column should contain SNP IDs.
#' @param id_col ID column name.
#' @param snp_col SNP column name in beta_wide.
#'
#' @return data.frame with ID and GxE adjustment score.
#' @export
compute_gxe_adjustment <- function(
    geno,
    env,
    beta_wide,
    id_col = "ID",
    snp_col = "SNP"
) {
  geno <- data.frame(geno, check.names = FALSE)
  env <- data.frame(env, check.names = FALSE)
  beta_wide <- data.frame(beta_wide, check.names = FALSE)

  if (!(id_col %in% names(geno))) {
    stop("id_col not found in geno: ", id_col, call. = FALSE)
  }

  if (!(id_col %in% names(env))) {
    stop("id_col not found in env: ", id_col, call. = FALSE)
  }

  if (!(snp_col %in% names(beta_wide))) {
    stop("snp_col not found in beta_wide: ", snp_col, call. = FALSE)
  }

  common_id <- intersect(geno[[id_col]], env[[id_col]])

  geno <- geno[match(common_id, geno[[id_col]]), , drop = FALSE]
  env <- env[match(common_id, env[[id_col]]), , drop = FALSE]

  beta_env_cols <- setdiff(names(beta_wide), snp_col)

  missing_env <- setdiff(beta_env_cols, names(env))
  if (length(missing_env) > 0) {
    stop("Environmental columns not found in env: ",
         paste(missing_env, collapse = ", "), call. = FALSE)
  }

  adj <- rep(0, length(common_id))

  for (e in beta_env_cols) {
    idx <- which(!is.na(beta_wide[[e]]) & beta_wide[[e]] != 0)

    if (length(idx) == 0) {
      next
    }

    snps <- beta_wide[[snp_col]][idx]
    betas <- beta_wide[[e]][idx]

    available <- snps %in% names(geno)

    if (!all(available)) {
      warning(
        "Some SNPs for environment ", e, " are missing in genotype data: ",
        paste(snps[!available], collapse = ", ")
      )
    }

    snps <- snps[available]
    betas <- betas[available]

    if (length(snps) == 0) {
      next
    }

    gmat <- as.matrix(geno[, snps, drop = FALSE])
    gmat[is.na(gmat)] <- 0

    env_vec <- env[[e]]
    env_vec[is.na(env_vec)] <- 0

    adj <- adj + as.numeric(gmat %*% betas) * env_vec
  }

  data.frame(
    ID = common_id,
    GxE_adjustment = adj
  )
}


#' Make GxE-adjusted exposure
#'
#' Computes:
#' \deqn{
#'   X' = X - A
#' }
#'
#' @param exposure data.frame with ID and exposure.
#' @param adjustment data.frame from gxe_compute_gxe_adjustment.
#' @param id_col ID column.
#' @param exposure_col Exposure column.
#' @param adjustment_col Adjustment column.
#'
#' @return data.frame with ID, exposure, GxE_adjustment, exposure_adjusted.
#' @export
make_adjusted_exposure <- function(
    exposure,
    adjustment,
    id_col = "ID",
    exposure_col = "exposure",
    adjustment_col = "GxE_adjustment"
) {
  exposure <- data.frame(exposure, check.names = FALSE)
  adjustment <- data.frame(adjustment, check.names = FALSE)

  dat <- dplyr::inner_join(
    exposure,
    adjustment,
    by = id_col
  )

  dat$exposure_adjusted <- dat[[exposure_col]] - dat[[adjustment_col]]

  dat
}
