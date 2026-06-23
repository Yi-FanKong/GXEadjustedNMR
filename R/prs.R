#' Compute main-effect PRS from genotype dosage and SNP weights
#'
#' @param geno Genotype data.frame with ID and SNP dosage columns.
#' @param weights data.frame with SNP and beta columns.
#' @param id_col ID column.
#' @param snp_col SNP column in weights.
#' @param beta_col beta column in weights.
#'
#' @return data.frame with ID and PRS.
#' @export
compute_main_prs <- function(
    geno,
    weights,
    id_col = "ID",
    snp_col = "SNP",
    beta_col = "beta"
) {
  geno <- data.frame(geno, check.names = FALSE)
  weights <- data.frame(weights, check.names = FALSE)

  snps <- weights[[snp_col]]
  betas <- weights[[beta_col]]

  available <- snps %in% names(geno)

  if (!all(available)) {
    warning(
      "Some PRS SNPs are missing in genotype data: ",
      paste(snps[!available], collapse = ", ")
    )
  }

  snps <- snps[available]
  betas <- betas[available]

  if (length(snps) == 0) {
    stop("No PRS SNPs are available in genotype data.", call. = FALSE)
  }

  gmat <- as.matrix(geno[, snps, drop = FALSE])
  gmat[is.na(gmat)] <- 0

  prs <- as.numeric(gmat %*% betas)

  data.frame(
    ID = geno[[id_col]],
    PRS = prs
  )
}
