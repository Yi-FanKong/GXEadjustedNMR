#' Select genome-wide significant main-effect SNPs
#'
#' @param main_res fastGxE main-effect result file or data.frame.
#' @param p_cutoff P-value threshold.
#' @param snp_col SNP column name.
#' @param beta_col Beta column name.
#' @param p_col P-value column name.
#'
#' @return A data.frame of significant SNPs.
#' @export
select_main_snps <- function(
    main_res,
    p_cutoff = 5e-8,
    snp_col = "SNP",
    beta_col = "beta",
    p_col = "p"
) {
  dat <- if (is.character(main_res)) {
    data.table::fread(main_res)
  } else {
    data.table::as.data.table(main_res)
  }

  dat <- dat[!is.na(dat[[beta_col]]) & !is.na(dat[[p_col]])]
  dat <- dat[dat[[p_col]] < p_cutoff]

  data.frame(dat)
}


#' Clump significant GxE SNPs
#'
#' @param gxe_res fastGxE GxE result file or data.frame.
#' @param p_cutoff P-value threshold for GxE.
#' @param clump_kb Clumping window.
#' @param clump_r2 LD r2 threshold.
#' @param pop Population for TwoSampleMR clumping.
#' @param bfile Reference PLINK bfile for clumping.
#' @param plink_bin PLINK binary.
#'
#' @return A clumped data.frame.
#' @export
clump_gxe_snps <- function(
    gxe_res,
    p_cutoff = 5e-8,
    clump_kb = 10000,
    clump_r2 = 0.001,
    pop = "EUR",
    bfile = NULL,
    plink_bin = "plink"
) {
  dat <- if (is.character(gxe_res)) {
    data.table::fread(gxe_res)
  } else {
    data.table::as.data.table(gxe_res)
  }

  dat <- dat[!is.na(p_gxe)]
  dat <- dat[order(p_gxe)]

  dat_format <- dat |>
    dplyr::rename(
      pval.exposure = p_gxe,
      chr_name = chrom,
      chrom_start = base
    )

  clumped <- TwoSampleMR::clump_data(
    dat = dat_format,
    clump_kb = clump_kb,
    clump_r2 = clump_r2,
    clump_p1 = p_cutoff,
    clump_p2 = p_cutoff,
    pop = pop,
    bfile = bfile,
    plink_bin = plink_bin
  )

  data.frame(clumped)
}
