#' Extract SNP genotypes using PLINK and recode to additive dosage
#'
#' @param bfile_prefixes Character vector of PLINK bfile prefixes, e.g. chr1 to chr22.
#' @param snps Character vector of SNP IDs.
#' @param keep_file Optional PLINK keep file.
#' @param out_prefix Output prefix.
#' @param plink_bin PLINK executable.
#' @param verbose Print command.
#'
#' @return Path to PLINK raw file.
#' @export
extract_genotypes_plink <- function(
    bfile_prefixes,
    snps,
    keep_file = NULL,
    out_prefix,
    plink_bin = "plink",
    verbose = TRUE
) {
  out_dir <- dirname(out_prefix)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  snp_file <- paste0(out_prefix, ".snplist.txt")
  utils::write.table(
    snps,
    snp_file,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )

  merge_list <- paste0(out_prefix, ".merge.list.txt")
  if (file.exists(merge_list)) {
    file.remove(merge_list)
  }

  tmp_prefixes <- character()

  for (i in seq_along(bfile_prefixes)) {
    tmp_out <- paste0(out_prefix, ".tmp", i)

    args <- c(
      "--bfile", bfile_prefixes[i],
      "--extract", snp_file,
      "--make-bed",
      "--out", tmp_out
    )

    if (!is.null(keep_file)) {
      args <- c(args, "--keep", keep_file)
    }

    gxe_run_cmd(plink_bin, args, verbose = verbose)

    if (file.exists(paste0(tmp_out, ".bed"))) {
      tmp_prefixes <- c(tmp_prefixes, tmp_out)
    }
  }

  if (length(tmp_prefixes) == 0) {
    stop("No genotype files were created. Check SNP IDs and bfile prefixes.", call. = FALSE)
  }

  if (length(tmp_prefixes) == 1) {
    merged_prefix <- tmp_prefixes[1]
  } else {
    base <- tmp_prefixes[1]
    rest <- tmp_prefixes[-1]

    utils::write.table(
      rest,
      merge_list,
      quote = FALSE,
      sep = "\t",
      row.names = FALSE,
      col.names = FALSE
    )

    merged_prefix <- paste0(out_prefix, ".merged")

    gxe_run_cmd(
      plink_bin,
      c(
        "--bfile", base,
        "--merge-list", merge_list,
        "--make-bed",
        "--out", merged_prefix
      ),
      verbose = verbose
    )
  }

  recode_prefix <- paste0(out_prefix, ".code")

  gxe_run_cmd(
    plink_bin,
    c(
      "--bfile", merged_prefix,
      "--recode", "A",
      "--out", recode_prefix
    ),
    verbose = verbose
  )

  raw_file <- paste0(recode_prefix, ".raw")

  if (!file.exists(raw_file)) {
    stop("PLINK raw file not found: ", raw_file, call. = FALSE)
  }

  raw_file
}
