#' Run fastGxE and optional mmSuSiE
#'
#' This function wraps the fastGxE workflow:
#' make GRM, process GRM, run GxE scan, optionally reformat GRM and run mmSuSiE.
#'
#' @param bfile PLINK bfile prefix or .bed file.
#' @param pheno Phenotype file path or data.frame.
#' @param trait Trait column name.
#' @param env_int Environmental variables used for GxE interaction. Supports ranges such as Age:sleep.
#' @param covar Covariates, e.g. "PC1:PC20".
#' @param class_covar Categorical covariates.
#' @param outdir Output directory.
#' @param prefix Output prefix.
#' @param fastgxe_path Path to fastGxE executable.
#' @param python_path Path to Python.
#' @param grm_prefix_existing Existing GRM prefix. If NULL, GRM will be built.
#' @param run_grm_reformat Whether to reformat GRM for mmSuSiE.
#' @param grm_npart Number of GRM partitions.
#' @param threads Number of threads.
#' @param cut_value Cut value for GRM grouping.
#' @param no_standardize_env Whether to disable environment standardization.
#' @param no_noisebye Whether to disable NoiseBye.
#' @param run_mmsusie Whether to run mmSuSiE.
#' @param mmsusie_snp User-specified SNP for mmSuSiE.
#' @param mmsusie_p_cutoff P-value cutoff for mmSuSiE lead SNP extraction.
#' @param mmsusie_ld_r2 LD r2 threshold.
#' @param verbose Print commands.
#'
#' @return A list of output file paths.
#' @export
run_fastgxe <- function(
    bfile,
    pheno,
    trait,
    env_int,
    covar = NULL,
    class_covar = NULL,
    outdir = "output",
    prefix = "test",
    fastgxe_path = "fastgxe",
    python_path = "python",
    grm_prefix_existing = NULL,
    run_grm_reformat = TRUE,
    grm_npart = NULL,
    threads = NULL,
    cut_value = 0.05,
    no_standardize_env = FALSE,
    no_noisebye = FALSE,
    run_mmsusie = TRUE,
    mmsusie_snp = NULL,
    mmsusie_p_cutoff = 5e-8,
    mmsusie_ld_r2 = 0.1,
    verbose = TRUE
) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  bfile <- gxe_strip_bed(bfile)

  if (is.data.frame(pheno)) {
    pheno_file <- file.path(outdir, paste0(prefix, "_pheno.txt"))
    utils::write.table(
      pheno,
      pheno_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE,
      na = "NA"
    )
    pheno_names <- names(pheno)
  } else {
    pheno_file <- pheno
    pheno_names <- names(utils::read.table(
      pheno_file,
      header = TRUE,
      nrows = 0,
      check.names = FALSE
    ))
  }

  env_int_mmsusie <- gxe_expand_col_range(env_int, pheno_names)

  if (is.null(grm_prefix_existing)) {
    grm_prefix <- file.path(outdir, basename(bfile))
  } else {
    grm_prefix <- grm_prefix_existing
  }

  if (!is.null(grm_npart)) {
    grm_npart <- as.integer(grm_npart)
    if (length(grm_npart) != 1 || is.na(grm_npart) || grm_npart < 2) {
      stop("grm_npart must be NULL or an integer >= 2.", call. = FALSE)
    }
  }

  fastgxe_threads_args <- if (is.null(threads)) {
    character(0)
  } else {
    c("--threads", as.character(threads))
  }

  gxe_prefix <- file.path(outdir, paste0(prefix, "_gxe"))
  mmsusie_prefix <- file.path(outdir, paste0(prefix, "_mmsusie"))

  if (is.null(grm_prefix_existing)) {
    if (is.null(grm_npart)) {
      gxe_run_cmd(
        fastgxe_path,
        c(
          "--make-grm",
          "--bfile", bfile,
          "--out", grm_prefix,
          fastgxe_threads_args
        ),
        verbose = verbose
      )
    } else {
      for (part in seq_len(grm_npart)) {
        gxe_run_cmd(
          fastgxe_path,
          c(
            "--make-grm",
            "--bfile", bfile,
            "--npart", as.character(grm_npart), as.character(part),
            "--out", grm_prefix,
            fastgxe_threads_args
          ),
          verbose = verbose
        )
      }

      gxe_run_cmd(
        fastgxe_path,
        c(
          "--process-grm",
          "--merge",
          "--grm", grm_prefix,
          "--npart", as.character(grm_npart),
          fastgxe_threads_args
        ),
        verbose = verbose
      )
    }

    gxe_run_cmd(
      fastgxe_path,
      c(
        "--process-grm",
        "--group",
        "--grm", grm_prefix,
        "--cut-value", as.character(cut_value),
        fastgxe_threads_args
      ),
      verbose = verbose
    )
  } else {
    if (verbose) {
      message("[INFO] Use existing processed GRM: ", grm_prefix)
    }
  }

  gxe_args <- c(
    "--test-gxe",
    "--grm", grm_prefix,
    "--bfile", bfile,
    "--data", pheno_file,
    "--trait", trait,
    fastgxe_threads_args
  )

  if (!is.null(covar)) {
    gxe_args <- c(gxe_args, "--covar", covar)
  }

  if (!is.null(class_covar)) {
    gxe_args <- c(gxe_args, "--class", class_covar)
  }

  gxe_args <- c(gxe_args, "--env-int", env_int)

  if (no_standardize_env) {
    gxe_args <- c(gxe_args, "--no-standardize-env")
  }

  if (no_noisebye) {
    gxe_args <- c(gxe_args, "--no-noisebye")
  }

  gxe_args <- c(gxe_args, "--out", gxe_prefix)

  gxe_run_cmd(fastgxe_path, gxe_args, verbose = verbose)

  if (run_mmsusie) {
    if (run_grm_reformat) {
      gxe_run_cmd(
        fastgxe_path,
        c(
          "--process-grm",
          "--reformat",
          "--grm", grm_prefix,
          "--out-fmt", "1",
          "--out", grm_prefix,
          fastgxe_threads_args
        ),
        verbose = verbose
      )
    }

    env_py <- paste0(
      "[",
      paste(vapply(env_int_mmsusie, gxe_py_quote, character(1)), collapse = ", "),
      "]"
    )

    snp_py <- if (is.null(mmsusie_snp)) {
      "None"
    } else {
      gxe_py_quote(mmsusie_snp)
    }

    py_file <- tempfile(fileext = ".py")

    py_code <- paste0(
      "from mmsusie import MMSuSiESp\n",
      "model = MMSuSiESp()\n\n",

      "assoc_file = ", gxe_py_quote(paste0(gxe_prefix, ".res")), "\n",
      "bed_file = ", gxe_py_quote(bfile), "\n",
      "grm_file = ", gxe_py_quote(grm_prefix), "\n",
      "pheno_file = ", gxe_py_quote(pheno_file), "\n",
      "trait = ", gxe_py_quote(trait), "\n",
      "env_int = ", env_py, "\n",
      "varcom_file = ", gxe_py_quote(paste0(gxe_prefix, ".var")), "\n",
      "out_file = ", gxe_py_quote(mmsusie_prefix), "\n",
      "snp_id = ", snp_py, "\n\n",

      "df_lead = model.ld_pure(\n",
      "    assoc_file,\n",
      "    bed_file,\n",
      "    ld_r2=", mmsusie_ld_r2, ",\n",
      "    snp='SNP',\n",
      "    p='p_gxe',\n",
      "    p_cutoff=", mmsusie_p_cutoff, "\n",
      ")\n\n",

      "df_lead.to_csv(out_file + '.leading_snps.txt', sep='\\t', index=False)\n\n",

      "if snp_id is None:\n",
      "    if df_lead.shape[0] == 0:\n",
      "        raise SystemExit('No lead SNP found. Try larger p_cutoff or specify mmsusie_snp.')\n",
      "    df_lead = df_lead.sort_values('p_gxe')\n",
      "    snp_id = str(df_lead.iloc[0]['SNP'])\n\n",

      "with open(out_file + '.selected_snp.txt', 'w') as f:\n",
      "    f.write(snp_id + '\\n')\n\n",

      "res = model.mmsusie_lead_gxe(\n",
      "    pheno_file,\n",
      "    trait,\n",
      "    env_int,\n",
      "    grm_file,\n",
      "    bed_file,\n",
      "    snp_id,\n",
      "    varcom_file,\n",
      "    out_file,\n",
      "    L=10,\n",
      "    maxiter=100,\n",
      "    tol=1e-3,\n",
      "    coverage=0.95,\n",
      "    min_abs_corr=0.5,\n",
      "    estimate_sigma=False\n",
      ")\n\n",

      "print('mmSuSiE finished')\n",
      "print('Selected SNP:', snp_id)\n"
    )

    writeLines(py_code, py_file)
    gxe_run_cmd(python_path, py_file, verbose = verbose)
  }

  list(
    pheno_file = pheno_file,

    grm_bin = paste0(grm_prefix, ".grm.bin"),
    grm_id = paste0(grm_prefix, ".grm.id"),
    grm_N = paste0(grm_prefix, ".grm.N"),
    grm_group = paste0(grm_prefix, ".grm.group"),
    grm_group_size = paste0(grm_prefix, ".grm.group.size"),
    grm_index_triplet = paste0(grm_prefix, ".grm.index_triplet"),

    gxe_res = paste0(gxe_prefix, ".res"),
    gxe_var = paste0(gxe_prefix, ".var"),
    gxe_main_random_res = paste0(gxe_prefix, ".main.random.res"),
    gxe_no_main_random_res = paste0(gxe_prefix, ".GxEnoMain.random.res"),
    gxe_random_res = paste0(gxe_prefix, ".GxE.random.res"),

    mmsusie_leading_snps = paste0(mmsusie_prefix, ".leading_snps.txt"),
    mmsusie_selected_snp = paste0(mmsusie_prefix, ".selected_snp.txt"),
    mmsusie_pip = paste0(mmsusie_prefix, ".pip.txt"),
    mmsusie_alpha = paste0(mmsusie_prefix, ".alpha.txt"),
    mmsusie_mu = paste0(mmsusie_prefix, ".mu.txt"),
    mmsusie_cs = paste0(mmsusie_prefix, ".cs.txt")
  )
}
