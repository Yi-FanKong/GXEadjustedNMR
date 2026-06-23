#' Run a system command safely
#'
#' @keywords internal
gxe_run_cmd <- function(cmd, args = character(), verbose = TRUE) {
  if (verbose) {
    message("[CMD] ", cmd, " ", paste(args, collapse = " "))
  }
  status <- system2(cmd, args = args)
  if (!identical(status, 0L)) {
    stop("Command failed: ", cmd, " ", paste(args, collapse = " "), call. = FALSE)
  }
  invisible(TRUE)
}

#' Strip .bed suffix from PLINK bfile prefix
#'
#' @keywords internal
gxe_strip_bed <- function(bfile) {
  sub("\\.bed$", "", bfile, ignore.case = TRUE)
}

#' Quote a string for Python code
#'
#' @keywords internal
gxe_py_quote <- function(x) {
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub('"', '\\"', x)
  paste0('"', x, '"')
}

#' Expand phenotype column range such as PC1:PC20 or Age:sleep
#'
#' @keywords internal
gxe_expand_col_range <- function(x, col_names) {
  x <- as.character(x)
  x <- unlist(strsplit(x, "\\s+"))
  x <- x[nzchar(x)]

  out <- unlist(lapply(x, function(z) {
    if (grepl(":", z, fixed = TRUE)) {
      parts <- strsplit(z, ":", fixed = TRUE)[[1]]
      if (length(parts) != 2) {
        stop("Invalid range specification: ", z, call. = FALSE)
      }

      start_col <- parts[1]
      end_col <- parts[2]

      i1 <- match(start_col, col_names)
      i2 <- match(end_col, col_names)

      if (is.na(i1)) {
        stop("Start column not found: ", start_col, call. = FALSE)
      }
      if (is.na(i2)) {
        stop("End column not found: ", end_col, call. = FALSE)
      }

      col_names[i1:i2]
    } else {
      if (!(z %in% col_names)) {
        stop("Column not found: ", z, call. = FALSE)
      }
      z
    }
  }), use.names = FALSE)

  unique(out)
}
