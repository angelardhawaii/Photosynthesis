# Code to clean RLC files created in Walz WinControl - to replace Python script
# By Angela Richards Donà
# August 20, 2025

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(tibble)
})

output_columns <- c(
  "date","time","id","f","f0","fm","fm_prime","epar","yield","etr","fv_fm","npq","delta_npq","retr"
)

# Constants to match your Python
etr_factor  <- 0.84
psII_factor <- 0.5

# ---- helpers ---------------------------------------------------------------

detect_sep <- function(lines) {
  # detect the first header-ish line and choose sep based on presence of ';'
  header_line <- lines[grepl('(^|")Date("|)\\b|(^|")Datetime("|)\\b', lines)][1]
  if (is.na(header_line)) {
    stop("At least one line in the raw PAM data file must contain headers like Datetime or Date.")
  }
  if (str_detect(header_line, ";")) ";" else ","
}

split_line <- function(line, sep) {
  # simple splitter that respects the chosen separator
  # (files from WinControl are simple enough that this works fine)
  str_split(line, fixed(sep), n = Inf)[[1]] |> str_replace_all('^"|"$', "") |> trimws()
}

build_col_positions <- function(lines, sep) {
  hdr_line <- lines[grepl(paste0("(?:^|[\"])(Date|Datetime)(?:[\"]|)", sep), lines)][1]
  if (is.na(hdr_line)) {
    stop("Could not find header line with Date/Datetime.")
  }
  hdr <- split_line(hdr_line, sep)
  setNames(seq_along(hdr), hdr)
}

get_value <- function(fields, col_pos, name) {
  if (!name %in% names(col_pos)) return("NA")
  p <- col_pos[[name]]
  if (p > length(fields)) stop(sprintf('Cannot find data for column %s in row %s', name, paste(fields, collapse="|")))
  val <- fields[p]
  val <- str_replace_all(val, '^"|"$', "")
  val
}

record_type_of <- function(fields, col_pos) {
  nm <- "Type"
  if (length(fields) <= 4) return("N/A")
  get_value(fields, col_pos, nm)
}

calc_r_etr <- function(par, yii) {
  if (yii == "-") return("-")
  round(as.numeric(par) * as.numeric(yii), 2)
}

# OutputRow equivalent
make_output_row <- function(fields, col_pos) {
  has_date_time <- all(c("Date", "Time") %in% names(col_pos))
  has_datetime  <- "Datetime" %in% names(col_pos)
  
  if (has_date_time) {
    date <- get_value(fields, col_pos, "Date")
    time <- get_value(fields, col_pos, "Time")
  } else if (has_datetime) {
    dt <- get_value(fields, col_pos, "Datetime")
    parts <- str_split(dt, "\\s+", n = 2)[[1]]
    date <- parts[1]
    time <- ifelse(length(parts) > 1, parts[2], "")
  } else {
    stop("Either DATE or DATETIME must be in the source data file.")
  }
  
  f        <- get_value(fields, col_pos, "1:F")
  fm_prime <- get_value(fields, col_pos, "1:Fm'")
  par      <- get_value(fields, col_pos, "1:PAR")
  yii      <- get_value(fields, col_pos, "1:Y (II)")
  etr      <- get_value(fields, col_pos, "1:ETR")
  npq      <- get_value(fields, col_pos, "1:NPQ")
  f0       <- get_value(fields, col_pos, "1:Fo'")
  fm       <- get_value(fields, col_pos, "1:Fm")
  fvfm_raw <- get_value(fields, col_pos, "1:Fv/Fm")
  
  # backfill Y(II) if '-' and Fm' != 0
  if (yii == "-" && !is.na(suppressWarnings(as.numeric(fm_prime))) && as.numeric(fm_prime) != 0) {
    yii <- round((as.numeric(fm_prime) - as.numeric(f)) / as.numeric(fm_prime), 3)
  }
  
  # backfill ETR if yii available and ETR was '-'
  if (!identical(yii, "-") && etr == "-") {
    etr <- round(as.numeric(par) * etr_factor * psII_factor * as.numeric(yii), 2)
  }
  
  r_etr <- if (identical(yii, "-")) "-" else round(as.numeric(par) * as.numeric(yii), 2)
  
  list(
    date = date, time = time, f = f, f0 = f0, fm = fm, fm_prime = fm_prime, par = par,
    yii = yii, etr = etr, fvfm_raw = fvfm_raw, npq = npq, r_etr = r_etr
  )
}

# ---- main per-file processor -----------------------------------------------

process_lines <- function(lines, warn_prefix = NULL) {
  sep  <- detect_sep(lines)
  cols <- build_col_positions(lines, sep)
  
  # Light-curve state
  lc_sample_id <- NULL
  lc_f         <- NA_character_
  lc_npq_zero  <- -1
  lc_delta_npq <- 0
  buffer       <- list()  # list of output rows (as character vectors)
  out_rows     <- list()  # list of buffers to be appended when SLCE ends
  
  flush_buffer <- function() {
    if (length(buffer) == 0) return(invisible(NULL))
    if (length(buffer) != 9) {
      # Mirror Python's warning behavior
      msg <- sprintf("Warning: less than 9 records (%d) found in light curve starting %s",
                     length(buffer),
                     if (length(buffer) > 0) paste(buffer[[1]][1], buffer[[1]][2]) else "N/A")
      if (!is.null(warn_prefix)) msg <- paste(warn_prefix, msg)
      message(msg)
    } else {
      out_rows <<- c(out_rows, buffer)
    }
    # reset per-curve
    lc_sample_id <<- NULL
    lc_f         <<- NA_character_
    lc_npq_zero  <<- -1
    lc_delta_npq <<- 0
    buffer       <<- list()
    invisible(NULL)
  }
  
  for (ln in lines) {
    fields <- split_line(ln, sep)
    # skip empty lines
    if (length(fields) == 0) next
    
    rtype <- record_type_of(fields, cols)
    
    if (rtype == "SLCS") {
      # start a new light curve, any prior one is discarded if not terminated
      lc_sample_id <- fields[length(fields)] |> str_replace_all('^"|"$', "") |> trimws()
      lc_f         <- NA_character_
      lc_npq_zero  <- -1
      lc_delta_npq <- 0
      buffer       <- list()
      
    } else if (rtype == "SLCE") {
      # end of curve -> write buffer if size==9 else warn
      flush_buffer()
      
    } else if (rtype == "FO") {
      o <- make_output_row(fields, cols)
      lc_f <- o$f
      # FO row: F and F0 are both 'o$f' per your Python
      row <- c(
        o$date, o$time, lc_sample_id %||% "NA",
        o$f, o$f, o$fm, o$fm_prime, o$par,
        as.character(o$yii), as.character(o$etr), o$fvfm_raw,
        o$npq, "0.0", as.character(o$r_etr)
      )
      names(row) <- output_columns
      buffer <- append(buffer, list(row))
      
    } else if (rtype == "F") {
      o <- make_output_row(fields, cols)
      
      if (lc_npq_zero == -1) {
        lc_npq_zero <- ifelse(o$npq == "-", -1, suppressWarnings(as.numeric(str_replace_all(o$npq, " ", ""))))
      }
      
      # mirror Python's condition: set delta_npq when current buffer size == 10 (before appending)
      if (length(buffer) == 10) {
        lc_delta_npq <- if (o$npq == "-") "-" else {
          cur <- suppressWarnings(as.numeric(str_replace_all(o$npq, " ", "")))
          if (is.na(cur) || lc_npq_zero == -1) "-" else round(cur - lc_npq_zero, 3)
        }
      }
      
      row <- c(
        o$date, o$time, lc_sample_id %||% "NA",
        o$f, lc_f %||% "NA", o$fm, o$fm_prime, o$par,
        as.character(o$yii), as.character(o$etr), o$fvfm_raw,
        o$npq, as.character(lc_delta_npq), as.character(o$r_etr)
      )
      names(row) <- output_columns
      buffer <- append(buffer, list(row))
    }
    # all other record types ignored (match Python behavior)
  }
  
  # if a file ends without SLCE, we don't flush in Python; we'll keep same behavior
  # Return combined tibble for rows captured
  if (length(out_rows) == 0) {
    tibble::as_tibble(setNames(rep(list(character(0)), length(output_columns)), output_columns))
  } else {
    as_tibble(do.call(rbind, out_rows), .name_repair = "minimal") |>
      setNames(output_columns)
  }
}

# ---- top-level function -----------------------------------------------------

#' Clean WinControl CSV files (Walz) and write a single tidy CSV
#' @param data_glob A glob pattern, e.g. "/path/to/*.csv"
#' @param output_file Output CSV path
pam_data_cleaner_tidy <- function(data_glob, output_file) {
  files <- Sys.glob(data_glob)
  if (length(files) == 0) stop("No files matched: ", data_glob)
  
  message("Found ", length(files), " files…")
  all_rows <-
    files |>
    map(function(f) {
      message("Processing: ", f)
      lines <- readr::read_lines(f, progress = FALSE)
      process_lines(lines, warn_prefix = paste0("[", basename(f), "] "))
    }) |>
    list_rbind()
  
  # enforce column order and write
  all_rows <- all_rows |>
    mutate(across(everything(), as.character)) |>
    select(all_of(output_columns))
  
  # create dir and write
  dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(all_rows, output_file, na = "")
  
  message("Wrote: ", output_file, " (", nrow(all_rows), " rows)")
  invisible(all_rows)
}

# small helper to coalesce NULL to value
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- Example usage (uncomment and edit) -------------------------------------
# result <- pam_data_cleaner_tidy("/path/to/wincontrol/*.csv", "cleaned/walz_clean.csv")
# View(result)

#copy this next line into R script to run the function
#source("/Users/angela/src/Photosynthesis/transform/scripts/pam_data_cleaner_tidy.R")

rlc_2025 <- pam_data_cleaner_tidy(
  data_glob = "/Users/angela/src/Photosynthesis/data/wastewater/input/pam_2025/*.csv",
  output_file = "/Users/angela/src/Photosynthesis/data/wastewater/transformed/pam_waste_p2_clean_tidy.csv")
glimpse(rlc_2025)
rlc_2023 <- pam_data_cleaner_tidy(
  data_glob = "/Users/angela/src/Photosynthesis/data/wastewater/input/pam_2023/*.csv",
  output_file = "/Users/angela/src/Photosynthesis/data/wastewater/transformed/pam_2023_clean_tidy.csv"
)
