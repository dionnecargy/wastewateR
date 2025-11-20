#' Summarise qPCR Data with mean, min, max, sd, difference and range.
#'
#' @param qPCR_results  Output from `readqPCR()`.
#' @param Samples Does your data contain samples or just controls? Yes or No.
#'
#' @returns A data frame of summarised data.
#' @export
#'
#' @import dplyr
#'
#' @author Dionne Argyropoulos
summariseqPCR <- function(qPCR_results, Samples = "n"){

  Samples <- tolower(trimws(Samples)) # Normalize input: lowercase, trim spaces
  Samples <- substr(Samples, 1, 1) # Only use first character (y/n)

  df <- getReplicates(qPCR_results, Samples = Samples)

  stats <- df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      vals   = list(c_across(starts_with("Rep"))),  # collect repeats
      mean   = if (all(is.na(vals))) NA_real_ else mean(vals, na.rm = TRUE),
      median = if (all(is.na(vals))) NA_real_ else median(vals, na.rm = TRUE),
      min    = if (all(is.na(vals))) NA_real_ else min(vals, na.rm = TRUE),
      max    = if (all(is.na(vals))) NA_real_ else max(vals, na.rm = TRUE),
      sd     = if (all(is.na(vals))) NA_real_ else sd(vals, na.rm = TRUE),
      difference  = max - min,
      range_score = dplyr::case_when(
        all(is.na(vals))        ~ NA_character_,
        difference > 0.5        ~ "Replicate range > 0.5",
        NullFlag!="PASS"        ~ NA_character_,
        TRUE                    ~ "PASS"
      )
    ) %>%
    dplyr::select(-vals) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Fluor, Target, Plate) %>%
    dplyr::arrange(desc(LogCopy), .by_group = TRUE) %>%
    dplyr::mutate(
      # 1. Calculate the absolute change in Cq
      mean_diff = abs(mean - dplyr::lead(mean)),        # absolute ΔCq
      # 2. Calculate the diference in log copy number
      log_diff  = LogCopy - dplyr::lead(LogCopy),
      # 3. Calculate expected ratio
      exp_ratio = log_diff * 3.32,
      # 4. Determine if ratio is "PASS" or "FAIL"
      dilution = dplyr::if_else(
        !is.na(mean_diff) & abs(mean_diff - exp_ratio) < 0.5,
        "PASS",
        dplyr::if_else(is.na(mean_diff), NA_character_, "CHECK")
      )
    ) %>%
    dplyr::ungroup()

  return(stats)
}
