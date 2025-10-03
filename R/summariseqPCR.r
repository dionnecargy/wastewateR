#' Summarise qPCR Data with mean, min, max, sd, difference and range.
#'
#' @param df  Output from `organiseData()`
#'
#' @returns A data frame of summarised data.
#' @export
#'
#' @import dplyr
#' @import tidyr
#'
#' @author Dionne Argyropoulos
summariseqPCR <- function(df){
  stats <- df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      vals = list(c_across(starts_with("Rep"))),  # collect repeats
      mean = if (all(is.na(vals))) NA_real_ else mean(vals, na.rm = TRUE),
      min  = if (all(is.na(vals))) NA_real_ else min(vals, na.rm = TRUE),
      max  = if (all(is.na(vals))) NA_real_ else max(vals, na.rm = TRUE),
      sd   = if (all(is.na(vals))) NA_real_ else sd(vals, na.rm = TRUE),
      cv   = sd * 100,
      difference  = max - min,
      range_score = dplyr::case_when(
        all(is.na(vals))        ~ NA_character_,
        difference > 0.5        ~ "Replicate range > 0.5",
        NullFlag!="PASS"        ~ NA_character_,
        TRUE                    ~ "PASS"
      )
    ) %>%
    dplyr::select(-vals) %>%
    dplyr::ungroup()
  return(stats)
}
