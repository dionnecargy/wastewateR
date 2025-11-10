#' Check whether Dilutions are 1:10 Cq
#'
#' @param stats_df  Output from `summariseqPCR()`
#'
#' @returns A data frame with pass/fail
#' @export
#'
#' @import dplyr
#'
#' @author Dionne Argyropoulos
checkDilutions <- function(stats_df) {
  stats_df %>%
    dplyr::group_by(Fluor, Target, Sample) %>%
    dplyr::arrange(desc(LogCopy), .by_group = TRUE) %>%
    dplyr::mutate(
      mean_diff = abs(mean - dplyr::lead(mean)),        # absolute ΔCq
      log_diff  = LogCopy - dplyr::lead(LogCopy),
      exp_ratio = log_diff * 3.32,
      ok_dilution = dplyr::if_else(
        !is.na(mean_diff) & abs(mean_diff - exp_ratio) < 0.5,
        "PASS",
        dplyr::if_else(is.na(mean_diff), NA_character_, "CHECK")
      )
    ) %>%
    dplyr::ungroup()
}
