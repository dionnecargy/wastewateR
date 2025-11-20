#' Summarise qPCR Data with mean, min, max, sd, difference and range.
#'
#' @param qPCR_results  Output from `readqPCR()`.
#' @param Samples Does your data contain samples or just controls? Yes or No.
#'
#' @returns A data frame of summarised data.
#' @export
#'
#' @import dplyr
#' @import tidyr
#'
#' @author Dionne Argyropoulos
getReplicates <- function(qPCR_results, Samples = "n"){

  # Load qPCR results
  df <- qPCR_results[[1]]

  # Are you also considering unique isolate samples?
  Samples <- tolower(trimws(Samples)) # Normalize input: lowercase, trim spaces
  Samples <- substr(Samples, 1, 1) # Only use first character (y/n)

  # Identify replicates of each target per fluor, per logcopy, per plex, per plate
  if(Samples == "y"){

    stds <- df %>%
      dplyr::filter(str_detect(Sample, "STD")) %>%
      dplyr::group_by(Fluor, LogCopy, Target, Plex, Plate) %>%
      dplyr::mutate(
        Cq = as.numeric(Cq),
        Rep = paste0("Rep", row_number())
      ) %>%
      tidyr::pivot_wider(
        id_cols = Fluor:Plate,
        names_from = Rep,
        values_from = Cq
      ) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        NullFlag = {
          repeat_values <- c_across(starts_with("Rep"))
          repeat_names <- names(across(starts_with("Rep")))

          # Count non-NA values
          n_non_na <- sum(!is.na(repeat_values))

          if (n_non_na >= 3) {
            "PASS"
          } else {
            # Flag all NAs, column names joined by ", " and a single " is NA"
            na_indices <- which(is.na(repeat_values))
            if (length(na_indices) > 0) {
              paste0(paste(repeat_names[na_indices], collapse = ", "), " is NA")
            } else {
              "PASS"
            }
          }
        }
      ) %>%
      dplyr::ungroup()

    samples <-  df %>%
      filter(!str_detect(Sample, "STD")) %>%
      dplyr::group_by(Fluor, Sample, Target, Plex, Plate) %>%
      drop_na(Cq) %>%
      dplyr::mutate(
        Cq = as.numeric(Cq),
        Rep = paste0("Rep", row_number())
      ) %>%
      tidyr::pivot_wider(
        id_cols = Fluor:Plate,
        names_from = Rep,
        values_from = Cq
      ) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        NullFlag = {
          repeat_values <- c_across(starts_with("Rep"))
          repeat_names <- names(across(starts_with("Rep")))

          # Count non-NA values
          n_non_na <- sum(!is.na(repeat_values))

          if (n_non_na >= 3) {
            "PASS"
          } else {
            # Flag all NAs, column names joined by ", " and a single " is NA"
            na_indices <- which(is.na(repeat_values))
            if (length(na_indices) > 0) {
              paste0(paste(repeat_names[na_indices], collapse = ", "), " is NA")
            } else {
              "PASS"
            }
          }
        }
      ) %>%
      dplyr::ungroup()

    replicates <- dplyr::bind_rows(stds, samples) %>%
      dplyr::select(Fluor:Plate, starts_with("Rep"), NullFlag)

  } else {

    replicates <- df %>%
      dplyr::group_by(Fluor, LogCopy, Target, Plex, Plate) %>%
      dplyr::mutate(
        Cq = as.numeric(Cq),
        Rep = paste0("Rep", row_number())
      ) %>%
      tidyr::pivot_wider(
        id_cols = Fluor:Plate,
        names_from = Rep,
        values_from = Cq
      ) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        NullFlag = {
          repeat_values <- c_across(starts_with("Rep"))
          repeat_names <- names(across(starts_with("Rep")))

          # Count non-NA values
          n_non_na <- sum(!is.na(repeat_values))

          if (n_non_na >= 3) {
            "PASS"
          } else {
            # Flag all NAs, column names joined by ", " and a single " is NA"
            na_indices <- which(is.na(repeat_values))
            if (length(na_indices) > 0) {
              paste0(paste(repeat_names[na_indices], collapse = ", "), " is NA")
            } else {
              "PASS"
            }
          }
        }
      ) %>%
      dplyr::ungroup()

  }

  return(replicates)
}
