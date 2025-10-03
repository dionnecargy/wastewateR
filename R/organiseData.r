#' Title
#'
#' @param filtered_data Output from `organiseData()`
#' @param target Does the "Target" column have data in it? Either supply another
#' data frame with the Sample and Target or leave blank. Default is NULL.
#'
#' @returns A data frame of organised and cleaned data
#' @export
#'
#' @import dplyr
#' @import tidyr
#'
#' @author Dionne Argyropoulos
organiseData <- function(filtered_data, target = NULL){

  df <- filtered_data %>%
    dplyr::mutate(
      CopyRxn = case_when(
        is.na(Biological_Set_Name) ~ NA_real_,
        TRUE ~ 10 ^ as.numeric(str_remove(Biological_Set_Name, "10\\^"))
      ),
      LogCopy = log10(CopyRxn),
      Target = na_if(trimws(Target), "")
    )

  # Check if "Target" column is empty in dataset and needs to be supplied
  target_col <- any(is.na(df$Target))

  if (target_col == TRUE){
    df <- df %>%
      dplyr::select(-Target) %>%
      dplyr::left_join(target, by = "Sample")
  } else {
    df <- df
  }

  df <- df %>%
    dplyr::select(Fluor, Target, LogCopy, CopyNum = Biological_Set_Name, CopyRxn, everything()) %>%
    dplyr::group_by(Fluor, Target, LogCopy, Sample) %>%
    dplyr::mutate(
      Cq = as.numeric(Cq),
      Rep = paste0("Rep", row_number())
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-Well) %>%
    tidyr::pivot_wider(
      id_cols = Fluor:Sample,
      names_from = Rep,
      values_from = Cq
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      NullFlag = {
        repeat_values <- c_across(starts_with("Rep"))
        repeat_names <- names(across(starts_with("Rep")))
        first_na <- which(is.na(repeat_values))[1]
        if (!is.na(first_na)) paste0(repeat_names[first_na], " is NA") else "PASS" # Flag which are NULL
      }
    ) %>%
    dplyr::ungroup()

  return(df)

}
