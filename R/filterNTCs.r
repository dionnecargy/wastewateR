#' Search and Filter No Template Controls (NTCs)
#'
#' This function specifically looks at the NTCs and confirm whether they are all "Na".
#' If they are "NA", then they are deleted from the data for further processing. If
#' they are not "NA", then an error message is displayed to check the NTCs.
#'
#' @param processed_qPCR Output list of dataframes from readqPCR
#'
#' @returns A cleaned results data frame with no NTCs
#' @export
#'
#' @import dplyr
#'
#' @author Dionne Argyropoulos
filterNTCs <- function(processed_qPCR) {

  df <- processed_qPCR[[3]]

  all_ntc_nan <- df %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Cq = as.numeric(ifelse(Cq %in% c("NaN", "", " "), NA, Cq))) %>%
    dplyr::filter(Content == "NTC" | Sample == "NTC") %>%
    dplyr::summarise(all_na = all(is.na(Cq))) %>%
    dplyr::pull(all_na)

  if (all_ntc_nan) {
    message("All NTCs are NaN")
  } else {
    stop("Check NTCs")
  }

  clean_df <- df %>% dplyr::filter(!(Content == "NTC" | Sample == "NTC"))

  return(clean_df)
}
