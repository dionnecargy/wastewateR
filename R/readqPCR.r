#' Read raw qPCR Data
#'
#' A function to read in raw qPCR data
#'
#' @param raw_data One file path (.csv)
#' @param MPlex Does your data contain multiplex data? "Yes" or "No". Default = "No.
#' @param Samples Does your data contain samples or just controls? "Yes" or "No".
#'
#' @returns A list of data frames containing (i) raw data (untouched),
#' (ii) run information on qPCR machine, (iii) results data frame
#' @export
#'
#' @import dplyr
#' @import tidyr
#' @importFrom utils read.csv
#'
#' @author Dionne Argyropoulos
readqPCR <- function(raw_data, MPlex = "No", Samples){

  # User inputs for MPlex and Samples
  MPlex <- tolower(trimws(MPlex)) # Normalize input: lowercase, trim spaces
  MPlex <- substr(MPlex, 1, 1) # Only use first character (y/n)
  Samples <- tolower(trimws(Samples)) # Normalize input: lowercase, trim spaces
  Samples <- substr(Samples, 1, 1) # Only use first character (y/n)

  # Data processing
  final_df <- list()
  df <- read.csv(here::here(raw_data), na.strings = c("NaN"))

  # Store row numbers
  run_info_end_row <- which(df$File.Name == "Well")
  final_row <- nrow(df)

  # Save Raw df information
  final_df$raw_data <- df

  # Save Run information
  final_df$run_info <- df[1:(run_info_end_row-1),]

  # Re-arrange Biological Set Name and Sample
  results <- df[run_info_end_row:final_row, ]

  # Save Results
  results <- results %>% janitor::row_to_names(row_number = 1)

  # test if column looks like numbers of form 10^N
  is_numeric_like <- any(str_detect(na.omit(results$`Biological Set Name`), "^10\\^[0-9]+$"))

  if (is_numeric_like == FALSE) {
    results <- results %>% dplyr::rename(Biological_Set_Name = Sample, Sample = `Biological Set Name`)
  } else {
    results <- results %>% dplyr::rename(Biological_Set_Name = `Biological Set Name`)
  }

  # Does your data contain multiplex or simplex ?
  if (MPlex == "y"){
    results <- results %>%
      tidyr::drop_na(Cq) %>%
      dplyr::group_by(Fluor, Target) %>%
      dplyr::mutate(
        # find the "alternate" Biological Set Name (non-MPlex) for the same Fluor/Target
        alt_name = Sample[Sample != "MPlex"][1],
        # if row is MPlex, append alt_name; else leave unchanged
        Sample = ifelse(Sample == "MPlex", paste0(alt_name, "_MPlex"), Sample)
      ) %>%
      dplyr::select(-alt_name )
  } else {
    results <- results
  }


  # Does your data contain samples and standards?
  if (Samples == "y"){
    results <- results %>%
      dplyr::mutate(
        Sample = case_when(
          str_detect(Biological_Set_Name, "^10\\^") ~ "Control",
          .default = Biological_Set_Name
        ),
        Biological_Set_Name = case_when(
          str_detect(Biological_Set_Name, "^10\\^") ~ Biological_Set_Name,
          .default = NA
        )
      ) %>%
      dplyr::group_by(Fluor, Target, Biological_Set_Name, Sample) %>%
      dplyr::filter(any(!is.na(Cq)))

  } else {
    results <- results
  }

  final_df$results <- results %>% dplyr::ungroup()

  return(final_df)

}
