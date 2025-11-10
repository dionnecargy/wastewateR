#' Read raw qPCR Data
#'
#' A function to read in raw qPCR data
#'
#' @param raw_data One file path (.csv)
#' @param MPlex Does your data contain multiplex data? Yes or No. Default: No.
#' @param Samples Does your data contain samples or just controls? Yes or No.
#' @param Target Does the "Target" column have data in it? Either supply another
#' df with the Sample and Target or leave blank. Default is NULL.
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
readqPCR <- function(raw_data, MPlex = "No", Samples, Target = NULL){

  # Write raw data file name and extract "plate"
  raw_data_filenames <- gsub("\\s+", "", tolower(basename(raw_data)))

  # User inputs for MPlex and Samples
  MPlex   <- tolower(trimws(MPlex)) # Normalize input: lowercase, trim spaces
  MPlex   <- substr(MPlex, 1, 1) # Only use first character (y/n)
  Samples <- tolower(trimws(Samples)) # Normalize input: lowercase, trim spaces
  Samples <- substr(Samples, 1, 1) # Only use first character (y/n)

  # Initialise master list to store files
  master_list <- list(
    results   = NULL,  # Placeholder for processed results combined
    ntcs      = NULL   # Placeholder for any ntcs data combined
  )

  # Loop through each file and process accordingly
  for (i in seq_along(raw_data)) {
    file <- raw_data[i]
    file_name <- raw_data_filenames[i]

    #############################################################################
    # Data processing
    #############################################################################
    final_df <- list()
    df <- read.csv(here::here(file), na.strings = c("NaN"))
    # Store row numbers
    run_info_end_row <- which(df$File.Name == "Well")
    final_row <- nrow(df)
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

    #############################################################################
    # Does your data contain a "Targets" column that is filled out?
    #############################################################################
    results <- results %>% dplyr::mutate(Target = na_if(trimws(Target), ""))
    # Check if "Target" column is empty in dataset and needs to be supplied
    target_col <- any(is.na(results$Target))
    if (target_col == TRUE){
      results <- results %>%
        dplyr::select(-Target) %>%
        dplyr::left_join(target, by = "Sample")
    } else {
      results <- results
    }

    #############################################################################
    # Filter NTCs
    #############################################################################
    ntcs <- results %>%
      dplyr::ungroup() %>%
      dplyr::mutate(Cq = as.numeric(ifelse(Cq %in% c("NaN", "", " "), NA, Cq))) %>%
      dplyr::filter(Content == "NTC" | Sample == "NTC")

    ntc_check <- ntcs %>%
      dplyr::summarise(all_na = all(is.na(Cq))) %>%
      dplyr::pull(all_na)

    if (ntc_check) {
      message("All NTCs are NaN")
    } else {
      stop("Check NTCs")
    }

    # Remove NTCs from results once confirmed they are all NaN
    results <- results %>%
      dplyr::filter(!(Content == "NTC" | Sample == "NTC"))

    #############################################################################
    # Does your data contain multiplex or simplex ?
    #############################################################################
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

    #############################################################################
    # Does your data contain samples and standards?
    #############################################################################
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

      results <- results %>%
        ungroup() %>%
        dplyr::mutate(
          SQ = format(10 ^ as.numeric(str_remove(Biological_Set_Name, "10\\^")), scientific = FALSE, trim = TRUE, drop0trailing = TRUE),
          SQ = ifelse(SQ == "NA", NA, SQ ),
          Sample = ifelse(Sample == "Control", paste0("STD_", SQ), Sample),
          LogCopy = log10(as.numeric(SQ))
        ) %>%
        dplyr::select(Well, Fluor, Sample, Cq, SQ, LogCopy, Target)

    } else {
      # Remove Sample column
      results <- results %>%
        dplyr::mutate(
          SQ = format(10 ^ as.numeric(str_remove(Biological_Set_Name, "10\\^")), scientific = FALSE, trim = TRUE, drop0trailing = TRUE),
          Sample = paste0("STD_", SQ),
          LogCopy = log10(as.numeric(SQ))
        ) %>%
        dplyr::select(Well, Fluor, Sample, Cq, SQ, LogCopy, Target) %>%
        dplyr::select(-Sample)

      ntcs <- ntcs %>% dplyr::select(-Sample)
    }

    #############################################################################
    # Add to master_list
    #############################################################################
    # Add 'plate_id' column to each dataframe
    raw_data_filenames_info <- tibble(file_name) %>%
      mutate(
        date  = str_extract(file_name, "\\d{8}"),
        plate = str_extract(file_name, "(?<=plate)(?=\\d)\\d+"),  # plate followed by digits
        test  = str_extract(file_name, "(?<=test)(?=\\d)\\d+"),   # test followed by digits
        plate_id = case_when(
          !is.na(date) & !is.na(plate) & !is.na(test) ~ paste0(date, "_plate", plate, "_test", test),
          !is.na(date) & !is.na(plate)                ~ paste0(date, "_plate", plate),
          !is.na(date) & !is.na(test)                 ~ paste0(date, "_test", test),
          !is.na(plate) & !is.na(test)                ~ paste0("plate", plate, "_test", test),
          !is.na(date)                                ~ date,
          !is.na(plate)                               ~ paste0("plate", plate),
          !is.na(test)                                ~ paste0("test", test),
          TRUE                                        ~ NA_character_
        )
      )
    results$Plate   <- raw_data_filenames_info$plate_id
    ntcs$Plate      <- raw_data_filenames_info$plate_id

    # Add processed file's tables to the master list
    master_list$results   <- dplyr::bind_rows(master_list$results, results)     # Combine processed results
    master_list$ntcs      <- dplyr::bind_rows(master_list$ntcs, ntcs)           # Combine ntcs

  }

 return(master_list)

}
