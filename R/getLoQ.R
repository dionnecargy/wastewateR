#' Get Limit of Quantification
#'
#' Based on the method described in Schang et al., 2021 Environ. Sci. Technol.
#' where the LOQ was estimated from the replicated standard curves and standard
#' dilutions. The coefficient of variation (CV) for log-normal distributed
#' datasets is calculated using the equation:
#'
#' \deqn{
#' CV = \sqrt{ (1 + E)^{SD^2 \times \ln(1+E)} - 1 }
#' }
#'
#' where E is the PCR efficiency and SD is the standard deviation across
#' replicate Cq values for a stnadard curve dilution.
#'
#' @param qPCR_results  Output from `readqPCR()`.
#' @param MPlex Does your data contain multiplex data? Yes or No. Default: No.
#' @param Samples Does your data contain samples or just controls? Yes or No.
#' @param cv_threshold Decimal number. Default = 0.35 (35%).
#'
#' @returns A table of LoQ results and a plot of the CV and log-copy number.
#' @export
#'
#' @author Dionne Argyropoulos
getLoQ <- function(qPCR_results, MPlex = "n", Samples = "n", cv_threshold = 0.35){

  # User inputs for MPlex and Samples
  MPlex   <- tolower(trimws(MPlex)) # Normalize input: lowercase, trim spaces
  MPlex   <- substr(MPlex, 1, 1) # Only use first character (y/n)
  Samples <- tolower(trimws(Samples)) # Normalize input: lowercase, trim spaces
  Samples <- substr(Samples, 1, 1) # Only use first character (y/n)

  # Load Data
  df      <- summariseqPCR(qPCR_results, Samples = Samples)
  std     <- plotStd(qPCR_results, MPlex = MPlex, Samples = Samples)

  # 2. Join PCR efficiency to SD table
  if (Samples == "y"){
    df_joined <- df %>%
      dplyr::filter(str_detect(Sample, "STD")) %>%
      dplyr::left_join(
        std$results_table %>% dplyr::select(Target, Plex, E = pcr_efficiency),
        by = c("Target", "Plex")
      )
  } else {
    df_joined <- df %>%
      dplyr::left_join(
        std$results_table %>% dplyr::select(Target, Plex, E = pcr_efficiency),
        by = c("Target", "Plex")
      )
  }

  # 3. Calculate CV
  df_joined <- df_joined %>%
    mutate(
      # Calculate CV using PCR efficiency
      CV = sqrt((1 + E)^(sd^2 * log(1 + E)) - 1),
      # Create a flag for passing CV threshold
      CV_pass = ifelse(CV <= cv_threshold, TRUE, FALSE)
    )

  LoQ_table <- df_joined %>%
    dplyr::arrange(Target, Plex, LogCopy) %>%
    dplyr::group_by(Target, Plex) %>%
    dplyr::group_modify(~{
      d <- .x

      # Identify failing and passing points
      d <- d %>% dplyr::mutate(pass = CV < cv_threshold)

      # If all pass, LOQ is the lowest concentration
      if (all(d$pass)) {
        return(tibble::tibble(
          LoQ = min(10^d$LogCopy),
          LoQ_LogCopy = min(d$LogCopy),
          note = "All CV < threshold"
        ))
      }

      # Find first passing point after a sequence of failing points
      fail_idx <- max(which(!d$pass))
      pass_idx <- fail_idx + 1

      # If no passing point afterwards → LOQ undefined
      if (pass_idx > nrow(d)) {
        return(tibble::tibble(
          LoQ = NA,
          LoQ_LogCopy = NA,
          note = "No passing CV above threshold"
        ))
      }

      # Extract points for interpolation
      x1 <- d$LogCopy[fail_idx]
      y1 <- d$CV[fail_idx]
      x2 <- d$LogCopy[pass_idx]
      y2 <- d$CV[pass_idx]

      # Linear interpolation to find CV = threshold
      log_loq <- x1 + (cv_threshold - y1) * (x2 - x1) / (y2 - y1)
      loq <- 10^log_loq

      tibble(
        LoQ = loq,
        LoQ_LogCopy = log_loq,
        note = "Interpolated"
      )
    }) %>%
    ungroup()

  # Plot LoQ
  df_plot <- df_joined %>% dplyr::left_join(LoQ_table, by = c("Target", "Plex"))

  LoQ_plot <- df_plot %>%
    ggplot2::ggplot(aes(x = LogCopy, y = CV*100)) +
    ggplot2::geom_point() +
    ggplot2::geom_hline(yintercept = 35, linetype = "dashed") +
    ggplot2::geom_vline(aes(xintercept = LoQ_LogCopy), colour = "red") +
    ggplot2::facet_grid(Target ~ Plex) +
    ggplot2::labs(
      x = "log10(Copy Number)",
      y = "Coefficient of Variation (CV, %)",
      title = "LOQ Estimation from CV Threshold"
    ) +
    ggplot2::theme_bw()

  return(
    list(
      LoQ_table = LoQ_table,
      LoQ_plot = LoQ_plot
    )
  )

}
