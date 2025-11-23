#' Limit of Detection
#'
#' Calculate the limit of detection at a threshold (discrete calculation).
#'
#' @param qPCR_results  Output from `readqPCR()`.
#' @param MPlex Does your data contain multiplex data? Yes or No. Default: No.
#' @param Samples Does your data contain samples or just controls? Yes or No.
#' @param lod_threshold Limit of detection threshold. Default = 0.95 (95%).
#'
#' @returns A table and a plot the LoD results.
#' @export
#'
#' @author Dionne Argyropoulos
getLoD <- function(qPCR_results, MPlex = "n", Samples = "n", lod_threshold = 0.95){

  # User inputs for MPlex and Samples
  MPlex   <- tolower(trimws(MPlex)) # Normalize input: lowercase, trim spaces
  MPlex   <- substr(MPlex, 1, 1) # Only use first character (y/n)
  Samples <- tolower(trimws(Samples)) # Normalize input: lowercase, trim spaces
  Samples <- substr(Samples, 1, 1) # Only use first character (y/n)

  # filter dataset
  df <- qPCR_results[[1]]

  # calculate detections: 1 = detected, 0 = not detected
  if(Samples == "y"){
    detections <- df %>%
      dplyr::filter(str_detect(Sample, "STD")) %>%
      dplyr::mutate(detected = ifelse(!is.na(Cq), 1, 0)) %>%
      dplyr::group_by(SQ, Target, LogCopy) %>%
      summarise(
        n_detected = sum(detected, na.rm = TRUE),
        n_total = n(),
        rate = n_detected / n_total,
        .groups = "drop"
      )

  } else if (MPlex == "y") {
    detections <- df %>%
      dplyr::mutate(detected = ifelse(!is.na(Cq), 1, 0)) %>%
      dplyr::group_by(SQ, Target, LogCopy, Plex) %>%
      summarise(
        n_detected = sum(detected, na.rm = TRUE),
        n_total = n(),
        rate = n_detected / n_total,
        .groups = "drop"
      )
  } else {
    detections <- df %>%
      dplyr::mutate(detected = ifelse(!is.na(Cq), 1, 0)) %>%
      dplyr::group_by(SQ, Target, LogCopy) %>%
      summarise(
        n_detected = sum(detected, na.rm = TRUE),
        n_total = n(),
        rate = n_detected / n_total,
        .groups = "drop"
      )
  }

  # calculate LoD for 95% analytical sensitivity
  if (MPlex == "y"){
    lod <- detections %>%
      group_by(Target, Plex) %>%
      arrange(LogCopy) %>%   # make sure it's in order
      summarise(
        LoD = LogCopy[which(rate >= lod_threshold)[1]],  # first LogCopy >= 0.95
        .groups = "drop"
      )
  } else {
    lod <- detections %>%
      group_by(Target) %>%
      arrange(LogCopy) %>%   # make sure it's in order
      summarise(
        LoD = LogCopy[which(rate >= lod_threshold)[1]],  # first LogCopy >= 0.95
        .groups = "drop"
      )
  }


  # plot LoD detection rates
  if (MPlex == "y"){
    LoD_plot <- detections %>%
      ggplot2::ggplot(aes(x = LogCopy, y = rate)) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(Plex~Target, scales = "free_x") +
      ggplot2::theme_bw() +
      # Horizontal line at 95% detection
      ggplot2::geom_hline(yintercept = lod_threshold, linetype = "dashed", color = "red") +
      # Vertical line for LoD per Target
      ggplot2::geom_vline(
        data = lod,
        aes(xintercept = LoD),
        linetype = "dashed",
        color = "red"
      ) +
      ggplot2::labs(
        x = "Log Copy Number",
        y = "Detection Rate",
        title = "Limit of Detection per Target"
      ) +
      ggplot2::lims(x = c(min(detections$LogCopy), max(detections$LogCopy)), y = c(0, 1))

  } else {
    LoD_plot <- detections %>%
      ggplot2::ggplot(aes(x = LogCopy, y = rate)) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(~Target, scales = "free_x") +
      ggplot2::theme_bw() +
      # Horizontal line at 95% detection
      ggplot2::geom_hline(yintercept = lod_threshold, linetype = "dashed", color = "red") +
      # Vertical line for LoD per Target
      ggplot2::geom_vline(
        data = lod,
        aes(xintercept = LoD),
        linetype = "dashed",
        color = "red"
      ) +
      ggplot2::labs(
        x = "Log Copy Number",
        y = "Detection Rate",
        title = "Limit of Detection per Target"
      ) +
      ggplot2::lims(x = c(min(detections$LogCopy), max(detections$LogCopy)), y = c(0, 1))

  }

  return(
    list(
      LoD_table = lod,
      LoD_plot = LoD_plot
    )
  )
}
