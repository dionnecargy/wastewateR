#' Plot the standard curve
#'
#' @param qPCR_results  Output from `readqPCR()`.
#' @param MPlex Does your data contain multiplex data? Yes or No. Default: No.
#' @param Samples Does your data contain samples or just controls? Yes or No.
#'
#' @returns A table of the linear model formula, slope, R², and PCR efficiency,
#' and a plot of the standards and fit line.
#' @export
#'
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @importFrom stats lm sd
#' @importFrom broom  tidy glance
#'
#' @author Dionne Argyropoulos
plotStd <- function(qPCR_results, MPlex = "n", Samples = "n"){

  # User inputs for MPlex and Samples
  MPlex   <- tolower(trimws(MPlex)) # Normalize input: lowercase, trim spaces
  MPlex   <- substr(MPlex, 1, 1) # Only use first character (y/n)
  Samples <- tolower(trimws(Samples)) # Normalize input: lowercase, trim spaces
  Samples <- substr(Samples, 1, 1) # Only use first character (y/n)

  # Load Data
  df <- summariseqPCR(qPCR_results, Samples = Samples)

  # Step 1: Calculate formulas, slope, R², and PCR efficiency
  if(Samples == "y"){
    results_table <- df %>%
      dplyr::filter(str_detect(Sample, "STD")) %>%
      dplyr::group_by(Target, Plex) %>%
      dplyr::do({
        model           <- lm(mean ~ LogCopy, data = .)
        tidy_coef       <- broom::tidy(model)
        glance_stats    <- broom::glance(model)

        slope           <- tidy_coef$estimate[tidy_coef$term == "LogCopy"]
        intercept       <- tidy_coef$estimate[tidy_coef$term == "(Intercept)"]
        r_squared       <- glance_stats$r.squared
        pcr_efficiency  <- 10^(-1 / slope) - 1

        data.frame(
          slope = slope,
          intercept = intercept,
          r_squared = r_squared,
          r_sq_pass = ifelse(r_squared >= 0.9, "PASS", "FAIL"),
          pcr_efficiency = pcr_efficiency,
          pcr_pass = ifelse(pcr_efficiency >= 0.9 & pcr_efficiency <= 1.1, "PASS", "FAIL"),
          formula = paste0(
            "y = ", round(intercept, 2),
            ifelse(slope >= 0, " + ", " - "),
            round(abs(slope), 2), " * x"
          )
        )
      }) %>%
      dplyr::ungroup()
  } else {
    results_table <- df %>%
      dplyr::group_by(Target, Plex) %>%
      dplyr::do({
        model           <- lm(mean ~ LogCopy, data = .)
        tidy_coef       <- broom::tidy(model)
        glance_stats    <- broom::glance(model)

        slope           <- tidy_coef$estimate[tidy_coef$term == "LogCopy"]
        intercept       <- tidy_coef$estimate[tidy_coef$term == "(Intercept)"]
        r_squared       <- glance_stats$r.squared
        pcr_efficiency  <- 10^(-1 / slope) - 1

        data.frame(
          slope = slope,
          intercept = intercept,
          r_squared = r_squared,
          r_sq_pass = ifelse(r_squared >= 0.9, "PASS", "FAIL"),
          pcr_efficiency = pcr_efficiency,
          pcr_pass = ifelse(pcr_efficiency >= 0.9 & pcr_efficiency <= 1.1, "PASS", "FAIL"),
          formula = paste0(
            "y = ", round(intercept, 2),
            ifelse(slope >= 0, " + ", " - "),
            round(abs(slope), 2), " * x"
          )
        )
      }) %>%
      dplyr::ungroup()
  }

  # Step 2: Plot Standard Curve with Trendline
  if(Samples == "y" & MPlex == "n"){

    stdcurve_plot <- df %>%
      dplyr::filter(str_detect(Sample, "STD")) %>%
      dplyr::select(Fluor:Plate, starts_with("Rep")) %>%
      tidyr::pivot_longer(
        -c(Fluor:Plate),
        names_to = "Rep",
        values_to = "Cq"
      ) %>%
      drop_na(Cq) %>%
      ggplot2::ggplot(aes(LogCopy, Cq, colour = Target)) +
      ggplot2::geom_point() +
      ggplot2::geom_smooth(aes(group = Target), method = "lm", se = TRUE) +
      ggplot2::facet_wrap(~Target) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none")

  } else if (Samples == "n" & MPlex == "y"){

    stdcurve_plot <- df %>%
      dplyr::select(Fluor:Plate, starts_with("Rep")) %>%
      tidyr::pivot_longer(
        -c(Fluor:Plate),
        names_to = "Rep",
        values_to = "Cq"
      ) %>%
      drop_na(Cq) %>%
      ggplot2::ggplot(aes(LogCopy, Cq, colour = Target)) +
      ggplot2::geom_point() +
      ggplot2::geom_smooth(aes(group = Target), method = "lm", se = TRUE) +
      ggplot2::facet_wrap(Plex~Target) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none")

  } else if (Samples == "n" & MPlex == "n"){

    stdcurve_plot <- df %>%
      dplyr::select(Fluor:Plate, starts_with("Rep")) %>%
      tidyr::pivot_longer(
        -c(Fluor:Plate),
        names_to = "Rep",
        values_to = "Cq"
      ) %>%
      drop_na(Cq) %>%
      ggplot2::ggplot(aes(LogCopy, Cq, colour = Target)) +
      ggplot2::geom_point() +
      ggplot2::geom_smooth(aes(group = Target), method = "lm", se = TRUE) +
      ggplot2::facet_wrap(~Target) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none")

  } else {
    message("Please choose Samples and MPlex type.")
  }

  return(list(
    results_table = results_table,
    stdcurve_plot = stdcurve_plot
  ))
}
