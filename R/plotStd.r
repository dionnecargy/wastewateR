#' Plot the standard curve
#'
#' @param stats_df  Output from `summariseqPCR()`
#'
#' @returns A table of the linear model formula, slope, R², and PCR efficiency,
#' and a plot of the standards and fit line.
#' @export
#'
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @importFrom stats lm na.omit sd
#'
#' @author Dionne Argyropoulos
plotStd <- function(stats_df){
  # Step 1: Calculate formulas, slope, R², and PCR efficiency
  results_table <- stats_df %>%
    tidyr::drop_na() %>%
    dplyr::group_by(Sample) %>%
    dplyr::do({
      model <- lm(mean ~ LogCopy, data = .)
      tidy_coef <- tidy(model)
      glance_stats <- glance(model)

      slope <- tidy_coef$estimate[tidy_coef$term == "LogCopy"]
      intercept <- tidy_coef$estimate[tidy_coef$term == "(Intercept)"]
      r_squared <- glance_stats$r.squared
      pcr_efficiency <- 10^(-1 / slope) - 1

      data.frame(
        slope = slope,
        intercept = intercept,
        r_squared = r_squared,
        r_sq_pass = ifelse(r_squared >= 0.9, "PASS", "FAIL"),
        pcr_efficiency = pcr_efficiency,
        pcr_pass = ifelse(pcr_efficiency >= 0.9 & pcr_efficiency <= 1.1, "PASS", "FAIL"),
        formula = paste0("y = ", round(intercept, 2),
                         ifelse(slope >= 0, " + ", " - "),
                         round(abs(slope), 2), " * x")
      )
    }) %>%
    dplyr::ungroup()

  # Step 2: Plot Standard Curve with Trendline
  stdcurveplot <- stats_df %>%
    ggplot2::ggplot(aes(LogCopy, mean, colour = Sample)) +
    ggplot2::geom_point() +
    # ggplot2::scale_colour_manual(values = melbourne::melb_trams(length(unique(stats_df$Sample)))) +
    ggplot2::geom_smooth(aes(group = Sample), method = "lm", se = TRUE) +
    ggplot2::facet_wrap(~Sample) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")

  return(list(stdcurve = stdcurveplot, model = results_table))
}
