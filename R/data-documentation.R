#' Example qPCR Dataset
#'
#' Test qPCR data for different pathogens.
#'
#' This file is stored in inst/extdata
#'
#' @format A data frame with X rows and variables:
#' \describe{
#'   \item{Well}{Location on a 96 well plate}
#'   \item{Fluor}{Fluorophore used}
#'   \item{Target}{Target organism or gene of interest}
#'   \item{Sample}{Unique identifier for the sample}
#'   \item{Cq}{Cycle quantification output from qPCR machine}
#'   \item{SQ}{Starting quantity}
#'   \item{Biological Set Name}{Can be the same as sample, primer/probe set, multiplex, treatment type}
#' }
#' @source Jex Lab, WEHI
#' @name test_qPCR_data_csv
NULL
