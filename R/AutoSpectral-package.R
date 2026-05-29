#' @title AutoSpectral: Tools for Unmixing Spectral Flow Cytometry Data
#'
#' @description
#' As both a refinement and advancement of the unmixing process of full spectrum cytometry data, `AutoSpectral` provides a suite of functions/tools that work in concert to provide the user with optimal unmixing results.
#'
#' Use `AutoSpectral` to:
#' \itemize{
#'   \item isolate and refine clean spectral signatures
#'   \item mitigate autofluorescence contamination in complex controls
#'   \item unmix fcs files using standard algorithms or the hallmark `AutoSpectral` approach
#' }
#'
#' @section `AutoSpectral` (v 1.5.7) Workflow:
#' To maximize the success of `AutoSpectral`, the workflow is organized into the following logical step-wise processes:
#' 1. [create.control.file]: Generate a .csv file that contains descriptive information about your single stained controls
#'     * [Control File](https://drcytometer.github.io/AutoSpectral/articles/01_Full_AutoSpectral_Workflow.html#creating-the-control-file) (see article)
#' 1. [define.flow.control]: uses the result of [create.control.file] to convert FCS files and their associated metadata into an optimized data structure for downstream workflows
#'    *  [Flow Control](https://drcytometer.github.io/AutoSpectral/articles/01_Full_AutoSpectral_Workflow.html#loading-the-data) (see article)
#' 1. ...
#' 1. ...
#'
#' @name AutoSpectral
"_PACKAGE"
