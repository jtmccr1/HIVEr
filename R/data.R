#' A small meta data file
#'
#' A small meta data file containing a few samples that
#' exihibit what we want to do in the code and a few intereseting
#' cases. Each row is a sample.
#'
#' @format A data frame with 20 rows and 17 variables:
#' \describe{
#'   \item{X}{Useless should remove}
#'   \item{HOUSE_ID}{House identifier of sample}
#'   \item{ENROLLILD}{personal identifier}
#'   \item{SPECID}{sample identifier}
#'   \item{onset}{Date of symptom onset}
#'   \item{collect}{date of sample collection}
#'   \item{vaccination_status}{vaccine status}
#'   \item{pcr_result}{strain H1N1 H3N2}
#'   \item{LAURING_ID}{Old sample id}
#'   \item{DPI}{Days post infection}
#'   \item{season}{Season of illness}
#'   \item{log_copy_num}{Log copy number of titer in qpcr well}
#'   \item{gc_ul}{genomes per ul in sample isolate}
#'   \item{HIGHSD}{Y or N for high standard deviation in qpcr}
#'   \item{sequenced}{T/F was the sample sequenced}
#'   \item{home_collected}{T/F was the sample taken at home}
#'   \item{snv_qualified}{T/F does the sample qualify for snv calling}
#'
#'
#'
#' }
#' @source data-raw/small_meta.R
"small_meta"