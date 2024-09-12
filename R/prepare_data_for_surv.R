#' @title Prepare RA Outcomes
#'
#' @param dx_code_counts_df tibble of diagnosis codes. See details.
#' @param seropositive_df tibble of sero-status. See details.
#' @param meds_df tibble of medications. See details.
#'
#' @details
#' dx_code_counts_df is a tibble with the four columns: person_id, ra_dx
#' ra_dx_count, date_first_dx. The values in ra_dx can be 'RA', 'SPRA', or
#' 'SNRA'. ra_dx_count is an integer indicating how many times that diagnosis
#' code was seen on that patient. date_first_dx is the date that diagnosis was
#' first seen on that patient.
#'
#' seropositive_df is a tibble with two columns: person_id, and seropos.
#' The values in seropos can be TRUE and FALSE.
#'
#' meds_df is a tibble with two columns: person_id and treated_w_ra_meds.
#' The values in treated_w_ra_meds can be TRUE and FALSE.
#'
#'
#' @return
#' @export
#'
prepare_ra_outcomes = function(dx_code_counts_df, seropositive_df, meds_df) {


}
