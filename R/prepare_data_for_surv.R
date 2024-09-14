#' @title Prepare RA Outcomes
#'
#' @param dx_code_counts_df tibble of diagnosis codes. See details.
#' @param serostatus_df tibble of sero-status. See details.
#' @param ra_meds_df tibble of medications. See details.
#' @param min_num_dx_codes minimum number of occurrences of a specific
#' diagnosis code for a patient to be considered to have that disease
#'
#' @details
#' dx_code_counts_df is a tibble with the four columns: person_id, dx,
#' dx_count, date_first_dx. The values in dx can be 'RA', 'SPRA', or
#' 'SNRA'. dx_count is an integer indicating how many times that diagnosis
#' code was seen on that patient. date_first_dx is the date that diagnosis was
#' first seen on that patient.
#'
#' serostatus_df is a tibble with two columns: person_id, and seropos.
#' The values in seropos can be TRUE and FALSE.
#'
#' ra_meds_df is a tibble with two columns: person_id and treated_w_ra_meds.
#' The values in treated_w_ra_meds can be TRUE and FALSE.
#'
#' @return tibble with ten columns: person_id and nine RA phenotype
#' definitions, three for RA, three for SPRA, and three for SNRA. The three are
#' (1) sensitive, based on diagnosis codes only, (2) moderate,based on
#' diagnosis code and either labs or meds, and (3) specific, based on diagnosis
#' codes, labs, and meds
#'
#' @export
#'
prepare_ra_outcomes = function(dx_code_counts_df,
															 serostatus_df,
															 ra_meds_df,
															 min_num_dx_codes) {

	dx_code_counts_df |>
		filter(dx_count >= min_num_dx_codes) |>
		select(-dx_count) ->
		result

	result |>
		pivot_wider(id_cols = 'person_id',
								names_from = 'dx',
								names_glue = 'date_first_{dx}_sensitive',
								values_from = 'date_first_dx') |>
		left_join(serostatus_df, by = 'person_id') |>
		left_join(ra_meds_df, by = 'person_id') |>
		mutate(
			date_first_ra_moderate = ifelse(!is.na(seropos) | treated_w_ra_meds, date_first_ra_sensitive, NA_Date_),
			date_first_spra_moderate = ifelse(seropos | treated_w_ra_meds, date_first_spra_sensitive, NA_Date_),
			date_first_snra_moderate = ifelse((!seropos) | treated_w_ra_meds, date_first_snra_sensitive, NA_Date_),
			date_first_ra_specific = ifelse(!is.na(seropos) & treated_w_ra_meds, date_first_ra_sensitive, NA_Date_),
			date_first_spra_specific = ifelse(seropos & treated_w_ra_meds, date_first_spra_sensitive, NA_Date_),
			date_first_snra_specific = ifelse((!seropos) & treated_w_ra_meds, date_first_snra_sensitive, NA_Date_)) |>
		select(-seropos, -treated_w_ra_meds) |>
		mutate(across(starts_with('date_first'), as.Date)) ->
		result

	return(result)
}


#' @title prepare baseline data
#'
#' @param demographics tibble with baseline demographic data. See details.
#' @param chip_calls tibble with chip call status. See details
#'
#' @return baseline data prepped for survival analysis
#'
#' @details demographics includes: person_id (character), sex (factor),
#' race (factor), ever_smoker (logical), biosample_date (date),
#' age_at_biosample (double), age2 (double), date_first_heme_ca
#' (date or NA), and date_last_dx (date or NA).
#'
#' chip_calls includes: person_id (character), chip_gene (string or NA), and
#' AF (double or NA).
#'
prepare_baseline_data = function(demographics, chip_calls) {

	demographics |>
		# filter out people with missing data
		filter(!is.na(person_id), !is.na(sex), !is.na(race)) |>
		filter(!is.na(ever_smoker), !is.na(biosample_date), !is.na(age_at_biosample)) |>
		filter(!is.na(age2), !is.na(date_last_dx)) |>
		# must have some observation time after biosample_date
		filter(date_last_dx > biosample_date) ->
		cohort_for_study

	cohort_for_study |>
		left_join(chip_calls, by = 'person_id') |>
		filter(is.na(date_first_heme_ca) | date_first_heme_ca > biosample_date) |>
		mutate(
			censor_date = date_last_dx,
			has_chip = !is.na(chip_gene),
			has_chip10 = has_chip & AF > 0.10,
			has_chip_d3a = has_chip & chip_gene == 'DNMT3A',
			has_chip_tet2 = has_chip & chip_gene == 'TET2',
			has_chip_asxl1 = has_chip & chip_gene == 'ASXL1') |>
		select(-date_last_dx) ->
		result

	return(result)
}
