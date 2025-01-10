case_when_na = function(...) {
	case_when(..., .default = NA_Date_)
}

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
prepare_ra_outcomes = function(dx_code_counts_df,
															 serostatus_df,
															 ra_meds_df,
															 min_num_dx_codes = 2) {

	dx_code_counts_df |>
		filter(dx_count >= min_num_dx_codes) |>
		select(-dx_count) ->
		result

	result |>
		pivot_wider(id_cols = 'person_id',
								names_from = 'dx',
								names_glue = 'date_first_{dx}',
								values_from = 'date_first_dx') |>
		# if SPRA or SNRA was diagnosed before RA, the date of first RA
		# is moved up to date of first SNRA or SPRA (since those are types of RA)
		mutate(date_first_ra = pmin(date_first_ra,
																date_first_spra,
																date_first_snra,
																na.rm = TRUE)) |>
		left_join(serostatus_df, by = 'person_id') |>
		left_join(ra_meds_df, by = 'person_id') |>
		mutate(
			date_first_ra_or = case_when_na(
				!is.na(seropos) | treated_w_ra_meds ~ date_first_ra),

			date_first_spra_or = case_when_na(
				# for pts who quality based on labs, any RA dx counts as SPRA
				seropos ~ date_first_ra,
				# for pts who quality based on meds, only SPRA dx counts as SPRA
				treated_w_ra_meds ~ date_first_spra),

			date_first_snra_or = case_when_na(
				# for pts who quality based on labs, any RA dx counts as SNRA
				!seropos ~ date_first_ra,
				# for pts who quality based on meds, only SNRA dx counts as SNRA
				treated_w_ra_meds ~ date_first_snra),

			date_first_ra_and = case_when_na(
				!is.na(seropos) & treated_w_ra_meds ~ date_first_ra),

			date_first_spra_and = case_when_na(
				seropos & treated_w_ra_meds ~ date_first_ra),

			date_first_snra_and = case_when_na(
				!seropos & treated_w_ra_meds ~ date_first_ra)
		) |>
		select(-seropos, -treated_w_ra_meds) |>
		mutate(across(starts_with('date_first'), as.Date)) ->
		result

	return(result)
}

# Function to check if mutation is present for a given gene
check_gene_mutation <- function() {
	gene_name <- str_remove(cur_column(), "has_") |> toupper()
	return(chip_gene == gene_name)
}

# Function to create threshold variants for a mutation indicator
create_threshold_variants <- function(has_mutation, af_value, threshold) {
	case_when(
		has_mutation & af_value > threshold ~ TRUE,
		has_mutation ~ NA,
		.default = FALSE
	)
}

# Function to generate threshold variants for a column
apply_thresholds <- function(has_mutation_col) {
	map(thresholds, function(t) create_threshold_variants(has_mutation_col, AF, t))
}


#' @title prepare baseline data
#'
#' @param demographics tibble with baseline demographic data. See details.
#' @param chip_calls tibble with chip call status. See details
#' @param min_num_dx minimum number of diagnoses a patient must have to be studied
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
prepare_baseline_data = function(demographics,
																 chip_calls,
																 min_num_dx = 5,
																 thresholds = c("05" = 0.05, "10" = 0.10, "15" = 0.15, "20" = 0.20),
																 genes = c("DNMT3A", "TET2")) {

	demographics |>
		# filter out people with missing data
		filter(!is.na(person_id), !is.na(sex), !is.na(race)) |>
		filter(!is.na(ever_smoker), !is.na(biosample_date), !is.na(age_at_biosample)) |>
		filter(!is.na(age2), !is.na(date_last_dx)) |>
		# must have at least 5 dx codes observed
		filter(num_dx >= min_num_dx) |>
		# must have some observation time after biosample_date
		filter(date_last_dx > biosample_date) ->
		cohort_for_study

	cohort_for_study |>
		left_join(chip_calls, by = 'person_id') |>
		filter(is.na(date_first_heme_ca) | date_first_heme_ca > biosample_date) |>
		mutate(
			censor_date = date_last_dx,
			has_chip = !is.na(chip_gene),
			# Create base gene indicators
			across(
				all_of(paste0("has_", tolower(genes))),
				check_gene_mutation,
				.names = "{.col}"
			)
		) |>
		# Add threshold variations for CHIP and each gene
		mutate(
			across(
				starts_with("has_"),
				apply_thresholds,
				.names = "{.col}_{names(thresholds)}"
			)
		) |>
		select(-date_last_dx) ->
		result

	return(result)
}
