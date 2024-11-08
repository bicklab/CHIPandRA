#' @title years between
#'
#' @param start start date
#' @param end end date
#'
#' @return the number of days between the two dates
#' @export
#'
years_between = function(start, end) {
	as.numeric(difftime(end, start, units = "days"))/365.24
}


#' @title CHIP to RA survival analysis
#'
#' @param bl_data baseline data
#' @param outcomes outcome data
#' @param min_num_events mnimum number of events (outcomes) needed to be observed to continue with analysis
#' @param min_num_events_w_chip minimum number of events (outcomes) needed to be onbserved among people with CHIP to continue analysis
#' @param ra_types types of RA to study
#' @param sensspecs spectrum of sensitive to specific to study
#' @param chip_types types of CHIP to study
#' @param debug flag to debug vs run normally
#'
#' @return results tibble
#' @export
#'
chip_to_ra_survival = function(bl_data,
															 outcomes,
															 min_num_events = 10,
															 min_num_events_w_chip = 1,
															 ra_types = c('ra', 'spra', 'snra'),
															 sensspecs = c('or', 'and'),
															 chip_types =  names(select(bl_data, starts_with('has_chip'))),
															 debug = FALSE) {

	bl_data |>
		left_join(outcomes, by = 'person_id') ->
		df_for_survprep

	results = list()

	for (ra_type in ra_types) {

		message('\n*** Starting ', ra_type, '...')

		for (sensspec in sensspecs) {

			message('** Starting ', sensspec, '...')
			outcome_str = glue('date_first_{ra_type}_{sensspec}')

			df_for_survprep |>
				# remove people with prevalent disease
				filter(is.na(.data[[outcome_str]]) | .data[[outcome_str]] > biosample_date) |>
				# event = had the disease before censoring
				mutate(event = !is.na(.data[[outcome_str]]) & .data[[outcome_str]] < censor_date) |>
				# time at risk definition
				mutate(time_at_risk = case_when(event ~ years_between(biosample_date, .data[[outcome_str]]),
																				.default = years_between(biosample_date, censor_date))) ->
				df_for_surv

			for(chip_type in chip_types) {

				message('Starting ', chip_type, '...')

				df_for_surv |>
					filter(if_any({{chip_type}}, ~!is.na(.))) |>
					group_by(across(chip_type)) |>
					summarise(
						num_pts = n(),
						total_fu = sum(time_at_risk),
						mean_fu = total_fu / num_pts,
						median_fu = median(time_at_risk),
						n_event = sum(!is.na(.data[[outcome_str]])),
						events_per_py = n_event / total_fu) |>
					arrange({{chip_type}}) ->
					event_counts

				if (sum(event_counts$n_event) < min_num_events) { next }
				if (min(event_counts$n_event) < min_num_events_w_chip) { next }

				if (debug) { browser() }

				str_form = glue('Surv(time = time_at_risk, event = event) ~ age_at_biosample + age2 + sex + race + ever_smoker + {chip_type}')
				fit = tryCatch(
					coxph(as.formula(str_form), data = df_for_surv),
					finally = NULL
				)

				if (debug) { browser() }

				tidy(fit) |>
					filter(str_detect(term, 'chip')) |>
					mutate(ra_type = ra_type,
								 sensspec = sensspec,
								 outcome = outcome_str,
								 num_prevalent_dz = nrow(df_for_survprep) - nrow(df_for_surv),

								 num_pts_no_chip = event_counts[['num_pts']][1],
								 num_pts_small_chip = event_counts[['num_pts']][2],
								 num_pts_big_chip = event_counts[['num_pts']][3],

								 total_fu_no_chip = event_counts[['total_fu']][1],
								 total_fu_small_chip = event_counts[['total_fu']][2],
								 total_fu_big_chip = event_counts[['total_fu']][3],

								 mean_fu_no_chip = event_counts[['mean_fu']][1],
								 mean_fu_small_chip = event_counts[['mean_fu']][2],
								 mean_fu_big_chip = event_counts[['mean_fu']][3],

								 median_fu_no_chip = event_counts[['median_fu']][1],
								 median_fu_small_chip = event_counts[['median_fu']][2],
								 median_fu_big_chip = event_counts[['median_fu']][3],

								 num_events_no_chip = event_counts[['n_event']][1],
								 num_events_small_chip = event_counts[['n_event']][2],
								 num_events_big_chip = event_counts[['n_event']][3],

								 event_rate_no_chip = event_counts[['events_per_py']][1],
								 event_rate_small_chip = event_counts[['events_per_py']][2],
								 event_rate_big_chip = event_counts[['events_per_py']][3])  ->
					toresult
				results[[paste0(outcome_str, '_', chip_type)]] = toresult
			}
		}
	}

	return(bind_rows(results))
}
