years_between = function(start, end) {
	as.numeric(difftime(end, start, units = "days"))/365.24
}

# TODO
chip_to_ra_survival = function(bl_data,
															 outcomes,
															 ra_types = c('ra', 'spra', 'snra'),
															 sensspecs = c('sensitive', 'moderate', 'specific'),
															 chip_types = c('has_chip', 'has_chip10',
															 							 'has_chip_d3a', 'has_chip_tet2',
															 							 'has_chip_asxl1'),
															 debug = FALSE) {

	bl_data |>
		left_join(outcomes, by = 'person_id') ->
		df_for_survprep

	results = list()

	for (ra_type in ra_types) {

		message('*** Starting ', ra_type, '...')

		for (sensspec in sensspecs) {

			message('**Starting ', sensspec, '...')
			outcome_str = glue('date_first_{ra_type}_{sensspec}')

			df_for_survprep |>
				# remove people with prevalent disease
				filter(is.na(.data[[outcome_str]]) | .data[[outcome_str]] > biosample_date) |>
				# event = had the disease before censoring (not really necessary)
				mutate(event = !is.na(.data[[outcome_str]]) & .data[[outcome_str]] < censor_date) |>
				# time at risk definition
				mutate(time_at_risk = case_when(event ~ years_between(biosample_date, .data[[outcome_str]]),
																				.default = years_between(biosample_date, censor_date))) |>
				# this shouldn't really be necessary...
				filter(time_at_risk > 0) ->
				df_for_surv

			for(chip_type in chip_types) {

				message('Starting ', chip_type, '...')

				# TODO: calculate time at risk, stratified by chip_type, and store that

				str_form = glue('Surv(time = time_at_risk, event = event) ~ age_at_biosample + age2 + sex + race + ever_smoker + {chip_type}')
				fit = coxph(as.formula(str_form), data = df_for_surv)
				if (debug) { browser() }
				tidy(fit) |>
					mutate(ra_type = ra_type,
								 sensspec = sensspec,
								 outcome = outcome_str,
								 n = nrow(df_for_surv),
								 n_event = sum(df_for_surv$event),
								 n_chip = sum(df_for_surv[[chip_type]])) ->
					toresult
				results[[paste0(outcome_str, '_', chip_type)]] = toresult
			}
		}
	}

	return(bind_rows(results))
}
