years_between = function(start, end) {
	as.numeric(difftime(end, start, units = "days"))/365.24
}


chip_to_ra_survival = function(bl_data,
															 outcomes,
															 min_num_events = 100,
															 min_num_events_w_chip = 10,
															 ra_types = c('ra', 'spra', 'snra'),
															 sensspecs = c('sensitive', 'moderate', 'specific'),
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
					group_by(across(chip_type)) |>
					summarise(n_event = sum(!is.na(.data[[outcome_str]]))) ->
					event_counts

				if (sum(event_counts$n_event) < min_num_events) { next }
				if (min(event_counts$n_event) < min_num_events_w_chip) { next }

				str_form = glue('time_at_risk ~ {chip_type}')
				py = pyears(as.formula(str_form), df_for_surv)$pyears

				str_form = glue('Surv(time = time_at_risk, event = event) ~ age_at_biosample + age2 + sex + race + ever_smoker + {chip_type}')
				fit = coxph(as.formula(str_form), data = df_for_surv)

				# str_form = glue('event ~ age_at_biosample + age2 + sex + race + ever_smoker + {chip_type}')
				# fit = glm(formula = as.formula(str_form), data = df_for_surv, family = 'binomial')

				if (debug) { browser() }
				tidy(fit) |>
					mutate(q.value = qvalue(p.value, lambda = 0)$qvalues,
								 ra_type = ra_type,
								 sensspec = sensspec,
								 outcome = outcome_str,
								 n = nrow(df_for_surv),
								 n_chip = sum(df_for_surv[[chip_type]]),
								 n_event = sum(event_counts$n_event),
								 n_events_w_chip = min(event_counts$n_event),
								 time_at_risk_nochip = py['FALSE'],
								 time_at_risk_chip = py['TRUE']

					) |>
					filter(str_detect(term, 'chip')) ->
					toresult
				results[[paste0(outcome_str, '_', chip_type)]] = toresult
			}
		}
	}

	return(bind_rows(results))
}
