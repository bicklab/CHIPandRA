#' @title Do metaanalysis
#'
#' @details
#' functions to carry out meta analysis combining three biobanks
#'
#' @return none, just side effects of tables and plots
#' @export
#'
do_meta_analysis = function() {

	data_dir = '~/OneDrive - VUMC/CHIP and RA/'

	read_csv(file.path(data_dir, 'ukb_results/ukb_chip_ra_results.csv')) |>
		mutate(cohort = 'ukb') |>
		glimpse() ->
		ukb

	ukb |> one_cohort_all_results_plot(title = 'UKB: hazard ratio for incident RA', p_or_q = 'p')
	ggsave(file.path(data_dir, 'ukb_results', 'ukb_chip_ra_forest_plot_pvals.pdf'),
				 height = 9, width = 9)

	ukb |> one_cohort_all_results_plot(title = 'UKB: hazard ratio for incident RA', p_or_q = 'q')
	ggsave(file.path(data_dir, 'ukb_results', 'ukb_chip_ra_forest_plot_qvals.pdf'),
				 height = 9, width = 9)

	ukb |> one_cohort_specific_results_plot(title = 'UKB: hazard ratio for incident RA')
	ggsave(file.path(data_dir, 'ukb_results', 'ukb_chip_ra_focused_forest_plot_pvals.pdf'),
				 height = 3, width = 9)



	read_csv(file.path(data_dir, 'aou_results/aou_chip_ra_results.csv')) |>
		mutate(cohort = 'aou') |>
		glimpse() ->
		aou

	aou |> one_cohort_all_results_plot(title = 'AOU: hazard ratio for incident RA', p_or_q = 'p')
	ggsave(file.path(data_dir, 'aou_results', 'aou_chip_ra_forest_plot_pvals.pdf'),
				 height = 9, width = 9)

	aou |> one_cohort_all_results_plot(title = 'AOU: hazard ratio for incident RA', p_or_q = 'q')
	ggsave(file.path(data_dir, 'aou_results', 'aou_chip_ra_forest_plot_qvals.pdf'),
				 height = 9, width = 9)

	aou |> one_cohort_specific_results_plot(title = 'AOU: hazard ratio for incident RA')
	ggsave(file.path(data_dir, 'aou_results', 'aou_chip_ra_focused_forest_plot_pvals.pdf'),
				 height = 3, width = 9)



	read_csv(file.path(data_dir, 'biovu_results/biovu_chip_ra_results.csv')) |>
		mutate(cohort = 'biovu') |>
		glimpse() ->
		biovu

	biovu |> one_cohort_all_results_plot(title = 'BioVU: hazard ratio for incident RA', p_or_q = 'p')
	ggsave(file.path(data_dir, 'biovu_results', 'biovu_chip_ra_forest_plot_pvals.pdf'),
				 height = 9, width = 9)

	biovu |> one_cohort_all_results_plot(title = 'BioVU: hazard ratio for incident RA', p_or_q = 'q')
	ggsave(file.path(data_dir, 'biovu_results', 'biovu_chip_ra_forest_plot_qvals.pdf'),
				 height = 9, width = 9)

	biovu |> one_cohort_specific_results_plot(title = 'BioVU: hazard ratio for incident RA')
	ggsave(file.path(data_dir, 'biovu_results', 'biovu_chip_ra_focused_forest_plot_pvals.pdf'),
				 height = 3, width = 9)


	bind_rows(ukb, aou, biovu) |>
		mutate(cohort = factor(cohort, levels = c('ukb', 'aou', 'biovu')),
					 ra_type = factor(ra_type, levels = c('ra', 'spra', 'snra')),
					 sensspec = factor(sensspec, levels = c('sensitive', 'moderate', 'specific'))) |>
		glimpse() ->
		all_results

	# all_results |>
	# 	filter(term == 'has_chipTRUE') |>
	# 	# filter(cohort == 'ukb') |>
	# 	arrange(cohort, ra_type, term, desc(sensspec)) |>
	# 	distinct(cohort, ra_type, term, .keep_all = TRUE) |>
	# 	# select(ra_type, sensspec, n:cohort) |>
	# 	mutate(fu_time = time_at_risk_nochip + time_at_risk_chip,
	# 				 event_rate = 100 *n_event / (time_at_risk_nochip + time_at_risk_chip),
	# 				 event_rate_nochip = (n_event - n_events_w_chip)/time_at_risk_nochip,
	# 				 event_rate_chip = (n_events_w_chip)/time_at_risk_chip) |>
	# 	select(cohort, ra_type, sensspec, term, starts_with('event_rate'))
	# # select(cohort, ra_type, sensspec, term, fu_time) |>
	# print(n = Inf)
	# #starts_with('time_at'))


	# group_by(cohort, ra_type) |>
	# 	summarise(num_studied = length(n),
	# 						n_event = unique(n_event),
	# 						p_event = round(100 * n_event / n(), 1)) |>
	# 	print(n = Inf)

	all_results |> arrange(p.value) |> print(n = 20)

	# formal meta analyses
	fig1_data = list()
	for (this_ra_type in unique(all_results$ra_type)) {
		for (this_sensspec in unique(all_results$sensspec)) {
			for (this_term in unique(all_results$term)) {

				this_name = paste0(this_ra_type, '_', this_sensspec, ': ', this_term)
				message(this_name)

				bind_rows(
					filter(ukb, ra_type == this_ra_type, sensspec == this_sensspec, term == this_term),
					filter(aou, ra_type == this_ra_type, sensspec == this_sensspec, term == this_term),
					filter(biovu, ra_type == this_ra_type, sensspec == this_sensspec, term == this_term),
				) |>
					mutate(upper_ci = exp(estimate + 1.96*std.error)) ->
					tometa

				if (nrow(tometa) == 0) { next }
				this_metaanalysis = rma(yi = tometa$estimate,
																sei = tometa$std.error)
				tidy_result = tidy(this_metaanalysis)
				fig1_data[[this_name]] = mutate(tidy_result,
																				ra_type = this_ra_type,
																				sensspec = this_sensspec,
																				term = this_term) #,
				# n = sum(tometa$n),
				# n_chip = sum(tometa$n_chip),
				# n_event = sum(tometa$n_event))
			}
		}
	}

	bind_rows(fig1_data) |> pull(term) |> table()
	bind_rows(fig1_data) |> glimpse()

	bind_rows(fig1_data) |> filter(sensspec == 'specific') |> arrange(p.value)

	bind_rows(fig1_data) |>
		mutate(q.value = qvalue(p.value)$qvalues,
					 ra_type = factor(ra_type,
					 								 levels = c('ra', 'spra', 'snra'),
					 								 labels = c('RA', 'SPRA', 'SNRA')),
					 sensspec = factor(sensspec,
					 									levels = c('sensitive', 'moderate', 'specific')),
					 term = factor(str_sub(term, 5, -5),
					 							levels = c('chip', 'chip05', 'chip10', 'chip15',
					 												 'chip_d3a', 'chip_tet2', 'chip_asxl1'),
					 							labels = c('CHIP', 'CHIP 5%+', 'CHIP 10%+', 'CHIP 15%+',
					 												 'D3A CHIP', 'TET2 CHIP', 'ASXL1 CHIP')),
					 stars = case_when(
					 	p.value < 0.01 ~ '**',
					 	p.value < 0.05 ~ '*',
					 	.default = '')) ->
		toplot


	toplot |>
		filter(sensspec == 'specific') |>
		# print(n = Inf)
		# filter(!is.na(term)) |>
		ggplot(aes(y = exp(estimate),
							 x = term,
							 color = sensspec,
							 group = term,
							 ymin = exp(estimate - 1.96*std.error),
							 ymax = exp(estimate + 1.96*std.error))) +
		geom_hline(yintercept = 1) +
		geom_pointrange(
			# position = position_dodge2(width = 0.4)
		) +
		geom_text(mapping = aes(label = stars),
							position = position_nudge(x = 0.1, y = 0.1),
							size = 8,
							# position = position_dodge2(width = 0.7),
							show.legend = FALSE) +
		facet_grid(row = vars(ra_type))
	labs(title = 'Meta-analysis: Hazard ratio of CHIP for incident RA, SPRA, and SNRA',
			 color = 'disease\nascertainment',
			 x = NULL,
			 y = 'Hazard ratio and 95% confidence interval')

	ggsave(file.path(data_dir, 'meta_forest_plot_main_fig.pdf'), height = 6, width = 7)


	bind_rows(fig1_data) |>
		mutate(ra_type = factor(ra_type,
														levels = c('ra', 'spra', 'snra'),
														labels = c('RA', 'SPRA', 'SNRA')),
					 sensspec = factor(sensspec,
					 									levels = c('sensitive', 'moderate', 'specific')),
					 stars = case_when(p.value < 0.01 ~ '**',
					 									p.value < 0.05 ~ '*',
					 									.default = '')) |>
		ggplot(aes(y = exp(estimate),
							 x = term,
							 color = sensspec,
							 group = term,
							 ymin = exp(estimate - 1.96*std.error),
							 ymax = exp(estimate + 1.96*std.error))) +
		geom_hline(yintercept = 1) +
		geom_pointrange(position = position_dodge2(width = 0.4)) +
		geom_text(mapping = aes(label = stars),
							position = position_dodge2(width = 0.7),
							show.legend = FALSE) +
		facet_grid(row = vars(ra_type)) +
		labs(title = 'Meta-analysis: Hazard ratio of CHIP for incident RA, SPRA, and SNRA',
				 color = 'disease\nascertainment',
				 x = NULL,
				 y = 'Hazard ratio and 95% confidence interval')

	ggsave(file.path(data_dir, 'meta_forest_plot_supplementary_fig.pdf'), height = 6, width = 12)

}


# # nicer plotting
# bind_rows(ukb, aou, biovu) |>
# 	arrange(ra_type, term) |>
# 	mutate(idx = 1:n()) ->
# 	toplot
#
# toplot |>
# 	ggplot(mapping = aes(x = estimate,
# 											 y = interaction(term, cohort, sep = ' : '),
# 											 color = term,
# 											 xmin = exp(log(estimate) - 1.96*std.error),
# 											 xmax = exp(log(estimate) + 1.96*std.error))) +
# 	geom_pointrange() +
# 	facet_grid(row = vars(ra_type)) +
# 	scale_color_manual(values = c('black', 'blue')) +
# 	geom_vline(xintercept = 1) +
# 	theme_minimal() +
# 	theme(strip.background = element_rect(fill = 'gray'))
#
# ggsave('meta_forest_plot_supp_fig.pdf', height = 4, width = 7)
