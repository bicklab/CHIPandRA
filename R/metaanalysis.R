# library(tidyverse)
# library(metafor)
# library(broom)
# library(glue)
# library(patchwork)
#
#
# data_dir = '~/OneDrive - VUMC/CHIP and RA'
#
# terms = c('CHIP', 'CHIP 10%+')
# ra_types = c('RA', 'SPRA', 'SNRA')
#
# read_csv(file.path(data_dir, 'ukb_chip_ra_results_toplot.csv')) |>
# 	select(-ra_def) |>
# 	filter(def_type == 'w meds & labs') |>
# 	filter(term %in% terms) |>
# 	mutate(def_type = 'spec.') |>
# 	mutate(cohort = 'ukb') |>
# 	glimpse() ->
# 	ukb
#
# read_csv(file.path(data_dir, 'aou_chip_ra_results_toplot.csv')) |>
# 	select(-ra_def) |>
# 	filter(def_type == 'w meds & labs') |>
# 	filter(term %in% terms) |>
# 	mutate(def_type = 'spec.') |>
# 	mutate(cohort = 'aou') |>
# 	glimpse() ->
# 	aou
#
# read_csv(file.path(data_dir, 'biovu_chip_ra_effects_toplot.csv')) |>
# 	select(term, estimate, std.error, p.value, n, n_event, n_chip,
# 			 ra_type, def_type, stars) |>
# 	filter(def_type == 'spec.') |>
# 	filter(term %in% terms) |>
# 	mutate(cohort = 'biovu') |>
# 	glimpse() ->
# 	biovu
#
#
# # formal meta analysis ( x 6 for each type of RA and CHIP vs big_chip)
# fig1_data = list()
# for (this_ra_type in ra_types) {
# 	for (this_term in terms) {
#
# 		this_name = paste0(this_ra_type, ': ', this_term)
# 		message(this_name)
#
# 		bind_rows(
# 			filter(ukb, ra_type == this_ra_type, term == this_term),
# 			filter(aou, ra_type == this_ra_type, term == this_term),
# 			filter(biovu, ra_type == this_ra_type, term == this_term),
# 		) |>
# 			mutate(upper_ci = exp(log(estimate) + 2*std.error)) ->
# 			tometa
#
# 		this_metaanalysis = rma(yi = log(tometa$estimate),
# 														sei = tometa$std.error)
# 		tidy_result = tidy(this_metaanalysis)
# 		fig1_data[[this_name]] = mutate(tidy_result, ra_type = this_ra_type, term = this_term,
# 																		n = sum(tometa$n),
# 																		n_chip = sum(tometa$n_chip),
# 																		n_event = sum(tometa$n_event))
# 		print(tidy_result)
#
# 		this_term_sanitized = str_remove_all(make.names(this_term), '\\.')
# 		xpos = switch(this_ra_type,
# 									'RA' = -1.5,
# 									'SPRA' = -3.5,
# 									'SNRA' = -3)
# 		r_col_pos = max(tometa$upper_ci)
# 		pval = tidy_result$p.value
#
# 		pdf(file = glue('meta_forest_plot_{this_ra_type}_{this_term_sanitized}.pdf'),
# 				width = 7, height = 4)
# 		forest(this_metaanalysis,
# 					 transf = exp,
# 					 refline = 1,
# 					 main = glue('Meta-analysis of Effect of {this_term} on {this_ra_type}'),
# 					 slab = c('UKB', 'AOU', 'BioVU'),
# 					 shade = 'zebra',
# 					 header = c('Cohort', 'Hazard Ratio [95% CI]'))
# 		dev.off()
# 	}
# }
#
# bind_rows(fig1_data) |>
# 	mutate(ra_type = factor(ra_type, levels = c('RA', 'SPRA', 'SNRA'))) |>
# 	ggplot(aes(x = exp(estimate),
# 						 y = term,
# 						 color = term,
# 						 xmin = exp(estimate - 1.96*std.error),
# 						 xmax = exp(estimate + 1.96*std.error))) +
# 	geom_vline(xintercept = 1) +
# 	geom_pointrange() +
# 	geom_text(aes(label = paste0('p = ', round(p.value, 3))),
# 						position = position_nudge(x = 0.1, y = 0.2),
# 						hjust = 'left',
# 						show.legend = FALSE,
# 						size = 3) +
# 	geom_text(aes(label = paste0(format(n_chip, big.mark = ','), ' with CHIP, ', format(n_event, big.mark = ','), ' cases')),
# 						position = position_nudge(x = 0.1, y = -0.2),
# 						hjust = 'left',
# 						show.legend = FALSE,
# 						size = 3) +
# 	facet_grid(row = vars(ra_type)) +
# 	scale_color_manual(values = c('black', 'blue')) +
# 	theme_minimal() +
# 	theme(strip.background = element_rect(fill = 'gray'),
# 				axis.title.y = element_blank(),
# 				legend.title = element_blank(),
# 				legend.background = element_rect(),
# 				legend.position = c(0.8,0.8)) +
# 	labs(title = 'Meta-analysis: Hazard ratio of CHIP for incident RA, SPRA, and SNRA',
# 			 x = 'Hazard ratio and 95% confidence interval')
#
# ggsave('meta_forest_plot_main_fig.pdf', height = 4, width = 7)
#
#
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
