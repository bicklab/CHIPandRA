# library(tidyverse)
# library(metafor)
# library(broom)
# library(glue)
#
# data_dir = '~/OneDrive - VUMC/CHIP and RA'
#
# # terms = c('CHIP', 'CHIP 10%+')
# # ra_types = c('RA', 'SPRA', 'SNRA')
#
# read_csv(file.path(data_dir, 'ukb_results/ukb_chip_ra_results.csv')) |>
# 	mutate(cohort = 'ukb') |>
# 	glimpse() ->
# 	ukb
#
# ukb |> one_cohort_all_results_plot(title = 'UKB: hazard ratio for incident RA', p_or_q = 'p')
# ggsave(file.path(data_dir, 'ukb_results', 'ukb_chip_ra_forest_plot_pvals.pdf'),
# 			 height = 9, width = 9)
#
# ukb |> one_cohort_all_results_plot(title = 'UKB: hazard ratio for incident RA', p_or_q = 'q')
# ggsave(file.path(data_dir, 'ukb_results', 'ukb_chip_ra_forest_plot_qvals.pdf'),
# 			 height = 9, width = 9)
#
# ukb |> one_cohort_specific_results_plot(title = 'UKB: hazard ratio for incident RA')
# ggsave(file.path(data_dir, 'ukb_results', 'ukb_chip_ra_focused_forest_plot_pvals.pdf'),
# 			 height = 3, width = 9)
#
#
# read_csv(file.path(data_dir, 'aou_results/aou_chip_ra_results.csv')) |>
# 	mutate(cohort = 'aou') |>
# 	glimpse() ->
# 	aou
#
# aou |> one_cohort_all_results_plot(title = 'AOU: hazard ratio for incident RA')
# aou |> one_cohort_specific_results_plot(title = 'AOU: hazard ratio for incident RA')
#
#
# read_csv(file.path(data_dir, 'biovu_results/biovu_chip_ra_results.csv')) |>
# 	mutate(cohort = 'biovu') |>
# 	glimpse() ->
# 	biovu
#
# bind_rows(ukb, aou, biovu) |>
# 	mutate(cohort = factor(cohort, levels = c('ukb', 'aou', 'biovu')),
# 				 ra_type = factor(ra_type, levels = c('ra', 'spra', 'snra')),
# 				 sensspec = factor(sensspec, levels = c('sensitive', 'moderate', 'specific'))) |>
# 	glimpse() ->
# 	all_results
#
# all_results |>
# 	group_by(cohort, ra_type, sensspec) |>
# 	summarise(n = unique(n),
# 						n_event = unique(n_event)) |>
# 	print(n = Inf)
#
#
# all_results |> arrange(p.value) |> print(n = 20)
#
# # formal meta analyses
# fig1_data = list()
# for (this_ra_type in unique(all_results$ra_type)) {
# 	for (this_sensspec in unique(all_results$sensspec)) {
# 		for (this_term in unique(all_results$term)[1:4]) {
#
# 			this_name = paste0(this_ra_type, '_', this_sensspec, ': ', this_term)
# 			message(this_name)
#
# 			bind_rows(
# 				filter(ukb, ra_type == this_ra_type, sensspec == this_sensspec, term == this_term),
# 				filter(aou, ra_type == this_ra_type, sensspec == this_sensspec, term == this_term),
# 				filter(biovu, ra_type == this_ra_type, sensspec == this_sensspec, term == this_term),
# 			) |>
# 				mutate(upper_ci = exp(estimate + 1.96*std.error)) ->
# 				tometa
#
# 			this_metaanalysis = rma(yi = tometa$estimate,
# 															sei = tometa$std.error)
# 			tidy_result = tidy(this_metaanalysis)
# 			fig1_data[[this_name]] = mutate(tidy_result,
# 																			ra_type = this_ra_type,
# 																			sensspec = this_sensspec,
# 																			term = this_term,
# 																			n = sum(tometa$n),
# 																			n_chip = sum(tometa$n_chip),
# 																			n_event = sum(tometa$n_event))
# 			print(tidy_result)
#
# 			# this_term_sanitized = str_remove_all(make.names(this_term), '\\.')
# 			# xpos = switch(this_ra_type,
# 			# 							'RA' = -1.5,
# 			# 							'SPRA' = -3.5,
# 			# 							'SNRA' = -3)
# 			# r_col_pos = max(tometa$upper_ci)
# 			# pval = tidy_result$p.value
# 			#
# 			# pdf(file = glue('meta_forest_plot_{this_ra_type}_{this_term_sanitized}.pdf'),
# 			# 		width = 7, height = 4)
# 			# forest(this_metaanalysis,
# 			# 			 transf = exp,
# 			# 			 refline = 1,
# 			# 			 main = glue('Meta-analysis of Effect of {this_term} on {this_ra_type}'),
# 			# 			 slab = c('UKB', 'AOU', 'BioVU'),
# 			# 			 shade = 'zebra',
# 			# 			 header = c('Cohort', 'Hazard Ratio [95% CI]'))
# 			# dev.off()
# 		}
# 	}
# }
#
# bind_rows(fig1_data) |>
# 	glimpse()
#
# bind_rows(fig1_data) |>
# 	mutate(ra_type = factor(ra_type,
# 													levels = c('ra', 'spra', 'snra'),
# 													labels = c('RA', 'SPRA', 'SNRA')),
# 				 sensspec = factor(sensspec,
# 				 									levels = c('sensitive', 'moderate', 'specific')),
# 				 term = factor(str_sub(term, 5, -5),
# 				 							levels = c('chip', 'chip05', 'chip10', 'chip15'),
# 				 							labels = c('CHIP', 'CHIP 5%+', 'CHIP 10%+', 'CHIP 15%+')),
# 				 stars = case_when(p.value < 0.01 ~ '**',
# 				 									p.value < 0.05 ~ '*',
# 				 									.default = '')) |>
# 	ggplot(aes(y = exp(estimate),
# 						 x = term,
# 						 color = sensspec,
# 						 group = term,
# 						 ymin = exp(estimate - 1.96*std.error),
# 						 ymax = exp(estimate + 1.96*std.error))) +
# 	geom_hline(yintercept = 1) +
# 	geom_pointrange(position = position_dodge2(width = 0.4)) +
# 	geom_text(mapping = aes(label = stars),
# 						position = position_dodge2(width = 0.7),
# 						show.legend = FALSE) +
# 	# geom_text(aes(label = paste0('p = ', round(p.value, 3))),
# 	# 					position = position_nudge(x = 0.1, y = 0.2),
# 	# 					hjust = 'left',
# 	# 					show.legend = FALSE,
# 	# 					size = 3) +
# 	# geom_text(aes(label = paste0(format(n_chip, big.mark = ','), ' with CHIP, ', format(n_event, big.mark = ','), ' cases')),
# 	# 					position = position_nudge(x = 0.1, y = -0.2),
# 	# 					hjust = 'left',
# 	# 					show.legend = FALSE,
# 	# 					size = 3) +
# 	facet_grid(row = vars(ra_type)) +
# # scale_color_manual(values = c('black', 'blue')) +
# # theme_minimal() +
# # theme(strip.background = element_rect(fill = 'gray'),
# # 			axis.title.y = element_blank(),
# # 			legend.title = element_blank(),
# # 			legend.background = element_rect(),
# # 			legend.position = c(0.8,0.8)) +
# 	labs(title = 'Meta-analysis: Hazard ratio of CHIP for incident RA, SPRA, and SNRA',
# 			 color = 'disease\nascertainment',
# 			 x = NULL,
# 			 y = 'Hazard ratio and 95% confidence interval')
#
# ggsave('meta_forest_plot_main_fig.pdf', height = 6, width = 7)
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
