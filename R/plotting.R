#' @title one cohort all results plot
#'
#' @param results the results to plot
#' @param title the title of the plot
#' @param p_or_q plot using p values or q values
#'
#' @return the plot
#' @export
#'
one_cohort_all_results_plot = function(results, title, p_or_q = 'p') {

	stopifnot(p_or_q %in% c('p', 'q'))
	p_or_q_name = paste0(p_or_q, ' value')

	results |>
		mutate(
			chip_type = factor(str_sub(term, 5, -5),
												 levels = c('chip', 'chip05', 'chip10', 'chip15',
												 					 'chip_d3a', 'chip_tet2', 'chip_asxl1'),
												 labels = c('CHIP', 'CHIP 5%+', 'CHIP 10%+', 'CHIP 15%+',
												 					 'DNMT3A CHIP', 'TET2 CHIP', 'ASXL1 CHIP')),
			sensspec = factor(sensspec, levels = c('sensitive', 'moderate', 'specific')),
			ra_type = factor(ra_type, levels = c('ra', 'spra', 'snra'),
											 labels = c('RA', 'SPRA', 'SNRA')),
			`p value` = factor(cut(p.value, breaks = c(0, 0.01, 0.05,  1)),
												 levels = c('(0,0.01]', '(0.01,0.05]', '(0.05,1]'),
												 labels = c('< 0.01', '< 0.05', 'n.s.')),
			`q value` = factor(cut(q.value, breaks = c(0, 0.01, 0.05,  1)),
												 levels = c('(0,0.01]', '(0.01,0.05]', '(0.05,1]'),
												 labels = c('< 0.01', '< 0.05', 'n.s.'))
		) |>
		ggplot(mapping = aes(
			x = exp(estimate),
			y = chip_type,
			xmin = exp(estimate - 1.96*std.error),
			xmax = exp(estimate + 1.96*std.error))
		) +
		geom_vline(xintercept = 1) +
		{ if (p_or_q == 'p') geom_pointrange(aes(color = `p value`)) } +
		{ if (p_or_q == 'q') geom_pointrange(aes(color = `q value`)) } +
		scale_color_manual(values = c('red', 'blue', 'black'), drop = FALSE) +
		theme_minimal() +
		labs(x = 'hazard ratio for RA',
				 y = NULL,
				 title = title) +
		facet_grid(rows = vars(ra_type, sensspec))
}

#' @title one cohort specific results plot
#'
#' @param results the results to plot
#' @param title the title of the plot
#' @param p_or_q plot using p values or q values
#'
#' @return the plot
#' @export
#'
one_cohort_specific_results_plot = function(results, title, p_or_q = 'p') {

	stopifnot(p_or_q %in% c('p', 'q'))
	p_or_q_name = paste0(p_or_q, ' value')

	results |>
		mutate(
			chip_type = factor(str_sub(term, 5, -5),
												 levels = c('chip', 'chip05', 'chip10', 'chip15'),
												 labels = c('CHIP', 'CHIP 5%+', 'CHIP 10%+', 'CHIP 15%+')),
			sensspec = factor(sensspec, levels = c('sensitive', 'moderate', 'specific')),
			ra_type = factor(ra_type, levels = c('ra', 'spra', 'snra'),
											 labels = c('RA', 'SPRA', 'SNRA')),
			`p value` = factor(cut(p.value, breaks = c(0, 0.01, 0.05,  1)),
												 levels = c('(0,0.01]', '(0.01,0.05]', '(0.05,1]'),
												 labels = c('< 0.01', '< 0.05', 'n.s.')),
			`q value` = factor(cut(q.value, breaks = c(0, 0.01, 0.05,  1)),
												 levels = c('(0,0.01]', '(0.01,0.05]', '(0.05,1]'),
												 labels = c('< 0.01', '< 0.05', 'n.s.'))
		) |>
		filter(!is.na(chip_type)) |>
		arrange(ra_type, chip_type, desc(sensspec)) |>
		distinct(ra_type, chip_type, .keep_all = TRUE) |>
		ggplot(mapping = aes(
			x = exp(estimate),
			y = chip_type,
			xmin = exp(estimate - 1.96*std.error),
			xmax = exp(estimate + 1.96*std.error))
		) +
		geom_vline(xintercept = 1) +
		{ if (p_or_q == 'p') geom_pointrange(aes(color = `p value`)) } +
		{ if (p_or_q == 'q') geom_pointrange(aes(color = `q value`)) } +
		scale_color_manual(values = c('red', 'blue', 'black'), drop = FALSE) +
		theme_minimal() +
		labs(x = 'hazard ratio for RA',
				 y = NULL,
				 title = title) +
		facet_grid(rows = vars(ra_type))
}
