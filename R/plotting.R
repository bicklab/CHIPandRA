#' @title one cohort results plot
#'
#' @param results the results to plot
#' @param title the title of the plot
#' @param andor whether to plot and-based definition or or-based definition of disease
#' @param p_or_q plot using p values or q values
#'
#' @return the plot
#' @export
#'
one_cohort_results_plot = function(results,
																	 title,
																	 andor = 'and',
																	 p_or_q = 'p') {

	results |>
		filter(sensspec == andor) |>
		separate_wider_delim(term, delim = '_', names = c(NA, 'gene', 'size')) |>
		mutate(ra_type = factor(ra_type, levels = c('ra', 'spra', 'snra'), labels = c('RA', 'SPRA', 'SNRA'))) |>
		mutate(gene = as.numeric(factor(gene, levels = c('chip', 'dnmt3a', 'tet2', 'asxl1')))) |>
		mutate(size = 0.2*as.numeric(factor(size, levels = c('small', 'big'))) - 0.3) |>
		mutate(stars = case_when(p_or_q == 'p' & p.value < 0.001 ~ ', **',
														 p_or_q == 'p' & p.value < 0.01 ~ ', *',
														 p_or_q == 'p' ~ '',
														 p_or_q == 'q' & q.value < 0.001 ~ ', **',
														 p_or_q == 'q' & q.value < 0.01 ~ ', *',
														 p_or_q == 'q' ~ '',
														 .default = 'ERROR')) |>
		ggplot(mapping = aes(
			x = gene + size,
			y = exp(estimate),
			ymin = exp(estimate - 1.96*std.error),
			ymax = exp(estimate + 1.96*std.error))
		) +
		geom_hline(yintercept = 1, color = 'blue') +
		geom_pointrange() +
		geom_text(mapping = aes(label = paste0(round(exp(estimate), 1),  stars)), position = position_nudge(x = 0.05, y = 0.3)) +
		facet_grid(rows = vars(ra_type)) +
		scale_x_continuous(
			breaks = 1:4,
			labels = c('any CHIP', 'DNMT3A', 'TET2', 'ASXL1')) +
		theme_bw() +
		xlab('CHIP type') +
		ylab('Hazard ratio') +
		ggtitle('Hazard Ratio of CHIP for Incident RA')
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
		filter(sensspec == 'and') |>
		mutate(
			chip_type = factor(str_sub(term, 5, -5),
												 levels = c('chip', 'chip05', 'chip10', 'chip15',
												 					 'chip_d3a', 'chip_tet2', 'chip_asxl1',
												 					 'chip_d3a_10', 'chip_tet2_10', 'chip_asxl1_10'),
												 labels = c('CHIP', 'CHIP 5%+', 'CHIP 10%+', 'CHIP 15%+',
												 					 'DNMT3A CHIP', 'TET2 CHIP', 'ASXL1 CHIP',
												 					 'DNMT3A CHIP 10%+', 'TET2 CHIP 10%+', 'ASXL1 CHIP 10%+')),
			sensspec = factor(sensspec,
												levels = c('or', 'and'),
												labels = c('OR', 'AND')),
			ra_type = factor(ra_type, levels = c('ra', 'spra', 'snra'),
											 labels = c('RA', 'SPRA', 'SNRA')),
			`p value` = factor(cut(p.value, breaks = c(0, 0.01, 0.05,  1)),
												 levels = c('(0,0.01]', '(0.01,0.05]', '(0.05,1]'),
												 labels = c('< 0.01', '< 0.05', 'n.s.'))
			# `q value` = factor(cut(q.value, breaks = c(0, 0.01, 0.05,  1)),
			# 									 levels = c('(0,0.01]', '(0.01,0.05]', '(0.05,1]'),
			# 									 labels = c('< 0.01', '< 0.05', 'n.s.'))
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
