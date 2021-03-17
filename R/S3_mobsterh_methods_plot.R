#' Plot a MOBSTERh fit.
#'
#' @param x An object of class \code{"dbpmmh"}.
#' @param show_na Leave in the plot points for which a clustering assignment
#' is not available (\code{NA}).
#' @param add_density Add density of the fit on top of the histogram.
#' @param ... Unused.
#'
#' @return A ggplot or a cowplot object for the plot (depends on parameters).
#'
#'
#' @exportS3Method plot dbpmmh
#' @export plot.dbpmmh
#'
#' @examples
plot.dbpmmh = function(x,
                       show_na = FALSE,
                       add_density = TRUE,
                       ...)
{
  #############################################
  # Auxiliary function(s) private to the plot #
  #############################################

  # Generate dpareto density
  df_powerlaw_density = function(shape = 1,
                                 scale = 0.05,
                                 mixing = 0.5)
  {
    domain_x = seq(0, 1, 0.01)
    line_points = sads::dpareto(x = domain_x,
                                shape = shape,
                                scale = scale) * mixing

    data.frame(x = domain_x, y = line_points)
  }

  # Generate dbeta density
  df_Beta_density = function(a = 10,
                             b = 10,
                             mixing = 0.5)
  {
    domain_x = seq(0, 1, 0.01)
    line_points = dbeta(x = domain_x,
                        shape1 = as.numeric(a),
                        shape2 = as.numeric(b)) * mixing

    data.frame(x = domain_x, y = line_points)
  }

  # Add drivers annotation
  add_drivers = function(x, drivers_table, plot)
  {
    # Missing drivers
    ndrivers_missing = x$table %>%
      filter(is_driver, is.na(cluster)) %>%
      nrow

    if (drivers_table %>% nrow > 0)
      plot = plot +
        geom_vline(
          data = drivers_table,
          aes(xintercept = VAF, color = cluster),
          size = .5,
          linetype = 'dashed',
          show.legend = FALSE
        ) +
        geom_text(
          data = drivers_table,
          aes(
            x = VAF - 0.02,
            y = Inf,
            color = cluster,
            label = driver_label
          ),
          angle = 90,
          size = 3,
          show.legend = FALSE,
          hjust = 1
        )

    if (ndrivers_missing > 0)
      plot = plot + labs(subtitle = paste0(ndrivers_missing, " driver annotated has not been analysed"))

    return(plot)
  }

  #############################################
  # Auxiliary function(s) private to the plot #
  #############################################

  cli::cli_alert("Generating metadata for plot.")

  data_table = x$table

  if (!show_na)
    data_table = data_table %>% filter(!is.na(cluster))

  # Inline model interpretation
  fit_interpreter = clonality_interpreter(x) %>% filter(what == 'Subclone')

  fit_caption = "No subclonal expansions detected"
  if (fit_interpreter %>% nrow > 0)
    fit_caption = fit_interpreter %>%
    group_by(cluster) %>%
    summarise(l = paste(karyotype, collapse = ', ')) %>%
    mutate(l = paste0(cluster, ' [', l, ']')) %>%
    pull(l) %>%
    paste(collapse = ', ')

  # Nonsense plot
  #
  # fit_caption = clonality_interpreter(x)
  # ggplot(fit_caption, aes(x = karyotype, y = cluster, fill = what)) +
  #   geom_tile() +
  #   CNAqc:::my_ggplot_theme()

  # Cluster colors
  cluster_colors = NULL

  nBeta = data_table$cluster %>% unique()
  nBeta = nBeta[!is.na(nBeta)]
  nBeta = nBeta[nBeta > 0] %>% max

  tail_color = 'gray'

  Beta_colors = suppressWarnings(RColorBrewer::brewer.pal(nBeta, 'Set1'))
  Beta_colors = Beta_colors[1:nBeta]

  cluster_colors = c(tail_color, Beta_colors, `Not used` = 'lightpink')
  names(cluster_colors)[1:(length(cluster_colors) - 1)] = c("Tail", paste0("Beta", 1:nBeta))

  # Transform cluster labels into proper intelligle names
  data_table = data_table %>%
    mutate(cluster = case_when(
      cluster == 0 ~ "Tail",
      is.na(cluster) ~ "Not used",
      TRUE ~ paste0("Beta", cluster)
    ))

  # Drivers table
  drivers_table = data_table %>%
    filter(is_driver, cluster != "Not used") %>%
    mutate(driver_label = ifelse(
      is.na(driver_label),
      paste(chr, from, paste0(ref, '>', alt), sep = ':'),
      driver_label
    ))

  # VAF plot, make a temporary plot to return if not densities are required
  VAF_binwidth = 0.01

  density_plot = ggplot(data_table,
                        aes(VAF)) +
    geom_histogram(
      aes(y = ..count.. / sum(..count..), fill = cluster %>% paste),
      binwidth = 0.01,
      alpha = 0.6
    ) +
    facet_wrap(~ karyotype, scales = 'free_y') +
    CNAqc:::my_ggplot_theme() +
    scale_fill_manual(values = cluster_colors) +
    scale_color_manual(values = cluster_colors) +
    guides(fill = guide_legend("Cluster")) +
    labs(y = "Density", caption = fit_caption)

  if (!add_density)
    return(add_drivers(x, drivers_table, density_plot))

  cli::cli_alert("Generating fit densities.")

  # Used karyotypes
  used_karyotypes = c("1:0", "1:1", "2:0", "2:1", "2:2")

  # Power law density per karyotype
  pareto_params = get_pareto(x)

  pareto_params_df = NULL
  for (i in 1:nrow(pareto_params))
    pareto_params_df = pareto_params_df %>%
    bind_rows(
      df_powerlaw_density(
        shape = pareto_params$shape[i],
        scale = pareto_params$scale[i],
        mixing = pareto_params$mixing[i]
      ) %>%
        mutate(karyotype = pareto_params$karyotype[i],
               cluster = "Tail")
    )

  # Beta density per karyotype
  Beta_params = get_beta(x)

  beta_params_df = NULL
  for (i in 1:nrow(Beta_params))
    beta_params_df = beta_params_df %>%
    bind_rows(
      df_Beta_density(
        a = Beta_params$a[i],
        b = Beta_params$b[i],
        mixing = Beta_params$mixing[i]
      ) %>%
        mutate(
          karyotype = Beta_params$karyotype[i],
          cluster = Beta_params$cluster[i]
        )
    )

  # Create one plot per karyotype

  fit_plots = NULL
  used_karyotypes_plot = x$model_parameters %>% names %>% seq_along()

  for (s_k in used_karyotypes_plot)
  {
    k = (x$model_parameters %>% names)[s_k]

    cli::cli_alert("Generating and assembly plot for {.field {k}}")

    density_plot = ggplot(data_table %>% filter(karyotype == k),
                          aes(VAF)) +
      geom_histogram(
        aes(y = ..count.. / sum(..count..), fill = cluster %>% paste),
        binwidth = 0.01,
        alpha = 0.5
      ) +
      facet_wrap(~ karyotype, scales = 'free_y') +
      CNAqc:::my_ggplot_theme() +
      scale_fill_manual(values = cluster_colors) +
      guides(fill = guide_legend("Cluster")) +
      geom_line(
        data = pareto_params_df %>% filter(karyotype == k),
        aes(
          x = x,
          y = y * VAF_binwidth,
          color = cluster
        ),
        size = 1,
        inherit.aes = FALSE,
        show.legend = FALSE
      ) +
      scale_color_manual(values = cluster_colors) +
      geom_line(
        data = beta_params_df %>% filter(karyotype == k),
        aes(
          x = x,
          y = y * VAF_binwidth,
          color = cluster
        ),
        inherit.aes = FALSE,
        show.legend = FALSE,
        size = 1
      )

    density_plot = add_drivers(x, drivers_table %>% filter(karyotype == k), density_plot)

    if (s_k == 1)
      density_plot = density_plot + labs(y = "Density")
    else
      density_plot = density_plot + labs(y = NULL)

    if (s_k == used_karyotypes_plot %>% length) {
      density_plot = density_plot +
        labs(caption = fit_caption)
    }

    fit_plots = append(fit_plots, density_plot %>% list)
  }

  cowplot::plot_grid(plotlist = fit_plots,
                     align = 'h',
                     axis = 'tb')
}
