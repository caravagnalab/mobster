#' Plot a MOBSTERh fit.
#'
#' @param x An object of class \code{"dbpmmh"}.
#' @param show_na Leave in the plot points for which a clustering assignment
#' is not available (\code{NA}).
#' @param add_density Add density of the fit on top of the histogram.
#' @param empty_plot Force the inclusion of an empty plot for a karyotype
#' that has not been used (from the list of supported karyotypes).
#' @param assembly_plot Turn the plot into a one-line strip, instead of returning
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
                       empty_plot = TRUE,
                       assembly_plot = TRUE,
                       ...)
{
  #############################################
  # Auxiliary function(s) private to the plot #
  #############################################

  # Generate dpareto density
  df_powerlaw_density = function(shape = 1,
                                 scale = 0.05,
                                 mixing = 0.5,
                                 location = Inf)
  {
    domain_x = seq(0, 1, 0.01)
    if (location == Inf)
      line_points = sads::dpareto(x = domain_x,
                                  shape = shape,
                                  scale = scale) * mixing
    else
      line_points = VGAM::dtruncpareto(
        x = domain_x,
        lower = scale,
        upper = location,
        shape = shape
      ) * mixing

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
  add_drivers = function(x, drivers_table, plot, facet = FALSE)
  {
    if ((drivers_table %>% nrow) == 0)
      return(plot)

    if (facet)
      plot = plot +
        scale_x_continuous(
          limits = c(0, 1.01),
          expand = c(0, 0),
          sec.axis = dup_axis(
            breaks = drivers_table$VAF,
            labels = drivers_table$driver_label,
            name = NULL
          )
        ) +
        geom_vline(
          data = drivers_table,
          aes(xintercept = VAF, color = cluster),
          size = .5,
          linetype = 'dashed',
          show.legend = FALSE
        ) +
        theme(axis.text.x.top = element_text(
          angle = 45,
          color = 'black',
          hjust = 0
        ))
    else
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

    return(plot)
  }

  # Has missing driver annotations?
  n_missing_drivers = function(x, drivers_table)
  {
    # Missing drivers
    x$data %>%
      filter(is_driver, is.na(cluster)) %>%
      nrow
  }

  #############################################
  # Auxiliary function(s) private to the plot #
  #############################################

  cli::cli_alert("Generating metadata for plot.")

  data_table = x$data

  if (!show_na)
    data_table = data_table %>% filter(!is.na(cluster))

  # Drivers table
  drivers_table = data_table %>%
    filter(is_driver, cluster != "Not used") %>%
    mutate(driver_label = ifelse(
      is.na(driver_label),
      paste(chr, from, paste0(ref, '>', alt), sep = ':'),
      driver_label
    ))

  # Inline model interpretation
  fit_interpreter = mobster:::clonality_interpreter(x) %>% filter(what == 'Subclone')

  fit_caption = "No subclonal expansions detected"

  if (fit_interpreter %>% nrow > 0)
    fit_caption = fit_interpreter %>%
    group_by(cluster) %>%
    summarise(l = paste(karyotype, collapse = ', ')) %>%
    mutate(l = paste0(cluster, ' [', l, ']')) %>%
    pull(l) %>%
    paste(collapse = ', ')

  fit_caption = paste0(fit_caption,
                       '; ',
                       n_missing_drivers(x, drivers_table),
                       ' driver(s) unassigned.')

  # Nonsense plot
  #
  # fit_caption = clonality_interpreter(x)
  # ggplot(fit_caption, aes(x = karyotype, y = cluster, fill = what)) +
  #   geom_tile() +
  #   CNAqc:::my_ggplot_theme()

  # Cluster colors
  cluster_colors = NULL

  nBeta = data_table %>%  filter(cluster != "Tail") %>%  pull() %>% unique() %>%  length()
  nBeta = nBeta[!is.na(nBeta)]
  #nBeta = nBeta[nBeta > 0] %>% max

  tail_color = 'gray'

  clonal_colors = suppressWarnings(RColorBrewer::brewer.pal(9, 'Set1'))[1:2]
  subclonal_colors = suppressWarnings(RColorBrewer::brewer.pal(7, 'Dark2'))[1:7]

  Beta_colors = c(clonal_colors, subclonal_colors)
  names(Beta_colors) = c("C1", "C2", "S1", "S2", "S3", "S4", "S5", "S6", "S7")

  cluster_colors = c("Tail" = tail_color, Beta_colors, `Not used` = 'lightpink')

  cluster_colors = cluster_colors[data_table$cluster %>% unique]

  # VAF plot, make a temporary plot to return if not densities are required
  VAF_binwidth = 0.01

  nkaryo = data_table %>%
    group_by(karyotype) %>%
    summarise(n = n()) %>%
    mutate(label_karyotype = paste0(karyotype, " (n = ", n, ')'))

  nkaryo_labels = nkaryo$label_karyotype
  names(nkaryo_labels) = nkaryo$karyotype

  # data_table$karyotype = nkaryo_labels[data_table$karyotype]


  density_plot = ggplot(data_table %>% mutate(karyotype = nkaryo_labels[karyotype]),
                        aes(x = VAF)) +
    geom_histogram(
      aes(y = ..count.. / sum(..count..), fill = cluster %>% paste),
      binwidth = 0.01,
      alpha = 0.6
    ) +
    facet_wrap( ~ karyotype, scales = 'free_y') +
    CNAqc:::my_ggplot_theme() +
    scale_fill_manual(values = cluster_colors) +
    scale_color_manual(values = cluster_colors) +
    guides(fill = guide_legend("Cluster")) +
    labs(y = "Density", caption = fit_caption)

  if (!add_density)
    return(add_drivers(x, drivers_table, density_plot, facet = FALSE))

  cli::cli_alert("Generating fit densities.")

  # Used karyotypes
  used_karyotypes = c("1:0", "1:1", "2:0", "2:1", "2:2")

  if (mobster:::has_tail(x)) {
    pareto_params = mobster:::get_pareto(x)
    # Power law density per karyotype

    pareto_params_df = NULL
    for (i in 1:nrow(pareto_params))
      pareto_params_df = pareto_params_df %>%
      bind_rows(
        df_powerlaw_density(
          shape = pareto_params$shape[i],
          scale = pareto_params$scale[i],
          mixing = pareto_params$mixing[i],
          location = pareto_params$location[i]
        ) %>%
          mutate(karyotype = pareto_params$karyotype[i],
                 cluster = "Tail")
      )

    pareto_params_df_low = NULL
    for (i in 1:nrow(pareto_params))
      pareto_params_df_low = pareto_params_df_low %>%
      bind_rows(
        df_powerlaw_density(
          shape = qlnorm(
            0.05,
            meanlog = log(pareto_params$shape[i]),
            sdlog = pareto_params$shape_noise[i]
          ),
          scale = pareto_params$scale[i],
          mixing = pareto_params$mixing[i],
          location = pareto_params$location[i]
        ) %>%
          mutate(karyotype = pareto_params$karyotype[i],
                 cluster = "Tail")
      )

    pareto_params_df_high = NULL
    for (i in 1:nrow(pareto_params))
      pareto_params_df_high = pareto_params_df_high %>%
      bind_rows(
        df_powerlaw_density(
          shape = qlnorm(
            0.95,
            meanlog = log(pareto_params$shape[i]),
            sdlog = pareto_params$shape_noise[i]
          ),
          scale = pareto_params$scale[i],
          mixing = pareto_params$mixing[i],
          location = pareto_params$location[i]
        ) %>%
          mutate(karyotype = pareto_params$karyotype[i],
                 cluster = "Tail")
      )
  }


  # Beta density per karyotype
  Beta_params = mobster:::get_beta(x)

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
  fit_plots = lapply(used_karyotypes, function(x) {
    ggplot(data = data.frame(x = 0, y = 0, label = "X"), aes(x = x,
                                                             y = y, label = label)) + CNAqc:::my_ggplot_theme() +
      theme(
        panel.border = element_rect(
          colour = "black",
          fill = NA,
          linetype = "dashed"
        ),
        panel.background = element_rect(fill = "gainsboro"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()
      )
  })
  names(fit_plots) = used_karyotypes

  used_karyotypes_plot = x$model_parameters %>% names %>% seq_along()

  for (s_k in used_karyotypes_plot)
  {
    k = (x$model_parameters %>% names)[s_k]

    cli::cli_alert("Generating plot for {.field {k}}")

    density_plot = ggplot(
      data_table %>% filter(karyotype == k) %>% mutate(karyotype = nkaryo_labels[karyotype]),
      aes(VAF)
    ) +
      geom_histogram(
        aes(y = ..count.. / sum(..count..), fill = cluster %>% paste),
        binwidth = 0.01,
        alpha = 0.5
      ) +
      facet_wrap( ~ karyotype, scales = 'free_y') +
      CNAqc:::my_ggplot_theme() +
      scale_fill_manual(values = cluster_colors) +
      guides(fill = guide_legend("Cluster")) +

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

    # Add tail density
    if (mobster:::has_tail(x)) {
      density_plot = density_plot + geom_line(
        data = pareto_params_df %>% filter(karyotype == k),
        aes(
          x = x,
          y = y * VAF_binwidth,
          color = cluster
        ),
        size = 1,
        inherit.aes = FALSE,
        show.legend = FALSE
      ) +  geom_line(
        data = pareto_params_df_low %>% filter(karyotype == k),
        aes(x = x,
            y = y * VAF_binwidth,),
        size = 0.5,
        color = "black",
        linetype = "dashed",
        inherit.aes = FALSE,
        show.legend = FALSE
      ) +  geom_line(
        data = pareto_params_df_high %>% filter(karyotype == k),
        aes(
          x = x,
          y = y * VAF_binwidth,
          color = cluster
        ),
        size = 0.5,
        color = "black",
        linetype = "dashed",
        inherit.aes = FALSE,
        show.legend = FALSE
      )
    }

    density_plot = add_drivers(x,
                               drivers_table %>% filter(karyotype == k),
                               density_plot,
                               facet = TRUE)

    if (s_k == 1)
      density_plot = density_plot + labs(y = "Density")
    else
      density_plot = density_plot + labs(y = NULL)

    if (s_k == (used_karyotypes_plot %>% length)) {
      density_plot = density_plot +
        labs(caption = fit_caption)
    }

    fit_plots[[k]] = density_plot
  }

  # Remove empty plots if required
  if (!empty_plot) {
    fit_plots = fit_plots[x$model_parameters %>% names]
  }

  if (!assembly_plot)
    return(fit_plots)

  cowplot::plot_grid(
    plotlist = fit_plots,
    nrow = 1,
    align = 'h',
    axis = 'tb'
  )
}
