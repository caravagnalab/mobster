

#' Title
#'
#' @param res
#' @param fig.lab
#' @param title
#' @param palette
#' @param TOP
#' @param cex
#' @param boxplot.ICL.range
#'
#' @return
#' @import ggpubr
#' @export
#'
#' @examples
mobster_plt_model_selection = function(res,
                               fig.lab = "",
                               title = 'MOBSTER best fit',
                               # palette = 'Set1',
                               TOP = 5,
                               cex = 1,
                               boxplot.ICL.range = NULL,
                               ...)
{
  # require(dbpmm)
  # require(ggplot2)
  # require(ggpubr)
  
  # Best model fit -- get plots
  best = plot.dbpmm(
    res$best,
    cex = 1 * cex,
    title = "",
    histogram.main = title,
    # palette = palette,
    silent = TRUE,
    annotation = paste0('Top-fit; N = ', nrow(res$best$data), ' mutations.'),
    ...
  )
  
  # Extract scores table
  scores = res$fits.table
  
  # Set format for columns with numbers
  scores.col = colnames(res$best$scores)
  scores[, scores.col] = apply(scores[, scores.col], 2, round, digit = 0)
  
  # Subset TOP scores and create a dataframe for plot
  if (nrow(scores) > TOP)
    scores = scores[1:TOP,]
  
  pio::pioTit("TOP scores")
  pio::pioDisp(scores)
  
  # Make a ggtexttable for `scores`
  df.scores = ggtexttable(
    scores,
    rows = NULL,
    theme = ttheme(
      base_size = 6 * cex,
      tbody.style = tbody_style(fill = rev(get_palette("RdBu", TOP)),
                                size = 6 * cex)
    )
  )
  
  for (s in scores.col) {
    idx = which.min(scores[, s]) + 1
    df.scores = table_cell_bg(
      df.scores,
      row = idx,
      column = which(colnames(scores) == s),
      linewidth = 2.5,
      fill = alpha("darkgreen", .5)
    )
  }
  
  # Model selection boxplot
  boxP = .plot.fit.summary(res,
                           alpha = .6,
                           silent = TRUE,
                           range = boxplot.ICL.range)$boxPlot
  
  # Other solutions ranked below top best -- maximum TOP 4 fixed
  TOP.plot = min(4, nrow(scores))
  
  # solutions = as.integer(rownames(scores)[2:(TOP.plot + 1)])
  # solutions = solutions[!is.na(solutions)]
  solutions = 
  
  # other.best = lapply(seq(solutions),
  other.best = lapply(2:(TOP.plot + 1),
                      function(w)
                        plot.dbpmm(
                          res$runs[[w]],
                          histogram.main = paste0("Solution #", w),
                          # annotation = paste0("Overall rank : ", solutions[w]),
                          silent = TRUE,
                          cex = .6 * cex
                        )$mainHist)
  
  # Final layout
  top_panel_best = ggarrange(
    best$mainHist,
    ggarrange(
      best$Initialization,
      best$Parameters,
      best$Proportions,
      ncol = 3,
      nrow = 1,
      labels = c("B", "C", "D")
    ),
    ggarrange(
      boxP,
      df.scores,
      ncol = 2,
      nrow = 1,
      labels = c("E", "F"),
      widths = c(1, 1.5)
    ),
    heights = c(2.5, 1),
    nrow = 3,
    ncol = 1,
    common.legend = TRUE,
    labels = 'A'
  )
  
  right_panel = ggarrange(
    plotlist = other.best,
    ncol = 1,
    nrow = length(other.best),
    labels = LETTERS[7:(7 + TOP.plot - 1)]
  )
  
  figure = ggarrange(
    top_panel_best,
    right_panel,
    widths = c(1.5, 1),
    nrow = 1,
    ncol = 2
  )
  
  figure = annotate_figure(
    figure,
    top = text_grob(
      "MOBSTER fit",
      color = "black",
      face = "bold",
      size = 18 * cex
    ),
    bottom = text_grob(
      " Panels: (A-D) Best fit; (E,F) Scores from model selection; (G-*) Lower-scoring fits.",
      color = "black",
      # face = "bold",
      hjust = 0,
      x = 0,
      size = 6 * cex
    ),
    fig.lab = fig.lab,
    fig.lab.face = "bold"
  )
  
  return(figure)
}
