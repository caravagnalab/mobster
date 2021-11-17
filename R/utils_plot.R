# Return smar colours
smart_colors = function(x, pl, colors)
{
  if(all(class(x) == 'dbpmm'))
  {
    # clusters in x
    wh_col = unique(x$data$cluster)
  }
  else wh_col = x
  
  # Missing colors
  wh_col_missing = !(wh_col %in% names(colors))
  wh_col = wh_col[wh_col_missing]
  
  # Complement colors
  mycolors = colors
  new_col = NULL
  if(length(wh_col) < 9) 
    new_col = suppressWarnings(RColorBrewer::brewer.pal(length(wh_col), 'Set1'))
  else
    new_col = rainbow(length(wh_col))
  
  names(new_col) = sort(wh_col)

  return(c(mycolors, new_col[!is.na(names(new_col))]))
}

# add colours to a ggplot plot
add_color_pl = function(x, pl, colors)
{
  pl + scale_color_manual(values = smart_colors(x, pl, colors))
}

# add fill to a ggplot plot
add_fill_pl = function(x, pl, colors)
{
  # if(!is.vector(colors) | any(is.na(colors))) return(pl)
  # 
  # # clusters in x
  # wh_col = unique(x$data$cluster)
  # stopifnot(all(wh_col %in% names(colors)))
  
  pl + scale_color_manual(values = smart_colors(x, pl, colors))
}

# add colours and fill to a ggplot plot
add_fill_color_pl = function(x, pl, colors)
{
  # if(!is.vector(colors) | any(is.na(colors))) return(pl)
  # 
  # # clusters in x
  # wh_col = unique(x$data$cluster)
  # stopifnot(all(wh_col %in% names(colors)))
  
  pl + 
    scale_color_manual(values = smart_colors(x, pl, colors)) +
    scale_fill_manual(values = smart_colors(x, pl, colors))
}

# Default ggplot theme for all plots of this package
my_ggplot_theme = function(cex = 1) 
{
  cex_opt = getOption('mobster_cex', default = 1)
  
  theme_light(base_size = 10 * cex_opt) +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(.3 * cex_opt, "cm"),
      panel.background = element_rect(fill = 'white')
    )
}