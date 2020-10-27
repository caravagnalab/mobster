#' Generate a random MOBSTER model and data.
#'
#' @description Generate a random MOBSTER model, its data and creates a plot for it.
#'
#' @param N Number of samples to generate (mutations).
#' @param K_betas Number of Beta components (subclones).
#' @param pi_tail_bounds 2D vector with min and max size of the tail's mutations (proportions).
#' @param pi_min Minimum mixing proportion for every component.
#' @param Betas_separation Minimum separation between the means of the Beta components.
#' @param Beta_variance_scaling The variance of the Beta is generated as U[0,1] and scaled by this value.
#' Values on the order of 1000 give low variance, 100 represents a dataset with quite some dispersion (
#' compared to a putative Binomial generative model).
#' @param Beta_bounds Range of values to sample the Beta means.
#' @param shape_bounds Range of values to sample the tail shape, default [1, 3],
#' @param scale Tail scale, default 0.05.
#' @param seed The seed to fix the process, default is 123.
#'
#' @return A list with the dataset in a tibble, the model parameters and a plot the data.
#'
#' @export
#'
#' @examples
#' x = random_dataset()
#' print(x)
random_dataset = function(N = 5000,
                          K_betas = 2,
                          pi_tail_bounds = c(.2, .4),
                          pi_min = 0.1,
                          Betas_separation = 0.1,
                          Beta_variance_scaling = 1e3,
                          Beta_bounds = c(.1, .9),
                          shape_bounds = c(1, 1, 3),
                          scale = 0.05,
                          seed = NULL)
{
  # require(MCMCpack)
  
  set.seed(seed)
  
  Bt = paste0('C', 1:K_betas)
  
  # Pi's
  repeat {
    # Avoid to import MCMCpack just for this..
    # pi = as.vector(MCMCpack::rdirichlet(1, rep(1, K_betas + 1)))
    pi = runif(K_betas + 1)
    pi = pi / sum(pi)
    
    names(pi) = c('Tail', Bt)
    
    if (pi['Tail'] > pi_tail_bounds[1] &
        pi['Tail'] < pi_tail_bounds[2] &
        all(pi > pi_min))
      break
  }
  
  # pio::pioStr("Mixing proportions", '\n')
  # print(pi)
  
  
  # Betas
  repeat {
    means = runif(K_betas)
    vars = runif(K_betas) / Beta_variance_scaling
    
    separations = abs(apply(expand.grid(means, means), 1, diff))
    separations = separations[separations > 0]
    
    if (all(separations > Betas_separation) &
        min(means) > Beta_bounds[1] &
        max(means) < Beta_bounds[2])
      break
    
  }
  
  a = b = NULL
  for (be in 1:K_betas)  {
    ab = mobster:::.estBetaParams(means[be], vars[be])
    a = c(a, ab$a)
    b = c(b, ab$b)
  }
  names(a) = names(b) = Bt
  
  # pio::pioStr("Betas (mean and variance)", '\n')
  # print(means)
  # print(vars)
  
  shape = runif(1, shape_bounds[1], shape_bounds[2])
  
  # pio::pioStr("Tail (shape and scale)", '\n')
  # print(shape)
  # print(scale)
  
  samples = rdbpmm(x = NULL, a, b, pi, shape, scale, n = N)
  names(samples) = NULL
  
  Nk = round(N * pi)
  ids_tail = rep('Tail', Nk['Tail'])
  ids_B = sapply(Bt, function(x)
    rep(x, Nk[x]))
  # ids_B = Reduce(c, ids_B)
  
  ids = c(unlist(ids_B), ids_tail)[1:N]
  
  plot = ggplot(data.frame(x = samples), aes(x, fill = ids)) +
    geom_histogram(binwidth = 0.01) +
    geom_vline(xintercept = means) +
    labs(title = 'MOBSTER synthetic dataset',
         subtitle = paste0('N = ', N)) +
    my_ggplot_theme() +
    guides(fill = guide_legend('Component'))
  
  input = data.frame(VAF = samples,
                     simulated_cluster = ids,
                     stringsAsFactors = FALSE) %>% as_tibble()
  input = input[!is.na(input$simulated_cluster),]
  
  return(list(
    data = input,
    model = list(
      a = a,
      b = b,
      shape = shape,
      scale = scale,
      pi = pi
    ),
    plot = plot
  ))
}
