#' Assign new observations to the clusters inside a MOBSTER fit.
#' 
#' @description For a new set of observations in the MOBSTER input format,
#' the parameters and clusters of a MOBSTER fit are used to determine the
#' assignments of the new points. MOBSTER density function is used to determine
#' which cluster has the highest density for each observation, and the cluster
#' label is computed accordingly. 
#'
#' @param x A MOBSTER fit.
#' @param y A MOBSTER input dataset, which has to have a VAF numeric column without NAs.
#'
#' @return The data in \code{y} is augmented with a colum per mixture component reporting
#' the corresponding density value. A final colum `cluster` is also added reporting the
#' component name for the hard clustering assignment of the point.
#' 
#' @export
#'
#' @examples
#' library(ggplot2)
#' data('fit_example', package = 'mobster')
#' 
#' # Generate some randome numbers and assign them to the most likely mixture component
#' new_assignments = Clusters_denovo(fit_example$best, data.frame(VAF = runif(1000))) 
#' print(new_assignments)
#' 
#' # Plot a histogram coloured according to the clusters
#' ggplot(new_assignments, aes(VAF, fill = cluster)) + geom_histogram(binwidth = 0.01)
Clusters_denovo = function(x, y)
{
  mobster:::is_mobster_fit(x)
  mobster:::is_mobster_input_dataset(y)
  
  # actual clusters
  components = names(mobster:::.params_Pi(x))
  
  # per-component density
  densities = y
  for(component in seq_along(components))
  {
    comp_density = ddbpmm(x, data = y$VAF, components = component, log = TRUE) %>% data.frame()
    colnames(comp_density) = components[component]
    
    densities = dplyr::bind_cols(
      densities,
      comp_density
    )
  }
  # colnames(densities) = components
  densities = densities %>% as_tibble()
  
  # Hard clustering assignments
  densities$cluster = components[apply(densities[, components, drop = FALSE], 1, which.max)]
  
  densities
  # dplyr::bind_cols(y, densities) %>% as_tibble()
}
