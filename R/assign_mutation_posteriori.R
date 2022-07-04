
#' Assign mutation a posteriori
#'
#' A function to assign filtered and non-subsetted mutations in used karyotypes
#' @param x mobsterh fit object
#'
#' @return a mobsterh object with new cluster annotations and a new field posteriori_measures with cluster probabilities
#' @export
#'
#' @examples
#' library(mobster)
#' library(dplyr)
#' assign_mutations_posteriori(mobster::fit_example_mobsterh$best)$data %>% filter(posteriori_annot) %>% head(5)


assign_mutations_posteriori <- function(x) {
  
  to_assign <- x$data %>% filter(is.na(cluster), karyotype %in% names(x$model_parameters)) 
  
  if(nrow(to_assign) == 0) {
    cli::cli_alert_info("No mutations to assign, returning input object!")
    return(x)
  }
  
  x$data$posteriori_annot <- FALSE
  
  by_karyo <- split_and_tolist(x = to_assign)
  
  ret <- lapply(1:length(by_karyo), function(i) assign_mutations_posteriori_aux(by_karyo[[i]], x, names(by_karyo)[i] )) 
  clusters <- lapply(ret, function(z) z$clusters) %>% do.call(c, .)
  assignment_probs <- lapply(ret, function(z) z$assignment_probs) %>% do.call(rbind, .)
  x$data[match(names(clusters), x$data$id), ]$cluster <- clusters
  x$data[match(names(clusters), x$data$id), ]$posteriori_annot <- TRUE
  x$posteriori_measures$assignment_probs <- assignment_probs
  
  return(x)
}


assign_mutations_posteriori_aux <- function(muts,fit_object, karyo) {
  
  betas <- get_beta(fit_object) %>% filter(karyotype == !!karyo)
  lk_clonal <- apply(betas,1, FUN = function(x) VGAM::dbetabinom.ab(muts[,1], size = muts[,2], shape1 = x[1] %>% as.numeric(), shape2 = x[2] %>% as.numeric(), log = F) * x[4] %>% as.numeric(), simplify = F)
  clonal_names <- paste0("C", 1:nrow(betas))
  
  if(fit_object$run_parameters$K >0){
    if(fit_object$run_parameters$subclonal_prior == "Beta"){
      betas_sub <- get_beta_sub(fit_object) %>% filter(karyotype == !!karyo)
      lk_sub <- apply(betas_sub,1, FUN = function(x) VGAM::dbetabinom.ab(muts[,1], size = muts[,2], shape1 = x[1] %>% as.numeric(), shape2 = x[2] %>% as.numeric(), log = F) * x[5] %>% as.numeric(), simplify = F)
      subclonal_names <- paste0("S", 1:nrow(betas_sub))
    } else {
      moyal_sub <- get_moyal_sub(fit_object) %>% filter(karyotype == !!karyo)
      lk_sub <- apply(moyal_sub,1, FUN = function(x) dtmoyalbin(muts, loc = x[1] %>% as.numeric(), scale = x[2] %>% as.numeric(), lower = 0, upper = min(betas$a / (betas$a + betas$b))) * x[5] %>% as.numeric(), simplify = F)
      subclonal_names <- paste0("S", 1:nrow(moyal_sub))
    }
    
  } else {
    lk_sub <- list(rep(0, nrow(muts)))
    subclonal_names <- "S1"
  }
  
  if(fit_object$run_parameters$tail > 0){
    tail <- get_pareto(fit_object) %>% filter(karyotype == !!karyo)
    tot <- as.numeric(substr(karyo,1,1)) + as.numeric(substr(karyo,3,3))
    lk_tail <- list( apply(muts, 1,function(mut)  integrate(pareto_beta_bin_mix, lower = tail$scale, upper = tail$location,
                           NV = mut[1], DP = mut[2],
                           shape = tail$shape * tot, scale = tail$scale)$value * tail$mixing, simplify = F) %>% unlist() )
  } else {
    lk_tail <- list(rep(0, nrow(muts)))
    
  }
    
  prob_ass <- do.call(rbind, c(lk_clonal, lk_sub, lk_tail))
  names_cluster <- c(clonal_names, subclonal_names, "Tail")
  cluster <- names_cluster[apply(prob_ass,2, which.max)]
  prob_ass <- apply(prob_ass,1, function(z) z / colSums(prob_ass), simplify = FALSE) %>% do.call(cbind,.)
  names(cluster) <- rownames(prob_ass)
  colnames(prob_ass) <- names_cluster
  
  return(list(clusters = cluster, assignment_probs = prob_ass) )
}


dtmoyalbin <- function(muts, loc, scale, lower, upper){
  apply(muts,1, function(mut) integrate(moyal_bin_mix, lower = lower, upper = upper,
            NV = mut[1], DP = mut[2],
            scale = scale, loc = loc, upper_moyal = upper, lower_moyal = lower)$value, simplify = F) %>% unlist() %>% return(.)
}

moyal_bin_mix <- function(p, NV, DP, loc, scale, lower_moyal = 0, upper_moyal = 1){
  
  dbinom(NV, DP, p) * dtruncmoyal(p,loc = loc, scale = scale,lower = lower_moyal, upper = upper_moyal, log = F)
  
}



