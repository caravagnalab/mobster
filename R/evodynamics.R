# Estimate mutation rate from MOBSTER fit
#
# Mutation rate is calculated based on population genetics model,
# for details see Williams et al. In summary the mutation rate scaled
# by the probability of lineage survival \eqn{\beta}, \eqn{\mu/\beta} is given by:
# \deqn{\mu/\beta = M / (1/fmin - 1/fmax)}
# Where \eqn{fmin} is the minimum VAF and \eqn{fmax} is the maximum, and
# \eqn{M} is the number of mutations between \eqn{fmin} and \eqn{fmax}.
#
# @param fit Mobster fit
# @param lq lower quantile of VAF (0.05)
# @param uq upper quantile of VAF (0.95)
# @param ploidy of mutations (2)
# @param ncells Number of cells that accumulate mutations at each division 1 or 2, default is 2
# @return Mutation rate per tumour doubling
# @examples
# mutationrate(mobsterfit)
# @export
mutationrate <-
  function(fit,
           lq = 0.05,
           uq = 0.95,
           ploidy = 2,
           ncells = 2) {
    VAFvec <- dplyr::filter(fit$best$data, cluster == 'Tail') %>%
      dplyr::filter(VAF < quantile(VAF, uq) &
                      VAF > quantile(VAF, lq)) %>%
      dplyr::pull(VAF)
    mu <-
      length(VAFvec) / (1 / (ploidy * min(VAFvec)) - 1 / (ploidy * max(VAFvec)))
    return(mu / ncells)
  }

# Extract relevent parameters from MOBSTER fit
#
# Extract the number of mutations in the subclone, the frequency of the subclone
# and calculate the time the subclone emerges.
#
# @param fit Mobster fit
# @param mu Mutation rate
# @param subclonenumber ID of subclone
# @return tibble with all paramteres
# @examples
# subcloneparameters(mobsterfit, mutationrate(mobsterfit), 1)
subcloneparameters <- function(fit, mu, subclonenumber = 1) {
  clusters <- fit$best$Clusters %>%
    dplyr::select(-init.value) %>%
    dplyr::filter(cluster != "Tail") %>%
    tidyr::spread(type, fit.value) %>%
    dplyr::arrange(Mean)

  subclonemutations <-
    clusters$`Mixing proportion`[subclonenumber] * fit$best$N
  subclonefrequency <-
    clusters$Mean[subclonenumber] * 2 #need ccf so times by 2
  time <-
    (subclonemutations / mu) / (2 * log(2)) - (-digamma(1) / log(2))

  res <- tibble::tibble(
    time = time,
    subclonefrequency = subclonefrequency,
    subclonemutations = subclonemutations,
    cluster = clusters$cluster[subclonenumber]
  )

  return(res)
}

# Calculate strength of selection of subclone
#
# Use properties of subclone fit to calculate selection intensity, selection
# is defined as the relative growth rates of host tumour cell
# populations (\eqn{\lambda h}) vs subclone (\eqn{\lambda s}):
# \deqn{1+s=\lambda h / \lambda s}
#
# @param time time subclone emerges (in tumour doublings)
# @param time_end Time when tumour is sampled (in tumour doublings)
# @param subclonefrequency Frequency of subclones
# @return s
# @examples
# subcloneparameters(mobsterfit, mutationrate(mobsterfit), 1)
selection <- function(time, time_end, subclonefrequency) {
  x1 <- log(2) * time
  x2 <- log(subclonefrequency / (1 - subclonefrequency))
  x3 <- log(2) * (time_end - time)

  s <- ((x1 + x2) / x3)
  return(s)
}

# Calculate strength of selection for 2 independent subclones
#
# Use properties of subclone fit to calculate selection intensity, selection
# is defined as the relative growth rates of host tumour cell
# populations (\eqn{\lambda h}) vs subclone (\eqn{\lambda s}):
# \deqn{1+s=\lambda h / \lambda s}
#
# @param time1 time subclone 1 emerges (in tumour doublings)
# @param time2 time subclone 2 emerges (in tumour doublings)
# @param time_end Time when tumour is sampled (in tumour doublings)
# @param subclonefrequency1 Frequency of subclone 1
# @param subclonefrequency2 Frequency of subclone 2
# @return s
# @examples
#
selection2clone <- function(time1,
                            time2,
                            time_end,
                            subclonefrequency1,
                            subclonefrequency2) {
  x1 <- log(2) * time1
  x2 <-
    log(subclonefrequency1 / (1 - subclonefrequency1 - subclonefrequency2))
  x3 <- log(2) * (time_end - time1)
  s1 <- ((x1 + x2) / x3)

  x1 <- log(2) * time2
  x2 <-
    log(subclonefrequency2 / (1 - subclonefrequency1 - subclonefrequency2))
  x3 <- log(2) * (time_end - time2)
  s2 <- ((x1 + x2) / x3)

  return(c(s1, s2))
}

# Calculate strength of selection for 2 nested subclones
#
# Use properties of subclone fit to calculate selection intensity,
# selection is defined as the relative growth rates of host tumour cell
# populations (\eqn{\lambda h}) vs subclone (\eqn{\lambda s}):
# \deqn{1+s=\lambda h / \lambda s}
#
# @param time1 time subclone 1 emerges (in tumour doublings)
# @param time2 time subclone 2 emerges (in tumour doublings)
# @param time_end Time when tumour is sampled (in tumour doublings)
# @param subclonefrequency1 Frequency of subclone 1
# @param subclonefrequency2 Frequency of subclone 2
# @return s
# @examples
#
selection2clonenested <- function(time1,
                                  time2,
                                  time_end,
                                  subclonefrequency1,
                                  subclonefrequency2) {
  x1 <- log(2) * time1
  x2 <- log((subclonefrequency1 - subclonefrequency2)
            / (1 - subclonefrequency1))
  x3 <- log(2) * (time_end - time)
  s1 <- ((x1 + x2) / x3)

  x1 <- log(2) * time2
  x2 <- log((subclonefrequency2)
            / (1 - subclonefrequency1))
  x3 <- log(2) * (time_end - time2)
  s2 <- ((x1 + x2) / x3)

  return(c(s1, s2))
}

#' Extract evolutionary parameters from a MOBSTER fit
#'
#' @description  Mutation rate, time of emergence and selection coefficient of subclones are calculated.
#' These values are calculated based on a population genetics model of tumour evolution
#' see Williams et al. 2016 and 2018 for more details (Nature Genetics).
#'
#' The mutation rate scaled by the probability of lineage survival \eqn{\beta}, \eqn{\mu/\beta} is given by:
#' \deqn{\mu/\beta = M / (1/fmin - 1/fmax)}
#' Where \eqn{fmin} is the minimum VAF and \eqn{fmax} is the maximum, and
#' \eqn{M} is the number of mutations between \eqn{fmin} and \eqn{fmax}.
#'
#' Selection is defined as the relative growth rates of host tumour cell
#' populations (\eqn{\lambda h}) vs subclone (\eqn{\lambda s}):
#' \deqn{1+s=\lambda h / \lambda s}
#'
#' @param x An object fit by MOBSTER.
#' @param Nmax Time when tumour is sampled (in tumour doublings)
#' @param lq lower quantile of VAF (0.05)
#' @param uq upper quantile of VAF (0.95)
#' @param ploidy of mutations (2)
#' @param ncells Number of cells that accumulate mutations at each division 1 or 2, default is 1
#' @return Mutation rate, time of emergence and selection coefficient of subclones.
#'
#' @export
#'
#' @examples
#' data('fit_example', package = 'mobster')
#' evolutionary_parameters(fit_example)
evolutionary_parameters <-
  function(x,
           Nmax = 10 ^ 10,
           lq = 0.1,
           uq = 0.9,
           ploidy = 2,
           ncells = 2){
    
    if(class(x$best) == "dbpmmh"){
      return(evolutionary_parameters_mobsterh(x))
    }
    
    is_list_mobster_fits(x)
    fit = x

    if (fit$best$fit.tail == FALSE)
      stop("No tail detected,
           evolutionary inference not possible")

    nsubclones <- fit$best$Kbeta - 1 #remove 1 for clonal mutations
    mu <-
      mutationrate(
        fit,
        lq = lq,
        uq = uq,
        ploidy = ploidy,
        ncells = ncells
      )
    powerlawexponent <- fit$best$shape + 1

    if (nsubclones == 0) {
      res <- tibble::tibble(mu = mu, exponent = powerlawexponent)

    } else if (nsubclones == 1) {
      scparams <- subcloneparameters(fit, mu, 1)
      time_end <-
        log(Nmax * (1 - scparams$subclonefrequency)) / log(2)
      scparams$s <-
        selection(scparams$time, time_end, scparams$subclonefrequency)
      res <-
        dplyr::bind_cols(tibble::tibble(mu = mu, exponent = powerlawexponent),
                         scparams)

    } else if (nsubclones == 2) {
      scparams1 <- subcloneparameters(fit, mu, 1)
      scparams2 <- subcloneparameters(fit, mu, 2)
      #check if subclones satisfy pigeon hole principle, if yes assume subclones are independent

      if (scparams1$subclonefrequency + scparams2$subclonefrequency < 1) {
        time_end <- log(Nmax * (
          1 - scparams1$subclonefrequency -
            scparams2$subclonefrequency
        )) / log(2)
        s <- selection2clone(
          scparams1$time,
          scparams2$time,
          time_end,
          scparams1$subclonefrequency,
          scparams2$subclonefrequency
        )
        scparams1$s <- s[1]
        scparams2$s <- s[2]
        scparams1$subclone <- "subclone1"
        scparams2$subclone <- "subclone2"
        res <- dplyr::bind_cols(tibble::tibble(mu = mu,
                                               exponent = powerlawexponent),
                                scparams1) %>%
          dplyr::bind_rows(dplyr::bind_cols(
            tibble::tibble(mu = mu, exponent = powerlawexponent),
            scparams2
          ))
      } else {
        largestsubclone <- max(scparams1$subclonefrequency,
                               scparams2$subclonefrequency)
        time_end <- log(Nmax * (1 - largestsubclone)) / log(2)
        s <- selection2clonenested(
          scparams1$time,
          scparams2$time,
          time_end,
          scparams1$subclonefrequency,
          scparams2$subclonefrequency
        )
        scparams1$s <- s[1]
        scparams2$s <- s[2]
        scparams1$subclone <- "subclone1"
        scparams2$subclone <- "subclone2"
        res <-
          dplyr::bind_cols(tibble::tibble(mu = mu, exponent = powerlawexponent),
                           scparams1) %>%
          dplyr::bind_rows(dplyr::bind_cols(
            tibble::tibble(mu = mu, exponent = powerlawexponent),
            scparams2
          ))
      }
    } else if (nsubclones > 2) {
      res <- tibble::tibble(mu = mu) #need to extend to 3+ subclones
    }
    return(res)
  }


subcloneparameters_mobsterh <- function(x, mu, subclonenumber = "S1"){
  

  
  subclonemutations <-
    x$data %>% group_by(cluster) %>% summarize(N = n()) %>% filter(cluster == !!subclonenumber) %>% pull(N)
  subclonefrequency <-
    get_ccf_subclones(x) %>% filter(cluster == !!subclonenumber) %>% pull(CCF)
  time <-
    (subclonemutations / mu) / (2 * log(2)) - (-digamma(1) / log(2))
  
  res <- tibble::tibble(
    time = time,
    subclonefrequency = subclonefrequency,
    subclonemutations = subclonemutations,
    cluster = subclonenumber
  )
  
  return(res)
}

evolutionary_parameters_mobsterh <- function(x,
                                             Nmax = 10 ^ 10,
                                             lq = 0.1,
                                             uq = 0.9,
                                             ncells = 2){
  
  if (!x$run_parameters$tail)
    stop("No tail detected,
           evolutionary inference not possible")
  
  mu_post <- mu_posterior(x, get_genome_length(x))
  mu <- mu_post$mean * sum((get_genome_length(x) %>% filter(karyotype %in% names(x$model_parameters)) %>%  pull(length)))
  
  powerlawexponent <- x$model_parameters[[1]]$tail_shape + 1
  nsubclones <- x$run_parameters$K
  
  if (nsubclones == 0) {
    res <- list(mu = mu_post, exponent = powerlawexponent)
    
  } else if (nsubclones == 1) {
    scparams <- subcloneparameters_mobsterh(x, mu, "S1")
    time_end <-
      log(Nmax * (1 - scparams$subclonefrequency)) / log(2)
    scparams$s <-
      selection(scparams$time, time_end, scparams$subclonefrequency)
    res <- list(mu = mu_post, exponent = powerlawexponent, selection_S1 = scparams)

    
  } else if (nsubclones == 2) {
    scparams1 <- subcloneparameters_mobsterh(x, mu, "S1")
    scparams2 <- subcloneparameters_mobsterh(x, mu, "S2")
    #check if subclones satisfy pigeon hole principle, if yes assume subclones are independent
    
    if (scparams1$subclonefrequency + scparams2$subclonefrequency < 1) {
      time_end <- log(Nmax * (
        1 - scparams1$subclonefrequency -
          scparams2$subclonefrequency
      )) / log(2)
      s <- selection2clone(
        scparams1$time,
        scparams2$time,
        time_end,
        scparams1$subclonefrequency,
        scparams2$subclonefrequency
      )
      scparams1$s <- s[1]
      scparams2$s <- s[2]
      scparams1$subclone <- "subclone1"
      scparams2$subclone <- "subclone2"
      res <- list(mu = mu_post, exponent = powerlawexponent, selection_S1 = scparams1,selection_S2 = scparams2)
      
    } else {
      largestsubclone <- max(scparams1$subclonefrequency,
                             scparams2$subclonefrequency)
      time_end <- log(Nmax * (1 - largestsubclone)) / log(2)
      s <- selection2clonenested(
        scparams1$time,
        scparams2$time,
        time_end,
        scparams1$subclonefrequency,
        scparams2$subclonefrequency
      )
      scparams1$s <- s[1]
      scparams2$s <- s[2]
      scparams1$subclone <- "subclone1"
      scparams2$subclone <- "subclone2"
      res <-
        list(mu = mu_post, exponent = powerlawexponent, selection_S1 = scparams1,selection_S2 = scparams2)
    }
  } else if (nsubclones > 2) {
    cli::cli_alert_info("Selection coefficient calculation is currently implemented only for max 2 subclones! ")
    res <- list(mu = mu_post, exponent = powerlawexponent)
  }
  return(res)
  
    
  
  
}


#' Estimate mutation rate posterior from MOBSTER fit
#'
#'  @description The parameters defining the posterior distribution of the mutation rate are computed, together with the plot
#'  of the distribution. The accumulation of subclonal mutations on the genome is described by a Poisson
#'  process with intensity \Lambda=\mu\sum_{i \in mathrm{karyotypes}}(1/f_{min;i}-1/f_{max;i})\ell_i,
#'  where \mu is the effective mutation rate, f_{min;i},f_{max;i} are respectively the minimum and maximum frequency
#'  of the mutations assigned to the tail and \ell_i is the length of the genomic region for a given karyotype.
#'  Assuming a gamma distribution with parameters alpha,beta as prior, the posterior distribution is again a
#'  gamma distribution with parameters \alpha^{\prime} = \alpha + \sum_{i in \mathrm{karyotypes}}M_i,
#'  \beta = \beta + \sum(1/f_{min;i}-1/f_{max;i})\ell_i, where M_i are the number of subclonal mutations for a given karyotype.
#
#'  @param fit Fit by MOBSTERh
#'  @param prior Object containing the alpha and beta parameters of the gamma prior
#'  @param genome_length
#'  @return a list containing the parameters of the posterior distribution and a density plot
#'  @examples
#'  data('fit_example_mobsterh', package = 'mobster')
#'  prior=list(alpha=10^-4,beta=10^-4)
#'  mu_posterior(fit_example_mobsterh,prior=prior)
#'  @export
mu_posterior <- function(fit,
                         genome_length = get_genome_length(fit),
                         quantiles = c(0.02,0.98),
                         prior = list(alpha = 1e-4, beta = 1e-4)){

  # check tail
  if (fit$run_parameters$tail == 0) {
    stop("No Tail")
  }

  required_karyotypes = names(fit$model_parameters)

  # check that genome_length contains all the names required
  stopifnot(is.data.frame(genome_length))
  stopifnot("karyotype" %in% names(genome_length))
  stopifnot("length" %in% names(genome_length))

  available_karyotypes = genome_length$karyotype

  missing_length = setdiff(required_karyotypes, available_karyotypes)

  if(missing_length %>% length > 0)
  {
    warning("Missing genome length information for ", paste(missing_length, collapse = ', '))
    stop("Will not compute mutation rate")
  }

  # get subloclonal mutations, min/max frequency and chromosome lenght for each karyotype

  subclonal_mutations = c()
  f_min = c()
  f_max = c()
  length_karyo = c()

  for (karyo in required_karyotypes)
  {
    tail_mutations = fit$model_parameters[karyo][[1]]$cluster_probs[1, ] %>% sum()
    subclonal_mutations = c(subclonal_mutations, tail_mutations %>% sum())

    f_min = c(f_min, fit$model_parameters[karyo][[1]]$tail_scale)
    alpha = fit$model_parameters[karyo][[1]]$beta_concentration1[1]
    beta = fit$model_parameters[karyo][[1]]$beta_concentration2[1]
    f_max = c(f_max, alpha / (alpha + beta))

    karyotype = fit$data %>% filter(karyotype == karyo)

    # L_karyo = genome_length %>% filter()

    # Old code - not required anymore
    # segment_ids <- karyotype$segment_id %>% unique()
    # segment_ids <-
    #   read.table(text = segment_ids,
    #              sep = ":",
    #              as.is = TRUE)
    # length_karyo = c(length_karyo, sum(segment_ids$V3 - segment_ids$V2))

    length_karyo = c(length_karyo,
                     genome_length %>% filter(karyotype == karyo) %>% pull(length))

  }

  # calculate alpha e beta mutation rate posterior
  ccf = get_ccf_subclones(fit)
  
  if(!is.null(ccf)){
    
     pi = ifelse(sum(ccf$CCF) > 1, max(ccf$CCF),sum(ccf$CCF))
    
  }else{ pi = 0 }
  
  alpha = prior$alpha + sum(subclonal_mutations)*(1-pi)
  beta =  prior$beta + sum((1 / f_min - 1 / f_max) * length_karyo)
  mean = alpha / beta
  var = alpha / (beta ^ 2)
  sampling = rgamma(10000, shape = alpha, rate = beta)
  q1 = quantile(sampling,quantiles[1])
  q2 = quantile(sampling,quantiles[2])

  #sampling from the posterior distribution

  plot = ggplot(data.frame(mu = sampling), aes(x = mu)) +
    geom_histogram(bins = 100, aes(y = ..density.., fill = "indianred")) +
    stat_function(
      fun = dgamma,
      colour = "cornflowerblue",
      size = 1,
      args = list(shape = alpha, rate = beta)
    ) + geom_vline(xintercept = mean, linetype = "dashed") +
    theme_bw() +
    theme(legend.position = "none")  +
    labs(title = "mutation rate posterior distribution", x = "mu", y = "Density")

  # return the results of the inference

  inference = tibble(
    alpha = alpha,
    beta = beta,
    mean = mean,
    var = var,
    lower_quantile = q1,
    upper_quantile = q2,
    plot = list(plot)
  )

  return(inference)

}


#' Estimate the parameters of the mutation rate prior from a mobster fit.
#'
#' @description It is assumed a gamma distribution as prior for the mutation rate. The parameters alpha and beta are estimated from
#' the data. We compute an estimator of mu for each karyotype, i.e. \mu_i= M_i/((1/f_{min;i}-1/f_{max;i})\ell_i) and take the mean and
#' variance accross different karyotypes. Using the MM formula for the gamma distribution we get an expression for
#' alpha and beta parameters.
#'
#' @param fit HMobster fit
#' @return a list containing the parameters of the prior distribution
#' @examples
#' data('fit_example_mobsterh', package = 'mobster')
#' estimate_prior(fit_example_mobsterh)

estimate_prior <- function(fit) {
  # get subloclonal mutations, min/max frequency and chromosome lenght for each karyotype

  subclonal_mutations = c()
  f_min = c()
  f_max = c()
  length_karyo = c()

  for (karyo in names(fit$best$model_parameters)) {
    tail_mutations = fit$best$model_parameters[karyo][[1]]$cluster_probs[1, ] %>% sum()
    subclonal_mutations = c(subclonal_mutations, tail_mutations %>% sum())
    f_min = c(f_min, fit$best$model_parameters[karyo][[1]]$tail_scale)
    alpha = fit$best$model_parameters[karyo][[1]]$beta_concentration1[1]
    beta = fit$best$model_parameters[karyo][[1]]$beta_concentration2[1]
    f_max = c(f_max, alpha / (alpha + beta))
    karyotype = fit$best$data %>% filter(karyotype == karyo)
    segment_ids <- karyotype$segment_id %>% unique()
    segment_ids <-
      read.table(text = segment_ids,
                 sep = ":",
                 as.is = TRUE)
    length_karyo = c(length_karyo, sum(segment_ids$V3 - segment_ids$V2))
  }

  #parameters estimates

  mu = mean(subclonal_mutations / ((1 / f_min - 1 / f_max) * length_karyo))
  var = sd(subclonal_mutations / ((1 / f_min - 1 / f_max) * length_karyo))^2

  alpha = (mu ^ 2) / var
  beta = mu / var

  prior = list(alpha = alpha, beta = beta)

  return(prior)

}


get_genome_length = function(fit, exome = FALSE, build = "hg38", karyotypes = NULL){
  
  
  id = fit$data %>% mutate(segment_id = paste0(segment_id,":",karyotype)) %>% 
    dplyr::select(segment_id) %>% unique()
  
  seg <- read.table(text = id$segment_id, sep = ":", as.is = TRUE) %>% mutate(V6 = paste0(V4,":",V5))
  
  seg <- seg %>% dplyr::select(V1,V2,V3,V6)
  
  colnames(seg) <- c("chr", "from", "to","karyotype")
  
  
  if(exome){
    if(build == "hg38"){
      if(!require(TxDb.Hsapiens.UCSC.hg38.knownGene)) stop("Plaease install TxDb.Hsapiens.UCSC.hg38.knownGene to use the exonic feature")
      exs <- exons(TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene) 
    } else if(build == "hg19") {
      if(!require(TxDb.Hsapiens.UCSC.hg19.knownGene)) stop("Plaease install TxDb.Hsapiens.UCSC.hg19.knownGene to use the exonic feature")
      exs <- exons(TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene) 
    } else {
      stop("Supporterd build are hg38 and hg19!")
    }
   
    seg_grang <- lapply(split(seg, seg$karyotype), FUN = function(x) GenomicRanges::makeGRangesFromDataFrame(x,start.field = "from", end.field = "to"))
    seg_int <- lapply(seg_grang, function(x) GenomicRanges::intersect(x, exs, ignore.strand = TRUE) %>% as.data.frame)
    seg_int <- mapply(names(seg_int), seg_int, FUN = function(x,y) {y$karyotype <- x 
    y
    }, SIMPLIFY = F) %>% do.call(rbind,.)
    seg <- seg_int %>% dplyr::select(seqnames, start,end, karyotype) %>% dplyr::rename(chr = seqnames, from = start, to = end)
  }
  
  if(!is.null(karyotypes)) seg %>% filter(karyotype %in% karyotypes)
  
  lengths = seg  %>% group_by(karyotype) %>% 
    summarize(length = sum(to-from)) 
  
  return(lengths)
  
}

#' Estimate selection coefficient posterior for a subclone from MOBSTER fit
#'
#'  @description Posterior distribution of time of origin t and selection coefficient s of a subclone.
#'  We assume a Poisson distribution as likelihood for the number of diploid heterozygous mutations of the subclonal cluster with mean M = 2\mu l\omega t,
#'  where l is the length of the diploid genome, \mu is the mutation rate and \omega is the growth rate of the tumour. 
#'  We assume a Beta distribution as likelihood for the mean VAF of the subclonal cluster, where the mean corresponds to the CCF of the subclone divided 
#'  by 2 and computed assuming exponential growth for both populations. The time of origin is expressed in tumour doublings.
#'  
#
#'  @param fit Fit by MOBSTERh
#'  @param prior_s Object containing range of values of s and the corresponding probabilities
#'  @param N_max Object containing 
#'  @return 
#'  @examples
#'  s_posterior(my_fit)
#'  @export

selection_posterior <- function(fit,
                        N_max = 10^10
                        ){
  
  # check subclone
  if (! "S1" %in% (fit$data$cluster %>% unique())){
    stop(paste0("No ",subclone," cluster"))
  }
  

library(stringr)
  
mu = mobster:::mu_posterior(fit = fit, genome_length = mobster:::get_genome_length(fit)) %>% pull(mean)

l = get_genome_length(fit)  %>% 
     mutate(ploidy = strsplit(karyotype,":")[[1]] %>% as.numeric() %>% sum())

n_subclones = grep(x = fit$data$cluster %>% unique(),pattern = "S") %>% length()

if(n_subclones == 1){
  
  M = fit$data %>% filter(!is.na(cluster),cluster == "S1") %>% as_tibble() %>% 
     group_by(karyotype) %>% 
      summarize(n = length(chr)) 
  
  
  ccf =    get_ccf_subclones(fit) %>% pull(CCF)
  
  nu = lapply(names(fit$model_parameters), 
              function(k){tibble(karyotype = k,n = fit$model_parameters[[k]]$n_trials_subclonal)}) %>% 
       bind_rows() %>% pull(n) %>% mean()
      
# likelihood calculation
  
t = rgamma(n=1e4, shape = M$n %>% sum(),rate = log(2)*mu*sum(l$length*l$ploidy)) 

f = rbeta(n=1e4, shape1 = nu*ccf, shape2 = (1-ccf)*nu)
 
Time = log((1-ccf)*N_max)/log(2) - digamma(1)/log(2)
 
s = (log(f/(1-f)) + log(2)*t)/(log(2)*(Time-t))

return(tibble(param = c("t","s"),mean = c(mean(t),mean(s)), var =  c(var(t),var(s)),
         inference = c(list(t),list(s))))

 }else{
    
   ccf = get_ccf_subclones(fit) 
   
   M = fit$data %>% filter(!is.na(cluster),cluster %in% c("S1","S2")) %>% as_tibble() %>% 
     group_by(cluster,karyotype) %>% summarize(n = length(chr))
   
   nu =   lapply(names(fit$model_parameters), 
          function(k){tibble(karyotype = k,
                             S1 = fit$model_parameters[[k]]$n_trials_subclonal[1],
                             S2 = fit$model_parameters[[k]]$n_trials_subclonal[2])}) %>% 
             bind_rows() 


   if(ccf %>% pull(CCF) %>% sum() > 1){
     
   
     t1 = rgamma(n=1e4, shape = M %>% filter(cluster == "S2") %>% pull(n) %>% sum(), 
                 rate = log(2)*mu*sum(l$length*l$ploidy))
     
     Time = log((1-ccf %>% pull(CCF) %>% sum())*N_max)/log(2) - digamma(1)/log(2)
     
     f1 = rbeta(n=1e4, shape1 = mean(nu$S2)*(ccf %>% filter(cluster == "S2") %>% pull(CCF)), 
                shape2 = mean(nu$S2)*(1 - ccf %>% filter(cluster == "S2") %>% pull(CCF)))
     
     f2 = rbeta(n=1e4, shape1 = mean(nu$S1)*(ccf %>% filter(cluster == "S1") %>% pull(CCF)), 
                shape2 = mean(nu$S1)*(1 - ccf %>% filter(cluster == "S1") %>% pull(CCF)))
     
     s1 = (log(max(f1-f2,0.01)/(1-f1)) + log(2)*t1)/(log(2)*(Time-t1))
     
     t2 = rgamma(n=1e4, shape = M %>% filter(cluster == "S1") %>% pull(n) %>% sum(), 
                 rate = log(2)*(1+s1)*mu*sum(l$length*l$ploidy)) + t1
     
     s2 = (log(f2/max(1-f1,0.01)) + log(2)*t2)/(log(2)*(Time-t2))
     
     
}else{
     
     t1 = rgamma(n=1e4, shape = M %>% filter(cluster == "S1") %>% pull(n) %>% sum(), 
                 rate = log(2)*mu*sum(l$length*l$ploidy)) 
     t2 = rgamma(n=1e4, shape = M %>% filter(cluster == "S2") %>% pull(n) %>% sum(), 
                 rate = log(2)*mu*sum(l$length*l$ploidy)) 
     
     f1 = rbeta(n=1e4, shape1 = mean(nu$S1)*(ccf %>% filter(cluster == "S1") %>% pull(CCF)), 
                       shape2 = mean(nu$S1)*(1 - ccf %>% filter(cluster == "S1") %>% pull(CCF)))
     
     f2 = rbeta(n=1e4, shape1 = mean(nu$S2)*(ccf %>% filter(cluster == "S2") %>% pull(CCF)), 
                       shape2 = mean(nu$S2)*(1 - ccf %>% filter(cluster == "S2") %>% pull(CCF)))
     
     Time = log((1-ccf %>% pull(CCF) %>% sum())*N_max)/log(2) - digamma(1)/log(2)
     
     s1 = (log(f1/max(1-f1-f2,0.01)) + log(2)*t1)/(log(2)*(Time-t1))
     
     s2 = (log(f2/max(1-f1-f2,0.01)) + log(2)*t2)/(log(2)*(Time-t2))
     
}
 
   return(tibble(param = c("t1","t2","s1","s2"),
                 mean = c(mean(t1),mean(t2),mean(s1),mean(s2)), 
                 var = c(var(t1),var(t2),var(s1),var(s2)),
                 inference = c(list(t1),list(t2),list(s1),list(s2))))  
 
   }


}

