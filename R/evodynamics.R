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
  x3 <- log(2) * (time_end - time1)
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
  function(fit,
           Nmax = 10 ^ 10,
           lq = 0.1,
           uq = 0.9,
           ploidy = 2,
           ncells = 2){
    
    if(class(fit$best) == "dbpmmh"){
      return(evolutionary_parameters_mobsterh(x$best))
    }
    
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
  
  mu_post <- mu_posterior(x, get_genome_length(x),ncells)
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
                         ncells = 1,
                          lq = 0.05,
                          uq = 0.95,
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
    # tail_mutations = fit$model_parameters[karyo][[1]]$cluster_probs[1, ] %>% sum()
    # subclonal_mutations = c(subclonal_mutations, tail_mutations)
    # f_min = c(f_min, fit$model_parameters[karyo][[1]]$tail_scale)
    # alpha = fit$model_parameters[karyo][[1]]$beta_concentration1[1]
    # beta = fit$model_parameters[karyo][[1]]$beta_concentration2[1]
    # f_max = c(f_max, alpha / (alpha + beta))
    
    tail_mutations = fit$data %>% filter(karyotype == karyo, cluster == "Tail")
    tail_mutations = tail_mutations %>% filter(VAF > quantile(VAF,lq) & VAF < quantile(VAF,uq))
    subclonal_mutations = c(subclonal_mutations, tail_mutations %>% nrow())
    f_min = c(f_min,  tail_mutations %>% pull(VAF) %>% min())
    f_max = c(f_max,  tail_mutations %>% pull(VAF) %>% max())

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
  
  good_karyo = (subclonal_mutations > 0) & (f_min < f_max)
  
  subclonal_mutations = subclonal_mutations[good_karyo]
  f_min = f_min[ good_karyo]
  f_max = f_max[ good_karyo]
  
  if(sum(good_karyo) == 0)
  {
    warning("Missing karyotype with enough mutations")
    stop("Will not compute mutation rate")
  }

  # calculate alpha e beta mutation rate posterior

  alpha = prior$alpha + sum(subclonal_mutations)
  beta =  prior$beta + sum((1 / f_min - 1 / f_max) * length_karyo * ncells)
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
    plot = list(plot),
    sampling = list(sampling)
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

estimate_prior <- function(fit,ncells = 1) {
  # get subloclonal mutations, min/max frequency and chromosome lenght for each karyotype

  subclonal_mutations = c()
  f_min = c()
  f_max = c()
  length_karyo = c()

  for (karyo in names(fit$model_parameters)) {
    tail_mutations = fit$model_parameters[karyo][[1]]$cluster_probs[1, ] %>% sum()
    subclonal_mutations = c(subclonal_mutations, tail_mutations %>% sum())
    f_min = c(f_min, fit$model_parameters[karyo][[1]]$tail_scale)
    alpha = fit$model_parameters[karyo][[1]]$beta_concentration1[1]
    beta = fit$model_parameters[karyo][[1]]$beta_concentration2[1]
    f_max = c(f_max, alpha / (alpha + beta))
    karyotype = fit$data %>% filter(karyotype == karyo)
    segment_ids <- karyotype$segment_id %>% unique()
    segment_ids <-
      read.table(text = segment_ids,
                 sep = ":",
                 as.is = TRUE)
    length_karyo = c(length_karyo, sum(segment_ids$V3 - segment_ids$V2))
  }

  #parameters estimates

  mu = mean(subclonal_mutations / ((1 / f_min - 1 / f_max) * length_karyo * ncells))
  var = mu/100

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
  
  lengths = seg  %>% dplyr::group_by(karyotype) %>% 
    dplyr::summarize(length = sum(to-from)) 
  
  return(lengths)
  
}

#' Estimate selection coefficient posterior for a subclone from MOBSTER fit
#'
#'  @description Posterior distribution of time of origin t, mrca and selection coefficient s of a subclone.
#'  We assume a Poisson distribution as likelihood for the number of diploid heterozygous mutations of the subclonal cluster with mean M = 2\mu l\omega t,
#'  where l is the length of the diploid genome, \mu is the mutation rate and \omega is the growth rate of the tumour. 
#'  We assume an exponential distribution as likelihood for the mean number of cells of the subclone. The rate of the exponential is the inverse of expected
#'  subclone size exp(-\omega*(T-t)), where T is the time of sample collection. Times and growth rate are expressed in tumour doubling units. 
#'  Posterior distributions are computed with HMC using RSTAN.
#'  
#
#'  @param fit Fit by MOBSTERh
#'  @param N_max Final tumour size
#'  @param ncells Model of cell division
#'  @param u = 0   mean of lognormal prior on s
#'  @param sigma = 1 sigma of lognormal prior on s
#'  @return stan inference and a plot of prior and posterior distributions
#'  @export
#'  @examples
#'  data('fit_example', package = 'mobster')
#'  evo_dynamics = s_posterior(fit_example,N_max = 10^9,ncells = 2, u = 0, sigma = 0.5)

s_posterior = function(fit,N_max = 10^9,ncells = 1,u1 = -0.5,sigma1 = 0.6,u2 = -0.5,sigma2 = 0.6){
  
  
  if(class(fit) == "dbpmmh"){
    
    order_karyo = c(names(fit$model_parameters)[names(fit$model_parameters) == "1:0"],
                    names(fit$model_parameters)[names(fit$model_parameters) == "1:1"],
                    names(fit$model_parameters)[names(fit$model_parameters) == "2:0"],
                    names(fit$model_parameters)[names(fit$model_parameters) == "2:1"],
                    names(fit$model_parameters)[names(fit$model_parameters) == "2:2"])
    
    karyo = order_karyo[1]
    
    if(length(karyo) == 0){
      warning("No common karyotypes")
      stop("Will not compute s")
    }
    
    n_alleles = as.numeric(strsplit(x = karyo, split = ":")[[1]][1]) + 
               as.numeric(strsplit(x = karyo, split = ":")[[1]][2])
    m1 = fit$data %>% filter(cluster == 'S1',karyotype == karyo) %>% nrow()
    m2 = fit$data %>% filter(cluster == 'S2',karyotype == karyo) %>% nrow()
    
    if(m1 == 0){
      warning("No subclonal mutations")
      stop("Will not compute s")
    }
    
    mu = mu_posterior(fit,ncells = ncells) %>% pull(mean)
    ccf1 = (fit$data %>% filter(cluster == 'S1',karyotype == karyo) %>% 
             pull(VAF) %>% mean())*n_alleles 
    ccf2 = (fit$data %>% filter(cluster == 'S2',karyotype == karyo) %>% 
              pull(VAF) %>% mean())*n_alleles 
    length_genome = get_genome_length(fit) %>% filter(karyotype == karyo) %>% 
      pull(length)
    
  }else{
  
    n_alleles = 2
    m1 = fit$data %>% filter(cluster == 'C2') %>% nrow()
    m2 = fit$data %>% filter(cluster == 'C3') %>% nrow()
    mu = mutationrate(fit,ncells = ncells)/(2*3*10^9)
    ccf1 = (fit$data %>% filter(cluster == 'C2') %>% pull(VAF) %>% mean())*n_alleles 
    ccf2 = (fit$data %>% filter(cluster == 'C3') %>% pull(VAF) %>% mean())*n_alleles 
    length_genome = 3*10^9
    
    
}

  library(rstan)
  # library(ggplot2)
  # library(ggpubr)
  
  data = list( m1 = m1,
               m2 = m2,
               Ntot = N_max,
               ccf1 = ccf1,
               ccf2 = ccf2,
               n = n_alleles,
               length_genome = length_genome,
               mu = mu,
               u1 = u1,
               sigma1 = sigma1,
               u2 = u2,
               sigma2 = sigma2
  )
  
  if(m2 == 0){

type = "single subclone"

model.stan = "

data {
  
  int<lower=0> m1;
  int <lower=0> m2;
  real<lower=0> Ntot;
  real<lower=0,upper = 1> ccf1;
  real<lower=0,upper = ccf1> ccf2;
  real<lower=0> length_genome;
  real<lower=0> mu;
  int n;
  real u1;
  real<lower=0> sigma1;
  real u2;
  real<lower=0> sigma2;
}

parameters {
  
  real<lower=0> t_mrca;
  real<lower=0,upper=t_mrca> t_driver;
  real<lower=0> s;
  
}

model{
  
  // priors
  
  t_mrca ~ uniform(0,log(Ntot*(1-ccf1))/log(2));
  t_driver ~ uniform(0,t_mrca);
  s ~ lognormal(u1,sigma1);
  
  //  Likelihood mutations
  
  m1 ~ poisson(n*mu*log(2)*length_genome*((t_driver) + (1+s)*(t_mrca-t_driver)));
  
  // Likelihood branching process
  
  target += exponential_lpdf(Ntot*ccf1 | exp(-log(2)*(1+s)*(log(Ntot*(1-ccf1))/log(2) - t_driver)));
  
}

"
  }else{
  
    if(ccf1 + ccf2 < 1){
      
      type = "two indipendent subclones"
      
      model.stan = "

data {
  
  int<lower=0> m1;
  int <lower=0> m2;
  real<lower=0> Ntot;
  real<lower=0,upper = 1> ccf1;
  real<lower=0,upper = ccf1> ccf2;
  real<lower=0> length_genome;
  real<lower=0> mu;
  int n;
  real u1;
  real<lower=0> sigma1;
  real u2;
  real<lower=0> sigma2;
}

parameters {
  
  real<lower=0> t_mrca1;
  real<lower=0,upper=t_mrca1> t_driver1;
  real<lower=0> t_mrca2;
  real<lower=0,upper=t_mrca2> t_driver2;
  real<lower=0> s1;
  real<lower=0> s2;
  
}

model{
  
  // priors
  
  t_mrca1 ~ uniform(0,log(Ntot*(1-ccf1-ccf2))/log(2));
  t_mrca2 ~ uniform(0,log(Ntot*(1-ccf1-ccf2))/log(2));
  t_driver1 ~ uniform(0,t_mrca1);
  t_driver2 ~ uniform(0,t_mrca2);
  s1 ~ lognormal(u1,sigma1);
  s2 ~ lognormal(u2,sigma2);
  
  //  Likelihood mutations
  
  m1 ~ poisson(n*mu*log(2)*length_genome*((t_driver1) + (1+s1)*(t_mrca1-t_driver1)));
  m2 ~ poisson(n*mu*log(2)*length_genome*((t_driver2) + (1+s2)*(t_mrca2-t_driver2)));
  
  // Likelihood branching process
  
  target += exponential_lpdf(Ntot*ccf1 | 
  exp(-log(2)*(1+s1)*(log(Ntot*(1-ccf1-ccf2))/log(2) - t_driver1)));
  target += exponential_lpdf(Ntot*ccf2 | 
  exp(-log(2)*(1+s2)*(log(Ntot*(1-ccf1-ccf2))/log(2) - t_driver2)));
  
}

"
      
    }else{
      
      type = "two nested subclones"
      
      model.stan = "

data {
  
  int<lower=0> m1;
  int <lower=0> m2;
  real<lower=0> Ntot;
  real<lower=0,upper = 1 > ccf1;
  real<lower=0,upper = ccf1> ccf2;
  real<lower=0> length_genome;
  real<lower=0> mu;
  int n;
  real u1;
  real<lower=0> sigma1;
  real u2;
  real<lower=0> sigma2;
}

parameters{
  
  real<lower=0> t_sp;
  real<lower=0,upper=t_sp> t_driver1;
  real<lower=t_sp> t_mrca2;
  real<lower=t_sp,upper=t_mrca2> t_driver2;
  real<lower=0> s1;
  real<lower=s1> s2;
  
}

model{
  
  // priors
  
  t_sp ~ uniform(0,log(Ntot*(1-ccf1))/log(2));
  t_driver1 ~ uniform(0,t_sp);
  t_mrca2 ~ uniform(t_sp,log(Ntot*(1-ccf1))/log(2));
  t_driver2 ~ uniform(t_sp,t_mrca2);
  s1 ~ lognormal(u1,sigma1);
  s2 ~ lognormal(u2,sigma2) T[s1,];
  
  //  Likelihood mutations
  
  m1 ~ poisson(n*mu*log(2)*length_genome*((t_driver1) + (1+s1)*(t_sp-t_driver)));
  m2 ~ poisson(n*mu*log(2)*length_genome*((t_driver2-t_sp) + (1+s2)*(t_mrca2-t_sp)));
  
  // Likelihood branching process
  
  target += exponential_lpdf(Ntot*ccf2 | 
  exp(-log(2)*(1+s2)*(log(Ntot*(1-ccf1))/log(2) - t_driver2)));
  
  target += exponential_lpdf(Ntot*(ccf1-ccf2) | 
  exp(-log(2)*(1+s1)*(log(Ntot*(1-ccf1))/log(2) - t_driver1)));
  
}

"
    }
     
}


fit <- stan(
  model_code = model.stan,  # Stan program
  data = data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 5000,            # total number of iterations per chain
  cores = 4,
  control = list(adapt_delta= 0.99),
  verbose = F
)



# post = as.data.frame(fit) %>% dplyr::select(c("t_mrca","t_driver","s")) %>% mutate(type = "post")
# pre = data.frame(t_mrca = runif(nrow(post),min = 0,max = log(N_max*(1-ccf))/log(2)),
#                  t_driver = runif(nrow(post),min = 0,max = log(N_max*(1-ccf))/log(2)),
#                  s = rlnorm(nrow(post),meanlog = u,sdlog = sigma), 
#                  type = "prior")

# post = rbind(post,pre)
# 
# p_mrca = ggplot(post  %>% reshape2::melt() %>% filter(variable == "t_mrca"))   + labs(x = "", title = "t_mrca") + 
#   CNAqc:::my_ggplot_theme() + 
#   scale_alpha_manual(values = c("prior" = 0.4, "post" = 1)) + geom_density(aes(x = value, alpha = type), fill = "indianred")  
# 
# p_driver = ggplot(post  %>% reshape2::melt() %>% filter(variable == "t_driver"))   + labs(x = "", title = "t_driver") + 
#   CNAqc:::my_ggplot_theme() + 
#   scale_alpha_manual(values = c("prior" = 0.4, "post" = 1)) + geom_density(aes(x = value, alpha = type), fill = "darkgreen")  
# 
# p_s = ggplot(post  %>% reshape2::melt() %>% filter(variable == "s"))   + labs(x = "", title = "s") + 
#   CNAqc:::my_ggplot_theme() + 
#   scale_alpha_manual(values = c("prior" = 0.4, "post" = 1)) + geom_density(aes(x = value, alpha = type), fill = "cornflowerblue") 
# 
# 
# p_dynamics = ggarrange(plotlist = list(p_mrca,p_driver,p_s),nrow = 1, ncol = 3, legend = "bottom")
 
# return(list(posterior = fit , plot = p_dynamics))

return(list(posterior = fit , model_type = type))

}

