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
  function(x,
           Nmax = 10^10,
           lq = 0.1,
           uq = 0.9,
           ploidy = 2,
           ncells = 2) 
  {
    
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
