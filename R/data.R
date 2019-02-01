#' Synthetic sequencing data generated from a evolutionary based cancer simulation.
#'
#' This data is generated from an evolutionary model where there is on subclonal population
#' and all other mutations are neutral passengers. The input parameters of the model are as follows,
#' MOBSTER can extract these parameters from the VAF distribution.
#' Mutation rate: 25.0
#' Number of subclones: 1
#' Subclone time: 9.0 (tumour doublings)
#' Subclone Mutations: 325
#' Subclone frequency: 0.53
#' Fitness advantage of subclone: 0.8
#' Tumour population size: 10^6
#' Clonal mutations: 500
#' 
#'
#' @format A dataframe containg a vector with variant allele frequencies (VAFs) ranging from 0 to 1
#' @source Generated using cancer sequencing simulation \url{https://github.com/marcjwilliams1/CancerSeqSim.jl}
"subclonedynamics"
