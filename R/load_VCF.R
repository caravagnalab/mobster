# https://github.com/openvax/varcode/blob/master/test/data/strelka-example.vcf
# 
# require(dplyr)
# 
# file = 'strelka-example.vcf'
# DP_column = 'DP'
# NV_column = 'NV'



#' Load data from a VCF file
#' 
#' @description 
#' 
#' This function parses a set of mutations in the Variant Calling
#' Format (VCF), using the \code{vcfR} R package. The function takes
#' in input the column names of the INFO field referring to 
#' allele depth counts, and number of reads with the variant. The
#' two columns are used to compute the variant allele frequency (VAF)
#' value; if any of the required columns are not found, an error
#' is returned. This function does not filter any of the variants
#' annotated in the input VCF file.
#'
#' @param file A VCF file name.
#' @param DP_column Column with the total depth (DP), i.e., sum of
#' reads at a locus.
#' @param NV_column Column with the number of reads with the
#' variant (NV).
#'
#' @return A tibble with the loaded data which contains, besides
#' all the content parsable from the VCF file, three columns
#' named DP, NV and VAF where VAF = NV/DP.
#' 
#' @export
#'
#' @examples
#' # Example VCF file in the https://github.com/openvax/varcode repository
#' file = 'https://raw.githubusercontent.com/openvax/varcode/master/test/data/strelka-example.vcf'
#' 
#' download.file(url = file, destfile = 'strelka-example.vcf')
#' 
#' load_vcf(file = 'strelka-example.vcf', DP_column)
load_vcf = function(file,
                    DP_column = 'DP',
                    NV_column = 'NV')
{
  if (!requireNamespace("vcfR", quietly = TRUE)) {
    stop("Package \"vcfR\" needed for this function to work. Please install it.",
      call. = FALSE)
  }
  if (!file.exists(file))
    stop("Input file not found: ", file)
  
  fi = file.info(file)
  fi_size = utils:::format.object_size(fi$size, "auto")
  
  cli::cli_alert_info("Loading VCF file {.field {file}} with package vcfR, file size {.field {fi_size}}.")
  
  vcf_file = vcfR::read.vcfR(file, verbose = FALSE)
  
  # Tidy, genotypes and genome coordinates -> all data
  info_tidy = vcfR::extract_info_tidy(vcf_file)
  genotype_tidy = vcfR::extract_gt_tidy(vcf_file, verbose = FALSE)
  
  all_data = dplyr::full_join(info_tidy, genotype_tidy, by = 'Key')
  all_data = cbind(vcf_file@fix, all_data) %>% tibble::as_tibble()
  
  cn = colnames(all_data) 
  
  # Check input format
  if(!(DP_column %in% cn)) stop("Depth column ", DP_column, " is not in the VCF columns: ", paste(cn, collapse = ', '), '.')
  if(!(NV_column %in% cn)) stop("Number of variants column ", NV_column, " is not in the VCF columns: ", paste(cn, collapse = ', '), '.')
  
  all_data %>%
    dplyr::mutate(DP = !!DP_column, NV = !!NV_column, VAF = NV/VAF)
}
