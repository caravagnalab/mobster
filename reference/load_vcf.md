# Load data from a VCF file

This function parses a set of mutations in the Variant Calling Format
(VCF), using the `vcfR` R package. The function takes in input the
column names of the INFO field referring to allele depth counts, and
number of reads with the variant. The two columns are used to compute
the variant allele frequency (VAF) value; if any of the required columns
are not found, an error is returned. This function does not filter any
of the variants annotated in the input VCF file.

## Usage

``` r
load_vcf(file, DP_column = "DP", NV_column = "NV")
```

## Arguments

- file:

  A VCF file name.

- DP_column:

  Column with the total depth (DP), i.e., sum of reads at a locus.

- NV_column:

  Column with the number of reads with the variant (NV).

## Value

A tibble with the loaded data which contains, besides all the content
parsable from the VCF file, three columns named DP, NV and VAF where VAF
= NV/DP.

## Examples

``` r
# Example VCF file in the https://github.com/openvax/varcode repository
file = 'https://raw.githubusercontent.com/openvax/varcode/master/test/data/strelka-example.vcf'

download.file(url = file, destfile = 'strelka-example.vcf')
#> Warning: downloaded length 0 != reported length 14
#> Warning: cannot open URL 'https://raw.githubusercontent.com/openvax/varcode/master/test/data/strelka-example.vcf': HTTP status was '404 Not Found'
#> Error in download.file(url = file, destfile = "strelka-example.vcf"): cannot open URL 'https://raw.githubusercontent.com/openvax/varcode/master/test/data/strelka-example.vcf'

# We pretend that the number of variants is the gt_DP2 field, which is wrong
# Anyway, this shows how you can load a VCF file.
load_vcf(file = 'strelka-example.vcf', DP_column = 'gt_DP', NV_column = 'gt_DP2')
#> Error in load_vcf(file = "strelka-example.vcf", DP_column = "gt_DP", NV_column = "gt_DP2"): Input file not found: strelka-example.vcf
```
