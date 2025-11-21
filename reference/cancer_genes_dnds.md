# List of cancer genes to compute dnds values.

List of gene ids (HuGO format) to compute dnds values. There are 2 lists
available of putative driver genes, and 2 lists of essential genes;
these lists have been subset to include only genes in the RefCDS
database.

Drivers have been compiled in \`Martincorena, et al. Cell 171.5 (2017):
1029-1041\`, and in \`Tarabichi, et al. Nature Genetics 50.12 (2018):
1630\`. Essential genes have been compiled using two different cell
lines in \`Wang et al. Science 350.6264 (2015): 1096-1101.\` and
\`Bloomen et al. Science 350.6264 (2015): 1092-1096\`. All lists are
available and named accordingly; use \`names(cancer_genes_dnds)\` to see
the available names.

## Usage

``` r
data(cancer_genes_dnds)
```

## Format

List of cancer genes to compute dnds values.

## Examples

``` r
data(cancer_genes_dnds)
print(lapply(cancer_genes_dnds, head))
#> $Martincorena_drivers
#> [1] "CCDC6"     "EIF1AX"    "HIST1H2BD" "MED12"     "POLE"      "SMARCB1"  
#> 
#> $Tarabichi_drivers
#> [1] "ACVR1"  "ACVR1B" "AKT1"   "ALK"    "AMER1"  "APC"   
#> 
#> $Wang_essentials
#> [1] "ABL1"    "RPL23A"  "AARS2"   "TRMT112" "FARSA"   "ABCB7"  
#> 
#> $Bloomen_essentials
#> [1] "AARS"     "AASDHPPT" "AATF"     "ABCB7"    "ABCE1"    "ABCF1"   
#> 
```
