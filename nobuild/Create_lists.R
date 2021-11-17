im = read.csv('~/Downloads/driver_genes_martincorena_2017.txt', stringsAsFactors = F) %>% unlist()
names(im) = NULL

mt = read.csv('~/Downloads/driver_genes_tarabichi_2018.txt', stringsAsFactors = F) %>% unlist()
names(mt) = NULL

w = read.csv('~/Downloads/wangK562essentialgenes.txt', stringsAsFactors = F) %>% unlist()
names(w) = NULL

w = w[w != 'AK6']

b = read.csv('~/Downloads/bloomenKBM7essentials.txt', stringsAsFactors = F) %>% unlist()
names(b) = NULL

b = b[!(b %in% c('C22orf28', 'CENPC1', 'FTSJD2', 'GTPBP5', 'LRRC33', 'MKI67IP', 'MLL', 'MNF1', 'SELRC1', 'TMEM48', 'UQCC'))]

data('cancer_genes_dnds', package = 'mobster')

cancer_genes_dnds = NULL
cancer_genes_dnds$Martincorena_drivers = im
cancer_genes_dnds$Tarabichi_drivers = mt
cancer_genes_dnds$Wang_essentials = w
cancer_genes_dnds$Bloomen_essentials = b

usethis::use_data(cancer_genes_dnds, overwrite = T)


PD4120a_breast_sample$best$data$chr = paste0('chr', PD4120a_breast_sample$best$data$chr)
PD4120a_breast_sample$runs = lapply(PD4120a_breast_sample$runs, function(x){ 
  x$data$chr = paste0('chr', x$data$chr)
  x
  } )

LUFF76_lung_sample$best$data$chr = paste0('chr', LUFF76_lung_sample$best$data$chr)
LUFF76_lung_sample$runs = lapply(LUFF76_lung_sample$runs, function(x){ 
  x$data$chr = paste0('chr', x$data$chr)
  x
} )

LU4_lung_sample$best$data$chr = paste0('chr', LU4_lung_sample$best$data$chr)
LU4_lung_sample$runs = lapply(LU4_lung_sample$runs, function(x){ 
  x$data$chr = paste0('chr', x$data$chr)
  x
} )

usethis::use_data(LU4_lung_sample, overwrite = T)
usethis::use_data(LUFF76_lung_sample, overwrite = T)
usethis::use_data(PD4120a_breast_sample, overwrite = T)


