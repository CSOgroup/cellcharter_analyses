library('DR.SC')
library('Seurat')
library('SeuratData')
suppressPackageStartupMessages(library(mclust))

samples = c(
    "151507"=7, 
    "151508"= 7, 
    "151509"= 7, 
    "151510"= 7, 
    "151669"= 5, 
    "151670"= 5, 
    "151671"= 5, 
    "151672"= 5, 
    "151673"= 7, 
    "151674"= 7, 
    "151675"= 7, 
    "151676"= 7
)
hvg = 2000
filter =  'svg'

set.seed(12345)
seeds = sample.int(32768, 10)

ARIs = list()
for (name in names(samples)) {
    ARIs_sample = c()
    
    for(i in seeds) {
        sample = as.Seurat(readRDS(paste('../../../data/Visium_DLPFC/preprocessed_rds/', name, '.rds', sep="")))
        sample = NormalizeData(sample, verbose = F)

        if (filter == 'hvg') {
            sample <- FindVariableFeatures(sample, nfeatures = hvg, verbose = F)
        }
        else {
            sample <- FindSVGs(sample, nfeatures = hvg)
        }
        
        set.seed(i)

        sample <- DR.SC(sample, K=samples[name], platform = 'Visium', verbose=T)
        ari = adjustedRandIndex(
            Idents(sample)[!is.na(sample@meta.data$layer_guess)], 
            sample[, !is.na(sample@meta.data$layer_guess)]@meta.data$layer_guess
        )
        ARIs_sample = c(ARIs_sample, ari)
        ARIs[[name]] = ARIs_sample
        write.csv(t(as.data.frame(do.call(cbind, ARIs))), paste('../../../results/benchmarking/individual/ARI_DR-SC_', filter, hvg, '_individual.csv', sep=""))
    }
}
