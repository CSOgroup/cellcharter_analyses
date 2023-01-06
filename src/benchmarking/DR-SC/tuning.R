library('DR.SC')
library('Seurat')
library('SeuratData')
suppressPackageStartupMessages(library(mclust))


samples = c(
    "151507"= 7, 
    "151672"= 5, 
    "151673"= 7
)

hvgs = c(
    500, 
    1000, 
    2000, 
    5000
)

filters = c('hvg', 'svg')

ARIs = list()
for (name in names(samples)) {
    for (filter in filters) {
        for (hvg in hvgs) {
            ARIs_sample = c()
            
            for(i in 1:5) {
                sample = as.Seurat(readRDS(paste('../../../data/Visium_DLPFC/preprocessed_rds/', name, '.rds', sep="")))

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
                ARIs[[paste(name, hvg, filter)]] = ARIs_sample
                write.csv(t(as.data.frame(do.call(cbind, ARIs))), paste('../../../results/benchmarking/hyperparm_tuning/ARI_DR-SC_tuning_', name, '.csv', sep=""))
            }
        }
    }
    
}
