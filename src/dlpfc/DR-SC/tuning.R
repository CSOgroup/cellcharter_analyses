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
#data("dlpfc151510", package = 'DR.SC')
#dlpfc151510 <- NormalizeData(dlpfc151510, verbose = F)
#seu <- FindVariableFeatures(dlpfc151510, nfeatures = 500, verbose = F)


#InstallData('ssHippo')
#data('ssHippo')
#locations = GetTissueCoordinates(ssHippo)[c('x', 'y')]
#ssHippo@meta.data = locations
#names(ssHippo@meta.data) = c('row', 'col')

#ssHippo <- NormalizeData(ssHippo, verbose = F)
#seuHippo <- FindVariableFeatures(ssHippo, nfeatures = 500, verbose = F)

# balbc = readRDS('/scratch/mvarrone/BALBc-3.rds')
# seu = CreateSeuratObject(balbc@raw_exprs, meta.data=balbc@spatial_locs)
# names(seu@meta.data) = c('orig.ident', 'nCount_RNA', 'nFeature_RNA', 'row', 'col', 'cell_ID')
# seu <- FindVariableFeatures(seu, nfeatures = 30, verbose = F)
# seu <- DR.SC(seu, K=7, platform = 'slide-seqv2', verbose=T)
# saveRDS(seu, '/scratch/mvarrone/BALBc-3_DRSC.rds')
# write.csv(seu@meta.data, '/scratch/mvarrone/BALBc-3_metadata.csv')
# write.csv(Idents(seu), '/scratch/mvarrone/BALBc-3_labels.csv')

ARIs = list()
for (name in names(samples)) {
    for (filter in filters) {
        for (hvg in hvgs) {
            ARIs_sample = c()
            
            for(i in 1:5) {
                sample = as.Seurat(readRDS(paste('/scratch/mvarrone/', name, '.rds', sep="")))

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
                write.csv(t(as.data.frame(do.call(cbind, ARIs))), paste('/work/FAC/FBM/DBC/gciriell/spacegene/Packages/cellcharter_analyses/results/dlpfc/DR-SC/ARI_tuning_', name, '.csv', sep=""))
            }
        }
    }
    
}
