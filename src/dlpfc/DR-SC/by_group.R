library('DR.SC')
library('Seurat')
library('SeuratData')
library(harmony)
suppressPackageStartupMessages(library(mclust))

samples = c(
    #"151507"=7, 
    "151508"= 7, 
    "151509"= 7, 
    "151510"= 7, 
    "151669"= 5, 
    "151670"= 5, 
    "151671"= 5, 
    #"151672"= 5, 
    #"151673"= 7, 
    "151674"= 7, 
    "151675"= 7, 
    "151676"= 7
)
hvg = 2000
filter =  'svg'

ARIs = list()
sample_list = c()

for (name in names(samples)) {
    sample_list = c(sample_list, as.Seurat(readRDS(paste('/scratch/mvarrone/', name, '.rds', sep=""))))
}

set.seed(12345)
seeds = sample.int(32768, 10)
for (j in 0:2) {
    sample = merge(sample_list[[(j*3+1)]], y = sample_list[(j*3+2):(j*3+3)], add.cell.ids = names(samples)[(j*3+1):(j*3+3)])
    sample$row[sample$sample_name == "151509"] = 
    150 + sample$row[sample$sample_name == "151509"]

    sample$col[sample$sample_name == "151510"] = 
    150 + sample$col[sample$sample_name == "151510"]

    sample$row[sample$sample_name == "151669"] = 
    150 + sample$row[sample$sample_name == "151669"]
    sample$col[sample$sample_name == "151669"] = 
    150 + sample$col[sample$sample_name == "151669"]

    sample$row[sample$sample_name == "151670"] = 
    300 + sample$row[sample$sample_name == "151670"]

    sample$col[sample$sample_name == "151671"] = 
    300 + sample$col[sample$sample_name == "151671"]

    sample$row[sample$sample_name == "151674"] = 
    300 + sample$row[sample$sample_name == "151674"]
    sample$col[sample$sample_name == "151674"] = 
    150 + sample$col[sample$sample_name == "151674"]


    sample$row[sample$sample_name == "151675"] = 
    150 + sample$row[sample$sample_name == "151675"]
    sample$col[sample$sample_name == "151675"] = 
    300 + sample$col[sample$sample_name == "151675"]

    sample$row[sample$sample_name == "151676"] = 
    300 + sample$row[sample$sample_name == "151676"]
    sample$col[sample$sample_name == "151676"] = 
    300 + sample$col[sample$sample_name == "151676"]

    if (filter == 'hvg') {
        sample <- FindVariableFeatures(sample, nfeatures = hvg, verbose = F)
    } else {
        sample <- FindSVGs(sample, nfeatures = hvg)
    }

    labels = data.frame(list('layer_guess'=sample@meta.data$layer_guess))

    for(i in seeds) {
        ARIs_sample = c()
        set.seed(i)

        sample <- DR.SC(sample, K=7, platform = 'Visium', verbose=T)

        labels[[paste('cluster_', i, sep="")]] = Idents(sample)
        write.csv(labels, paste('/work/FAC/FBM/DBC/gciriell/spacegene/Packages/cellcharter_analyses/results/dlpfc/DR-SC/labels/ARI_', filter, hvg, '_group_combined.csv', sep=""))

        for (name in names(samples)[(j*3+1):(j*3+3)]) {
            sample_single = sample[, sample$sample_name == name]
            ari = adjustedRandIndex(
                Idents(sample_single)[!is.na(sample_single@meta.data$layer_guess)], 
                sample_single[, !is.na(sample_single@meta.data$layer_guess)]@meta.data$layer_guess
            )
            ARIs_sample[[name]] = ari
        }
        ARIs = c(ARIs, ARIs_sample)
        write.csv(t(as.data.frame(do.call(cbind, ARIs))), paste('/work/FAC/FBM/DBC/gciriell/spacegene/Packages/cellcharter_analyses/results/dlpfc/DR-SC/accuracy/ARI_', filter, hvg, '_group_combined.csv', sep=""))
    }
}


