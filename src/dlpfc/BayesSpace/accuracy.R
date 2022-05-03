suppressPackageStartupMessages(library(BayesSpace))
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

nrep = 2000
hvg = 1000
pc = 15
gamma = 3

set.seed(12345)
seeds = sample.int(32768, 10)

ARIs = list()
for (name in names(samples)) {
    ARIs_sample = c()

    sample = readRDS(paste('/scratch/mvarrone/', name, '.rds', sep=""))
    dec <- scran::modelGeneVar(sample)
    top <- scran::getTopHVGs(dec, n = hvg)

    sample <- scater::runPCA(sample, subset_row=top)

    sample <- spatialPreprocess(sample, platform="Visium", skip.PCA=TRUE)
    
    for(i in seeds) {
        set.seed(i)

        sample <- spatialCluster(sample, q=samples[name], d=pc, platform='Visium',
                                                nrep=nrep, gamma=gamma, save.chain=TRUE)
        ari = adjustedRandIndex(
            sample[, !is.na(colData(sample)$layer_guess)]$layer_guess,
            sample[, !is.na(colData(sample)$layer_guess)]$spatial.cluster
        )
        ARIs_sample = c(ARIs_sample, ari)
        ARIs[[name]] = ARIs_sample
        write.csv(as.data.frame(do.call(rbind, ARIs)), paste('/work/FAC/FBM/DBC/gciriell/spacegene/Packages/cellcharter_analyses/results/dlpfc/BayesSpace/accuracy/ARI_hvg', hvg, '_pc', pc, '_gamma', gamma, '_nrep', nrep, '.csv', sep=""))
    }
}
