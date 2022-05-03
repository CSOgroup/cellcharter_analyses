suppressPackageStartupMessages(library(BayesSpace))
suppressPackageStartupMessages(library(mclust))

samples = c(
    "151507"= 7#, 
    #"151672"= 5#, 
    #"151673"= 7
)

hvgs = c(
    #500, 
    #1000#, 
    #2000#, 
    5000
)
pcs = c(
    5, 
    7, 
    10, 
    15, 
    20, 
    25, 
    30)
nreps = c(
    #2000#, 
    5000
    )
gammas = c(1,2,3,4)

ARIs = list()
times = list()

for (name in names(samples)) {
    
    sample = readRDS(paste('/scratch/mvarrone/', name, '.rds', sep=""))
    dec <- scran::modelGeneVar(sample)

    for (nrep in nreps) {
        for (hvg in hvgs) {
            top <- scran::getTopHVGs(dec, n = hvg)
            for (pc in pcs) {
                for (gamma in gammas) {
                    ARIs_sample = c()
                    for(i in 1:5) {
                        set.seed(i)
                        scater::runPCA(sample, subset_row=top)
                        sample <- spatialPreprocess(sample, platform="Visium", skip.PCA=TRUE)

                        q <- samples[name]  # Number of clusters
                        d <- pc  # Number of PCs

                        sample <- spatialCluster(sample, q=q, d=d, platform='Visium',
                                                nrep=nrep, gamma=gamma, save.chain=TRUE)
                        

                        ari = adjustedRandIndex(
                            sample[, !is.na(colData(sample)$layer_guess)]$layer_guess, 
                            sample[, !is.na(colData(sample)$layer_guess)]$spatial.cluster
                        )

                        ARIs_sample = c(ARIs_sample, ari)
                    }
                    ARIs[[paste(name, nrep, hvg, pc, gamma)]] = ARIs_sample
                    write.csv(t(as.data.frame(do.call(cbind, ARIs))), paste('/work/FAC/FBM/DBC/gciriell/spacegene/Packages/cellcharter_analyses/results/dlpfc/BayesSpace/tuning/ARI_tuning_', name, '.csv', sep=""))
                }
            }    
        }    
    }
}