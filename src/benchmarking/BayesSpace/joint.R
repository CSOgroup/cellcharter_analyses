suppressPackageStartupMessages({
	library(BayesSpace)
	library(mclust)
	library(harmony)
})

args <- commandArgs(trailingOnly=T);
if (length(args) > 0) {
    harmony_key = args[1]
} else {
    harmony_key = NULL
}



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

groups = c(
	"151507"=0, 
    "151508"= 0, 
    "151509"= 0, 
    "151510"= 0, 
    "151669"= 1, 
    "151670"= 1, 
    "151671"= 1, 
    "151672"= 1, 
    "151673"= 2, 
    "151674"= 2, 
    "151675"= 2, 
    "151676"= 2
)

nrep = 2000
hvg = 1000
pc = 15
gamma = 3

ARIs = list()

sample_list = c()
for (name in names(samples)) {
    s = readRDS(paste('../../../data/Visium_DLPFC/preprocessed_rds/', name, '.rds', sep=""))
    s$group = groups[name]
    rowData(s)$is.HVG = NULL 
    sample_list = c(sample_list, s)
}

sample = do.call(cbind, sample_list)


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

dec <- scran::modelGeneVar(sample)
top <- scran::getTopHVGs(dec, n = hvg)

sample <- scater::runPCA(sample, subset_row=top)

sample <- spatialPreprocess(sample, platform="Visium", skip.PCA=TRUE)


labels = data.frame(list('layer_guess'=sample$layer_guess, 'sample'=sample$sample_name))

set.seed(12345)
seeds = sample.int(32768, 10)

harmony_text = ''
if (!is.null(harmony_key)) {
	harmony_text = paste('_harmony-', harmony_key, sep="")
}

for(i in seeds) {
	ARIs_sample = c()
    set.seed(i)

    if (!is.null(harmony_key)) {
		sample <- RunHarmony(sample, harmony_key, verbose = F)
		sample <- spatialCluster(sample, q=7, d=pc, platform='Visium', use.dimred = "HARMONY",
                                                nrep=nrep, gamma=gamma, save.chain=TRUE)
	} else {
		sample <- spatialCluster(sample, q=7, d=pc, platform='Visium',
                                                nrep=nrep, gamma=gamma, save.chain=TRUE)
	}

    labels[[paste('cluster_', i, sep="")]] = sample$spatial.cluster
    write.csv(labels, paste('../../../results/benchmarking/joint/labels_BayesSpace_hvg', hvg, '_pc', pc, '_gamma', gamma, '_nrep', nrep, harmony_text, '_joint.csv', sep=""))

    for (name in names(samples)) {
        sample_single = sample[, sample$sample_name == name]
        ari = adjustedRandIndex(
            sample_single[, !is.na(colData(sample_single)$layer_guess)]$layer_guess,
            sample_single[, !is.na(colData(sample_single)$layer_guess)]$spatial.cluster
        )
        ARIs_sample[[name]] = ari
    }
	ARIs = c(ARIs, ARIs_sample)
	write.csv(as.data.frame(do.call(rbind, ARIs)), paste('../../../results/benchmarking/joint/ARI_BayesSpace_hvg', hvg, '_pc', pc, '_gamma', gamma, '_nrep', nrep, harmony_text, '_joint.csv', sep=""))
}

