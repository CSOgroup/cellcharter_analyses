suppressPackageStartupMessages({
	library(BayesSpace)
	library(mclust)
	library(harmony)
})

args = commandArgs(trailingOnly=TRUE)

n_samples = args[1]
n_cpus = args[2]
if (length(args) > 2) {
    harmony_key = args[3]
} else {
    harmony_key = NULL
}

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

sample_list = c()
for (name in names(samples)[1:n_samples]) {
	s = readRDS(paste('../../../data/Visium_DLPFC/preprocessed_rds/', name, '.rds', sep=""))
	s$group = groups[name]
	rowData(s)$is.HVG = NULL 
	sample_list = c(sample_list, s)
}

if (n_samples == 1) {
	sample = sample_list[[1]]
} else {
	sample = do.call(cbind, sample_list)
}

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

sample$row[sample$sample_name == "151507"] = 
450 + sample$row[sample$sample_name == "151507"]

sample$col[sample$sample_name == "151672"] = 
450 + sample$col[sample$sample_name == "151672"]

sample$row[sample$sample_name == "151673"] = 
450 + sample$row[sample$sample_name == "151673"]
sample$col[sample$sample_name == "151673"] = 
150 + sample$col[sample$sample_name == "151673"]

harmony_text = ''
if (!is.null(harmony_key)) {
	harmony_text = paste('_harmony-', harmony_key, sep="")
}
time_path = paste('../../../results/benchmarking/time/time_BayesSpace_hvg', hvg, '_pc', pc, '_gamma', gamma, '_nrep', nrep, harmony_text, '_ncpus', n_cpus,'.csv', sep='')

set.seed(12345)
seeds = sample.int(32768, 10)

for (i in seeds) {
    set.seed(i)
    start_hvg = Sys.time()
	dec <- scran::modelGeneVar(sample)
	top <- scran::getTopHVGs(dec, n = hvg)
	time_hvg = difftime(Sys.time(), start_hvg, units = "secs")

	start_pca = Sys.time()
	sample <- scater::runPCA(sample, subset_row=top)
	time_pca = difftime(Sys.time(), start_pca, units = "secs")

	start_preprocess = Sys.time()
	sample <- spatialPreprocess(sample, platform="Visium", skip.PCA=TRUE)
	time_preprocess = difftime(Sys.time(), start_preprocess, units = "secs")
	
	time_harmony = 0
	if (!is.null(harmony_key)) {
		start_harmony = Sys.time()
		sample <- RunHarmony(sample, harmony_key, verbose = F)
		time_harmony = difftime(Sys.time(), start_harmony, units = "secs")

		start_cluster = Sys.time()
		sample <- spatialCluster(sample, q=7, d=pc, platform='Visium', use.dimred = "HARMONY",
                                                nrep=nrep, gamma=gamma, save.chain=TRUE)
		time_cluster = difftime(Sys.time(), start_cluster, units = "secs")
	} else {
		start_cluster = Sys.time()
		sample <- spatialCluster(sample, q=7, d=pc, platform='Visium',
                                                nrep=nrep, gamma=gamma, save.chain=TRUE)
		time_cluster = difftime(Sys.time(), start_cluster, units = "secs")
	}
	
	row = data.frame(list('n_samples'=n_samples, 'HVGs'=time_hvg, 'PCA'=time_pca, 'preprocess'=time_preprocess, 'harmony'=time_harmony, 'clustering'=time_cluster))
	if (!file.exists(time_path)) {
		times = row
	} else {
		times = read.csv(time_path, row.names=1)
		print(times)
		times = rbind(times, row)
	}
	write.csv(times, time_path)
}
