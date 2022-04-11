library('DR.SC')
library('Seurat')
library('SeuratData')
suppressPackageStartupMessages(library(mclust))

args = commandArgs(trailingOnly=TRUE)

n_samples = args[1]
n_cpus = args[2]

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

sample_list = c()
for (name in names(samples)[1:n_samples]) {
    sample_list = c(sample_list, as.Seurat(readRDS(paste('/scratch/mvarrone/', name, '.rds', sep=""))))
}

if (n_samples == 1) {
    sample = sample_list[[1]]
} else {
    sample = merge(sample_list[[1]], y = sample_list[2:length(sample_list)], add.cell.ids = names(samples)[1:n_samples])
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

set.seed(12345)
seeds = sample.int(32768, 10)
time_path = paste('/work/FAC/FBM/DBC/gciriell/spacegene/Packages/cellcharter_analyses/results/dlpfc/DR-SC/time/time_', filter, hvg, '_cpu_ncpus', n_cpus,'.csv', sep='')
for (i in seeds) {
    set.seed(i)
    start_hvg = Sys.time()
	if (filter == 'hvg') {
		sample <- FindVariableFeatures(sample, nfeatures = hvg, verbose = F)
	} else {
		sample <- FindSVGs(sample, nfeatures = hvg)
	}
	time_hvg = difftime(Sys.time(), start_hvg, units = "secs")

	start_cluster = Sys.time()
	sample <- DR.SC(sample, K=7, platform = 'Visium', verbose=T)
	time_cluster = difftime(Sys.time(), start_cluster, units = "secs")
	
	row = data.frame(list('n_samples'=n_samples, 'HVGs'=time_hvg, 'clustering'=time_cluster))
	if (!file.exists(time_path)) {
		times = row
	} else {
		times = read.csv(time_path, row.names=1)
		print(times)
		times = rbind(times, row)
	}
	write.csv(times, time_path)
}
