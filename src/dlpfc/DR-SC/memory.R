library(peakRAM)
library('DR.SC')
library('Seurat')
library('SeuratData')
library(harmony)
suppressPackageStartupMessages(library(mclust))

args = commandArgs(trailingOnly=TRUE)

n_samples = args[1]


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

run = function(N_SAMPLES) {
	sample_list = c()
	for (name in names(samples)[1:N_SAMPLES]) {
		sample_list = c(sample_list, as.Seurat(readRDS(paste('/scratch/mvarrone/', name, '.rds', sep=""))))
	}

	if (N_SAMPLES == 1) {
		sample = sample_list[[1]]
	} else {
		sample = merge(sample_list[[1]], y = sample_list[2:length(sample_list)], add.cell.ids = names(samples)[1:N_SAMPLES])
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

	if (filter == 'hvg') {
		sample <- FindVariableFeatures(sample, nfeatures = hvg, verbose = F)
	} else {
		sample <- FindSVGs(sample, nfeatures = hvg)
	}
	sample <- DR.SC(sample, K=7, platform = 'Visium', verbose=T)
}

gc(reset=TRUE)

run(n_samples)

    
memory_info = gc()
max_mem = memory_info[12]

row = data.frame(max_mem, row.names=c(n_samples))


mem_path = paste('/work/FAC/FBM/DBC/gciriell/spacegene/Packages/cellcharter_analyses/results/dlpfc/DR-SC/memory/memory_', filter, hvg, '.csv', sep='')


if (!file.exists(mem_path)) {
	mems = row
} else {
	mems = read.csv(mem_path, row.names=1)
	mems = rbind(mems, row)
}

write.csv(mems, mem_path)
