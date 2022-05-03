suppressPackageStartupMessages({
	library(BayesSpace)
	library(mclust)
	library(harmony)
})
    

args = commandArgs(trailingOnly=TRUE)

n_samples = args[1]
if (length(args) > 1) {
    harmony_key = args[2]
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

run = function(N_SAMPLES) {
	sample_list = c()
	for (name in names(samples)[1:N_SAMPLES]) {
        s = readRDS(paste('/scratch/mvarrone/', name, '.rds', sep=""))
		s$group = groups[name]
        rowData(s)$is.HVG = NULL 
		sample_list = c(sample_list, s)
	}

	if (N_SAMPLES == 1) {
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

	dec <- scran::modelGeneVar(sample)
    top <- scran::getTopHVGs(dec, n = hvg)

	sample <- scater::runPCA(sample, subset_row=top)

	sample <- spatialPreprocess(sample, platform="Visium", skip.PCA=TRUE)

	if (!is.null(harmony_key)) {
		sample <- RunHarmony(sample, harmony_key, verbose = F)
		sample <- spatialCluster(sample, q=7, d=pc, platform='Visium', use.dimred = "HARMONY",
                                                nrep=nrep, gamma=gamma, save.chain=TRUE)
	} else {
		sample <- spatialCluster(sample, q=7, d=pc, platform='Visium',
                                                nrep=nrep, gamma=gamma, save.chain=TRUE)
	}

    
}

gc(reset=TRUE)

run(n_samples)

    
memory_info = gc()
max_mem = memory_info[12]

row = data.frame(max_mem, row.names=c(n_samples))

harmony_text = ''
if (!is.null(harmony_key)) {
	harmony_text = paste('_harmony-', harmony_key, sep="")
}
mem_path = paste('/work/FAC/FBM/DBC/gciriell/spacegene/Packages/cellcharter_analyses/results/dlpfc/BayesSpace/memory/memory_hvg', hvg, '_pc', pc, '_gamma', gamma, '_nrep', nrep, harmony_text, '.csv', sep='')


if (!file.exists(mem_path)) {
	mems = row
} else {
	mems = read.csv(mem_path, row.names=1)
	mems = rbind(mems, row)
}

write.csv(mems, mem_path)
