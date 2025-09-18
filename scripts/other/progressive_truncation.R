# Progressive truncation CUB
# JN and JLW 2025

# Load Packages ----------------------------------------------------------------

library(dplyr)
library(gRodon)
library(parallel)
library(Biostrings)
library(coRdon)
library(matrixStats)


# Calculate Codon Usage Statistics ---------------------------------------------

genomes_path <- "/gpfs/projects/WeissmanGroup/gRodon_training/madin_genomes/"
trimlens <- seq(75,450,75)

cu <- list()
for(i in 1:length(trimlens)){
  print(paste("truncate forward: ",trimlens[i]))
  cu[[i]] <- gRodon:::getStatisticsBatch(genomes_path,
                                         mc.cores = 28,
                                         fragments = NA,
                                         bg = "individual",
                                         trimlen = trimlens[i],
                                         trimside = "start",
                                         all_metrics = T)
}

cu_start <- list(cu=cu,trimlens=trimlens)

setwd("/gpfs/projects/WeissmanGroup/gRodon_training/")
save(cu_start, file = "CodonStatistics_truncated_genes_start_i.rda")

genomes_path <- "/gpfs/projects/WeissmanGroup/gRodon_training/madin_genomes/"
trimlens <- seq(75,450,75)

cu <- list()
for(i in 1:length(trimlens)){
  print(paste("truncate backward: ",trimlens[i]))
  cu[[i]] <- gRodon:::getStatisticsBatch(genomes_path,
                                         mc.cores = 28,
                                         fragments = NA,
                                         bg = "individual",
                                         trimlen = trimlens[i],
                                         trimside = "end",
                                         all_metrics = T)
}

cu_end <- list(cu=cu,trimlens=trimlens)

setwd("/gpfs/projects/WeissmanGroup/gRodon_training/")
save(cu_end, file = "CodonStatistics_truncated_genes_end_i.rda")