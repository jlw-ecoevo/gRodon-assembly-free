#JN and JLW 2025

# Load Packages ----------------------------------------------------------------

library(dplyr)
library(gRodon)
library(parallel)
library(Biostrings)
library(coRdon)
library(matrixStats)


# Calculate Codon Usage Statistics ---------------------------------------------


setwd("/gpfs/projects/WeissmanGroup/gRodon_training/madin_genomes/")
sim_names <- readLines("../sim.names")


statsMeta <- function(simname,bg="all"){
  genes <- readDNAStringSet(paste0(simname,"_fraggenescan.ffn"))
  ribo_names <- readLines(paste0(simname,".genes"))
  highly_expressed <- names(genes) %in% ribo_names
  
  
  growth150 <- try(predictGrowth(genes,
                                 highly_expressed,
                                 mode="metagenome_150bp"))
  growth <- try(predictGrowth(genes,
                              highly_expressed,
                              mode="metagenome_v2"))
  
  if(inherits(growth150,"try-error") | inherits(growth,"try-error")){
    mix_df_i <- data.frame(SIM=simname,
                           MILC=NA,
                           ENCprime=NA,
                           SCUO=NA,
                           B=NA,
                           MCB=NA,
                           CUBHE=NA,
                           GC=NA,
                           GCdiv=NA,
                           CUB=NA,
                           nHE=NA,
                           dCUB=NA,
                           d.150=NA,
                           d.mmv2=NA,
                           Consistency=NA)
    return(mix_df_i)
  }   else {
    
    mix_df_i <- data.frame(SIM=simname,
                           MILC=growth150$MILC,
                           ENCprime=growth150$ENCprime,
                           SCUO=growth150$SCUO,
                           B=growth150$B,
                           MCB=growth150$MCB,
                           CUBHE=growth$CUBHE,
                           GC=growth$GC,
                           GCdiv=growth$GCdiv,
                           CUB=growth$CUB,
                           nHE=growth$nHE,
                           dCUB=growth$dCUB,
                           d.150=growth150$d,
                           d.mmv2=growth$d,
                           Consistency=growth$ConsistencyHE)
    
    return(mix_df_i)}
}



mixture_list <- mclapply(sim_names,
                         statsMeta,
                         bg="individual",
                         mc.cores=28)


sim_df <- do.call("rbind",mixture_list) %>%
  as.data.frame(stringsAsFactors=F)

setwd("/gpfs/projects/WeissmanGroup/gRodon_training/")
save(sim_df, file = "CodonStatistics_simmeta150.rda")