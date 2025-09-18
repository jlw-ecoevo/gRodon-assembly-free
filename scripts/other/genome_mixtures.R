#JN and JLW 2025

# Load Packages ----------------------------------------------------------------

library(dplyr)
library(gRodon)
library(parallel)
library(Biostrings)
library(coRdon)
library(matrixStats)


# Calculate Codon Usage Statistics ---------------------------------------------


setwd("/gpfs/projects/WeissmanGroup/gRodon_training/madin_genomes_only/")
mix_df <- read.delim("../genome_mixtures_broad2d.tsv",sep=" ")


statsMixture <- function(i,mix_df,bg="all"){
  print(paste("mixture: ",i))
  print(mix_df$G1)
  print(mix_df$G10)
  
  G1 <- readDNAStringSet(mix_df$G1[i])
  G2 <- readDNAStringSet(mix_df$G2[i])
  G3 <- readDNAStringSet(mix_df$G3[i])
  G4 <- readDNAStringSet(mix_df$G4[i])
  G5 <- readDNAStringSet(mix_df$G5[i])
  G6 <- readDNAStringSet(mix_df$G6[i])
  G7 <- readDNAStringSet(mix_df$G7[i])
  G8 <- readDNAStringSet(mix_df$G8[i])
  G9 <- readDNAStringSet(mix_df$G9[i])
  G10 <- readDNAStringSet(mix_df$G10[i])
  
  genes <- c(G1,G2,G3,G4,G5,G6,G7,G8,G9,G10)
  
  cu150 <- gRodon:::getStatistics(gene_file=NA,
                                  genes=genes,
                                  fragments = NA,
                                  bg = bg,
                                  trimlen = 150,
                                  trimside = "start",
                                  all_metrics = T)
  cu <- gRodon:::getStatistics(gene_file=NA,
                               genes=genes,
                               fragments = NA,
                               bg = bg,
                               trimlen = 510,
                               trimside = "start",
                               all_metrics = T)
  highly_expressed <- grepl("^(?!.*(methyl|hydroxy)).*0S ribosomal protein",names(genes),ignore.case = T, perl = TRUE)
  growth <- predictGrowth(genes,
                          highly_expressed,
                          mode="metagenome_v2")
  
  print("cu:")
  print(cu)
  print("cu150:")
  print(cu150)
  print("growth")
  print(growth)
  
  mix_df_i <- data.frame(SIM=mix_df$SIM[i],
                         G1=mix_df$G1[i],
                         G2=mix_df$G2[i],
                         G3=mix_df$G3[i],
                         G4=mix_df$G4[i],
                         G5=mix_df$G5[i],
                         G6=mix_df$G6[i],
                         G7=mix_df$G7[i],
                         G8=mix_df$G8[i],
                         G9=mix_df$G9[i],
                         G10=mix_df$G10[i],
                         MILC=cu$MILC,
                         ENCprime=cu$ENCprime,
                         SCUO=cu$SCUO,
                         B=cu$B,
                         MCB=cu$MCB,
                         MILC.150=cu150$MILC,
                         ENCprime.150=cu150$ENCprime,
                         SCUO.150=cu150$SCUO,
                         B.150=cu150$B,
                         MCB.150=cu150$MCB,
                         CUBHE=growth$CUBHE,
                         GC=growth$GC,
                         GCdiv=growth$GCdiv,
                         CUB=growth$CUB,
                         nHE=growth$nHE,
                         dCUB=growth$dCUB,
                         d.mmv2=growth$d,
                         Consistency=growth$ConsistencyHE)
  
  return(mix_df_i)
}


lapply(1:nrow(mix_df),statsMixture,mix_df=mix_df,bg="individual")


mixture_list <- mclapply(1:nrow(mix_df),
                         statsMixture,
                         mix_df=mix_df,
                         bg="individual",
                         mc.cores=28)

mixture_df <- do.call("rbind",mixture_list) %>%
  as.data.frame(stringsAsFactors=F)

setwd("/gpfs/projects/WeissmanGroup/gRodon_training/")
save(mixture_df, file = "CodonStatistics_mixtures_broad2.rda")