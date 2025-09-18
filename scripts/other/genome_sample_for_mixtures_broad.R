library(dplyr)
library(data.table)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}

rgrep <- function(big,small_vec){
  small_vec[lapply(small_vec,grepl,x=big) %>% unlist()]
}

setwd("~/../../Documents/gRodon-assembly-free/data")
load("GrowthRates_Madin.rda")
load("Accession2Species_Madin.rda")
load("CodonStatistics_Madin.rda")


cu <- cu %>% mutate_all(unlist)
cu$Accession <- cu$File %>% gsub(pattern="[.].*",replace="")
names(d)[1] <- "Species"
d <- d %>% as.data.frame(stringsAsFactors=F)
rownames(spp_acc) <- spp_acc$V1 %>% gsub(pattern="[.].*",replace="")
cu$Spp <- spp_acc[cu$Accession,"V2"]
cu$Species <- lapply(cu$Spp,rgrep,small_vec=d$Species) %>%
  lapply("[",1) %>% unlist()
cu$Species[cu$Spp %in% d$Species] <- cu$Spp[cu$Spp %in% d$Species]
cu <- merge.easy(cu,d,key="Species") %>% subset(!is.na(Species))
cu$OGT <- cu$OptTemp
cu$OGT[cu$OptTemp=="NaN"] <- cu$GrowthTemp[cu$OptTemp=="NaN"] 
cu <- cu %>% subset(!is.na(OGT))

x <- readLines("cds.files")

x <- x[gsub(pattern="[.].*",replace="",x) %in% cu$Accession[cu$OGT<42 & cu$OGT>16]]

nsim <- 100
df <- data.frame()
counter <- 1
cutoffs <- seq(-1,log10(5),0.03)
for(j in 1:length(cutoffs)){
  x_cut <- x[gsub(pattern="[.].*",replace="",x) %in% cu$Accession[cu$d>cutoffs[j]]]
  for(i in 1:nsim){
    samp <- sample(x_cut,10)
    df <- rbind(df,data.frame(
      SIM=paste0("SIM",counter),
      G1=samp[1],
      G2=samp[2],
      G3=samp[3],
      G4=samp[4],
      G5=samp[5],
      G6=samp[6],
      G7=samp[7],
      G8=samp[8],
      G9=samp[9],
      G10=samp[10]))
    counter <- counter+1
  }
}

write.table(df,
            file = "genome_mixtures_broad2.tsv",
            row.names = F,
            col.names = T,
            quote = F)
