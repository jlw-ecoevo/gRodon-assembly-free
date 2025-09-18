
# JLW and JN 2005 --------------------------------------------------------------
#Benchmarking for gRodon modes metagenome_150bp and metagenome_250bp


# Packages and Helper Functions ------------------------------------------------

library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(caret)
library(pROC)
library(ggfortify)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}

boxcoxTransform <- function(x, lambda, back_transform = F) {
  if (back_transform == TRUE) {
    (x*lambda +1)^(1/lambda)  %>% return()
  } else {
    (((x^lambda) - 1) / lambda) %>% return()
  }
}

rgrep <- function(big,small_vec){
  small_vec[lapply(small_vec,grepl,x=big) %>% unlist()]
}

# CUB Metrics Under Truncation -------------------------------------------------

setwd("~/../../Documents/gRodon-assembly-free/data")
load("CodonStatistics_truncated_genes_start_i.rda")

pullCUB <- function(cu_start,metric){
  i <- 1
  cu <- cu_start$cu[[i]][,c("File",metric)] %>%
    subset(!grepl("Error",File))
  col_new <- paste0("CUBHE.",i)
  cu[,c(col_new)] <- unlist(cu[,c(metric)])
  cu$File = unlist(cu$File)
  cu <- cu[,-c(2)]
  for(i in 2:length(cu_start$trimlens)){
    if(nrow(cu_start$cu[[i]])>0){
      cu_i <- cu_start$cu[[i]][,c("File",metric)] %>%
        subset(!grepl("Error",File))
      col_new <- paste0("CUBHE.",i)
      cu_i[,c(col_new)] <- unlist(cu_i[,c(metric)])
      cu_i$File = unlist(cu_i$File)
      cu_i <- cu_i[,-c(2)]
      cu <- merge.easy(cu,cu_i,key="File")
    }
  }
  return(cu)
}

meltCUB <- function(cu,trimlens){
  n <- nrow(trimlens)
  cu_melt <- cbind(cu$File,100*(cu[,2:(n+1)] - as.data.frame(cu)[,n+1])/as.data.frame(cu)[,n+1]) %>%
    melt(id.vars="V1")
  cu_melt <- merge.easy(cu_melt,trimlens,key="variable")
  return(cu_melt)
}

n <- length(cu_start$trimlens)
trimlens <- data.frame(variable=paste0("CUBHE.",1:n),
                       trimlens=cu_start$trimlens)

MILC <- pullCUB(cu_start,"MILC") %>% meltCUB(.,trimlens)
ENCprime <- pullCUB(cu_start,"ENCprime") %>% meltCUB(.,trimlens)
B <- pullCUB(cu_start,"B") %>% meltCUB(.,trimlens)
SCUO <- pullCUB(cu_start,"SCUO") %>% meltCUB(.,trimlens)
MCB <- pullCUB(cu_start,"MCB") %>% meltCUB(.,trimlens)

pMILC <- ggplot(MILC,aes(x=trimlens,y=value,group=V1)) +
  geom_line(alpha=0.005) +
  theme_pubclean() +
  ylab("% Change in MILC") +
  xlab("Length of Truncated Genes")
pENCprime <- ggplot(ENCprime,aes(x=trimlens,y=value,group=V1)) +
  geom_line(alpha=0.005) +
  theme_pubclean() +
  ylab("% Change in ENC'") +
  xlab("Length of Truncated Genes")
pB <- ggplot(B,aes(x=trimlens,y=value,group=V1)) +
  geom_line(alpha=0.005) +
  theme_pubclean() +
  ylab("% Change in B") +
  xlab("Length of Truncated Genes")
pSCUO <- ggplot(SCUO,aes(x=trimlens,y=value,group=V1)) +
  geom_line(alpha=0.005) +
  theme_pubclean() +
  ylab("% Change in SCUO") +
  xlab("Length of Truncated Genes")
pMCB <- ggplot(MCB,aes(x=trimlens,y=value,group=V1)) +
  geom_line(alpha=0.005) +
  theme_pubclean() +
  ylab("% Change in MCB") +
  xlab("Length of Truncated Genes")

setwd("~/../../Documents/gRodon-assembly-free/figures")
png(file="truncate_genes.png",width=9,height=8,units="in",res=600)
ggarrange(ggarrange(pMILC,pENCprime,nrow=2,
                    labels=c("(a)","(b)")),
          ggarrange(pB,pSCUO,pMCB,nrow=3,
                    labels=c("(c)","(d)","(e)")),
          ncol=2)
dev.off()


# Add Growth Data --------------------------------------------------------------

setwd("~/../../Documents/gRodon-assembly-free/data")
load("GrowthRates_Madin.rda")
load("Accession2Species_Madin.rda")
load("CodonStatistics_Madin.rda")

cub_list <- list()
for(i in 1:nrow(trimlens)){
  cub_col <- paste0("CUBHE.",i)
  cub_list[[i]] <-
    data.frame(File = as.data.frame(pullCUB(cu_start,"MILC"))[,"File"],
               # PC1 = as.data.frame(PCCUB(cu_start,pccub,1))[,cub_col],
               # PC2 = as.data.frame(PCCUB(cu_start,pccub,2))[,cub_col],
               # PC3 = as.data.frame(PCCUB(cu_start,pccub,3))[,cub_col],
               MILC = as.data.frame(pullCUB(cu_start,"MILC"))[,cub_col],
               ENCprime = as.data.frame(pullCUB(cu_start,"ENCprime"))[,cub_col],
               B = as.data.frame(pullCUB(cu_start,"B"))[,cub_col],
               SCUO = as.data.frame(pullCUB(cu_start,"SCUO"))[,cub_col],
               MCB = as.data.frame(pullCUB(cu_start,"MCB"))[,cub_col])
}


names(d)[1] <- "Species"
d <- d %>% as.data.frame(stringsAsFactors=F)

# Merge datasets
rownames(spp_acc) <- spp_acc$V1 %>% gsub(pattern="[.].*",replace="")
cu <- cu %>% mutate_all(unlist)
cu$Accession <- cu$File %>% gsub(pattern="[.].*",replace="")
cu$Spp <- spp_acc[cu$Accession,"V2"]
cu$Species <- lapply(cu$Spp,rgrep,small_vec=d$Species) %>%
  lapply("[",1) %>% unlist()
cu$Species[cu$Spp %in% d$Species] <- cu$Spp[cu$Spp %in% d$Species]
cu <- merge.easy(cu,d,key="Species") %>% subset(!is.na(Species))

stat_data_list <- list()
stat_data_extremo_list <- list()
for(i in 1:nrow(trimlens)){
  cu_i <- cub_list[[i]]
  cu_i$Accession <- cu_i$File %>% gsub(pattern="[.].*",replace="")
  cu_i <- merge.easy(cu,cu_i,key="Accession") %>%
    subset(!is.na(MILC))

  stat_data_list[[i]] <- cu_i %>%
    subset(Extremophile == FALSE) %>%
    group_by(Species) %>%
    summarise_all(mean,na.rm=T) %>%
    subset(!is.na(Species))

  stat_data_extremo_list[[i]] <- cu_i %>%
    group_by(Species) %>%
    summarise_all(mean,na.rm=T) %>%
    subset(!is.na(Species))
  stat_data_extremo_list[[i]]$OGT <- stat_data_extremo_list[[i]]$OptTemp
  stat_data_extremo_list[[i]]$OGT[is.na(stat_data_extremo_list[[i]]$OGT)] <-
    stat_data_extremo_list[[i]]$GrowthTemp[is.na(stat_data_extremo_list[[i]]$OGT)]
}


# Cross-Test Length Specific Models --------------------------------------------

lengthPairTestC <- function(cub_train_len_df,
                           cub_test_len_df,
                           k = 5,
                           pred = F){
  #Function breaks data into k folds
  #then trains on training length CUBHE leaving each fold out iteratively
  #and tests using the testing length CUBHE for each left out fold

  fold_ids <- createFolds(1:nrow(cub_train_len_df),k=k,list=F)
  MSE_vec <- numeric(k)
  pred_df <- data.frame(d=numeric(nrow(cub_train_len_df)),
                        d.pred=numeric(nrow(cub_train_len_df)))
  for(i in 1:k){

    df.train <- cub_train_len_df[fold_ids!=i,] %>% as.data.frame()
    df.test <- cub_test_len_df[fold_ids==i,] %>% as.data.frame()

    mod_i <- lm(log(d)~MILC+ENCprime+B+SCUO+MCB,
                 data=df.train)

    # print(summary(mod_i))
    pred_i <- predict(mod_i,df.test)
    MSE_vec[i] <- mean((pred_i-log(df.test$d))^2)
    pred_df$d[fold_ids==i] <- df.test$d
    pred_df$d.pred[fold_ids==i] <- exp(pred_i)
  }
  if(pred==F){
    return(mean(MSE_vec))
  } else {
    return(pred_df)
  }
}

crossTrainModelC <- function(stat_data_list, trimlens){
  #Trains models for different train/test set truncation length combinations
  #To make heatmap
  trimlens_combos <- expand.grid(1:length(trimlens),1:length(trimlens))
  MSE_df <- data.frame(trimlen.train = trimlens[trimlens_combos[,1]],
                       trimlen.test = trimlens[trimlens_combos[,2]],
                       ind.train = trimlens_combos[,1],
                       ind.test = trimlens_combos[,2],
                       MSE = numeric(nrow(trimlens_combos)))
  for(i in 1:nrow(MSE_df)){
    cub_train_len_df <- stat_data_list[[MSE_df$ind.train[i]]]
    cub_test_len_df <- stat_data_list[[MSE_df$ind.test[i]]]
    MSE_df$MSE[i] <- lengthPairTestC(cub_train_len_df,
                                    cub_test_len_df)
  }

  return(MSE_df)
}

crossPredictModelC <- function(stat_data_list, trimlens,testlen,trainlen){
  #Trains model on trainlen and tests performance on testlen
  ind.train <- which(trimlens==trainlen)
  ind.test <- which(trimlens==testlen)
  cub_train_len_df <- stat_data_list[[ind.train]]
  cub_test_len_df <- stat_data_list[[ind.test]]
  pred_df <- lengthPairTestC(cub_train_len_df,
                            cub_test_len_df,
                            pred=T)
}


crosspred.l <- crossPredictModelC(stat_data_list,
                                 trimlens$trimlens,
                                 testlen=450,
                                 trainlen=450)

crosspred.s.l <- crossPredictModelC(stat_data_list,
                                   trimlens$trimlens,
                                   testlen=150,
                                   trainlen=450)

crosspred.s.s <- crossPredictModelC(stat_data_list,
                                   trimlens$trimlens,
                                   testlen=150,
                                   trainlen=150)

pL1 <- ggplot(crosspred.l,aes(x=log(2)/d,y=log(2)/d.pred))+
  scale_x_log10(limits=c(0.001,4)) +
  scale_y_log10(limits=c(0.001,4)) +
  geom_point(size=2,alpha=0.75,pch=21,fill="gray50") +
  theme_pubclean()  +
  geom_abline(slope=1,intercept=0,lty=2,color="#357BA2FF",lwd=1.2) +
  annotate("text",label="X=Y",angle=52,x=0.002,y=0.003,size=3,color="#357BA2FF") +
  xlab("Actual Max. Growth Rate (1/Hr)") +
  ylab("Predicted Max. Growth Rate (1/Hr)") +
  # geom_smooth(method="gam",color="white",fill="black") +
  ggtitle("Train 450bp/Test 450bp")

pL2 <- ggplot(crosspred.s.l,aes(x=log(2)/d,y=log(2)/d.pred))+
  scale_x_log10(limits=c(0.001,4)) +
  scale_y_log10(limits=c(0.001,4)) +
  geom_point(size=2,alpha=0.75,pch=21,fill="gray50") +
  theme_pubclean()  +
  geom_abline(slope=1,intercept=0,lty=2,color="#357BA2FF",lwd=1.2) +
  annotate("text",label="X=Y",angle=52,x=0.002,y=0.003,size=3,color="#357BA2FF") +
  xlab("Actual Max. Growth Rate (1/Hr)") +
  ylab("Predicted Max. Growth Rate (1/Hr)") +
  # geom_smooth(method="gam",color="white",fill="black") +
  ggtitle("Train 450bp/Test 150bp")

pL3 <- ggplot(crosspred.s.s,aes(x=log(2)/d,y=log(2)/d.pred))+
  scale_x_log10(limits=c(0.001,4)) +
  scale_y_log10(limits=c(0.001,4)) +
  geom_point(size=2,alpha=0.75,pch=21,fill="gray50") +
  theme_pubclean()  +
  geom_abline(slope=1,intercept=0,lty=2,color="#357BA2FF",lwd=1.2) +
  annotate("text",label="X=Y",angle=52,x=0.002,y=0.003,size=3,color="#357BA2FF") +
  xlab("Actual Max. Growth Rate (1/Hr)") +
  ylab("Predicted Max. Growth Rate (1/Hr)") +
  # geom_smooth(method="gam",color="white",fill="black") +
  ggtitle("Train 150bp/Test 150bp")

ct_mse <- crossTrainModelC(stat_data_list,trimlens$trimlens)
pHeat <- ggplot(ct_mse,
       aes(y=trimlen.test,
           x=trimlen.train,
           fill=log10(MSE))) +
  geom_tile() +
  scale_fill_viridis_c(direction = -1)+
  geom_vline(xintercept=240,lty=2,col="black") +
  geom_hline(yintercept=240,lty=2,col="black") +
  geom_vline(xintercept=150,lty=2,col="red") +
  geom_hline(yintercept=150,lty=2,col="red")+
  xlab("Train Gene Length") +
  ylab("Test Gene Length") +
  labs(fill=expression("CV "*log[10]*"(MSE)")) +
  theme_bw()



pMarg <- ggplot(ct_mse %>% subset(trimlen.test %in% c(150,300,450)),
       aes(x=trimlen.train,
           y=MSE,
           fill=factor(trimlen.test),
           group=factor(trimlen.test))) +
  geom_line(lwd=1.5) +
  geom_line(aes(color=factor(trimlen.test)),lwd=1) +
  geom_point(size=5,pch=21) +
  scale_y_log10() +
  scale_color_viridis_d(option="B") +
  scale_fill_viridis_d(option="B") +
  # geom_smooth(aes(color=factor(trimlen.test)),fill="black",method='gam') +
  theme_pubclean() +
  labs(color="Test Gene Length", fill="Test Gene Length") +
  theme(legend.position=c(0.75,0.85),
        legend.background = element_rect(fill = NA)) +
  xlab("Train Gene Length") +
  ylab("Cross-Validated MSE")


setwd("~/../../Documents/gRodon-assembly-free/figures")
png(file="cross_train_linear.png",width=9,height=8,units="in",res=600)
ggarrange(ggarrange(pHeat,
                    pMarg,
                    ncol=2,
                    labels = c("(a)","(b)"),
                    widths=c(3,2),
                    hjust=0),
          ggarrange(pL1,pL2,pL3,
                    ncol=3,
                    labels=c("(c)","(d)","(e)")),
          nrow=2)
dev.off()


# Training Full Model ----------------------------------------------------------

#Trains on all data (no CV)
mod_list <- list()
for(i in 1:length(trimlens$trimlens)){
  mod_list[[i]] <- lm(log(d)~MILC+ENCprime+B+SCUO+MCB,
                      data=stat_data_list[[i]])
}


# Genome Mixtures --------------------------------------------------------------

setwd("~/../../Documents/gRodon-assembly-free/data")
load("GrowthRates_Madin.rda")
load("Accession2Species_Madin.rda")
load("CodonStatistics_Madin.rda")
load("CodonStatistics_mixtures_broad2.rda")

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

cu$Fast <- cu$d<1
mix_df <- mixture_df %>%
  reshape2::melt(id.vars=c("SIM","d.mmv2","CUBHE","GC","GCdiv","CUB","nHE","dCUB","Consistency",
                 "MILC","ENCprime","B","SCUO","MCB",
                 "MILC.150","ENCprime.150","B.150","SCUO.150","MCB.150")) %>%
  mutate(Accession=substr(value,1,13)) %>%
  merge.easy(.,cu,key="Accession") %>%
  mutate(mu=log(2)/d,
         mu.mmv2=log(2)/d.mmv2) %>%
  subset(select=c(MILC,ENCprime,B,SCUO,MCB,
                  MILC.150,ENCprime.150,B.150,SCUO.150,MCB.150,
                  d,d.mmv2,Fast,SIM,mu.mmv2,mu)) %>%
  mutate(d=log(d)) %>%
  group_by(SIM) %>%
  summarise_all(mean) %>%
  mutate(d=exp(d)) %>%
  as.data.frame()

mix_df_150 <- mix_df %>%
  subset(select=c(MILC.150,ENCprime.150,B.150,SCUO.150,MCB.150))
names(mix_df_150) <- c("MILC","ENCprime","B","SCUO","MCB")
pred.150 <- matrix(ncol=length(trimlens$trimlens)+1,
                   nrow=5700)
pred.150[,1] <- mix_df$d
err.150 <- numeric()
for(i in 1:length(trimlens$trimlens)){
  pred.150[,i+1] <- predict(mod_list[[i]],mix_df_150)
  err.150[i] <- sum((log(mix_df$d)-pred.150[,i+1])^2)
}
colnames(pred.150) <- c("d",paste0("d.",trimlens$trimlens))

pred.450 <- matrix(ncol=length(trimlens$trimlens)+1,
                   nrow=5700)
pred.450[,1] <- mix_df$d
err.450 <- numeric()
for(i in 1:length(trimlens$trimlens)){
  pred.450[,i+1] <- predict(mod_list[[i]],mix_df)
  err.450[i] <- sum((log(mix_df$d)-pred.450[,i+1])^2)
}
colnames(pred.450) <- c("d",paste0("d.",trimlens$trimlens))

marg_df <- data.frame(err.150=err.150,
                      err.450=err.450,
                      trimlens=trimlens$trimlens) %>%
  reshape2::melt(id.vars="trimlens")
pCM <- ggplot(marg_df,
       aes(x=trimlens,
           y=value,
           fill=variable,
           group=variable)) +
  geom_line(lwd=1.5) +
  geom_line(aes(color=variable),lwd=1) +
  geom_point(size=3,pch=21) +
  scale_y_log10() +
  scale_color_viridis_d(option="B",labels=c(150,450)) +
  scale_fill_viridis_d(option="B",labels=c(150,450)) +
  # geom_smooth(aes(color=variable),fill="black",method='gam') +
  theme_pubclean() +
  labs(color="Test Gene Length", fill="Test Gene Length") +
  theme(legend.position=c(0.8,0.5),
        legend.background = element_rect(fill = NA)) +
  xlab("Train Gene Length") +
  ylab("Mean Square Prediction Error")

pCa <- ggplot(pred.150 %>% as.data.frame(),
       aes(x=log(2)/d,y=log(2)/exp(d.150))) +
  geom_point(size=2,alpha=0.25,pch=21,fill="gray50") +
  scale_x_log10(limits=c(0.01,3.5)) +
  scale_y_log10(limits=c(0.01,3.5)) +
  geom_abline(slope=1,intercept=0,lty=2,color="#357BA2FF",lwd=1.2) +
  annotate("text",label="X=Y",angle=52,x=3,y=2,size=3.5,color="#357BA2FF") +
  theme_pubclean() +
  xlab("Actual Avg. Max. Growth Rate") +
  ylab("Predicted Avg. Max. Growth Rate (1/Hr)") +
  ggtitle("Train 150bp/Test 150bp")

pCb <- ggplot(pred.150 %>% as.data.frame(),
       aes(x=log(2)/d,y=log(2)/exp(d.450))) +
  geom_point(size=2,alpha=0.25,pch=21,fill="gray50") +
  scale_x_log10(limits=c(0.01,3.5)) +
  scale_y_log10(limits=c(0.01,3.5)) +
  geom_abline(slope=1,intercept=0,lty=2,color="#357BA2FF",lwd=1.2) +
  annotate("text",label="X=Y",angle=52,x=3,y=2,size=3,color="#357BA2FF") +
  theme_pubclean() +
  xlab("Actual Avg. Max. Growth Rate") +
  ylab("Predicted Avg. Max. Growth Rate (1/Hr)") +
  ggtitle("Train 450bp/Test 150bp")


pCc <- ggplot(pred.450 %>% as.data.frame(),
       aes(x=log(2)/d,y=log(2)/exp(d.450))) +
  geom_point(size=2,alpha=0.25,pch=21,fill="gray50") +
  scale_x_log10(limits=c(0.01,3.5)) +
  scale_y_log10(limits=c(0.01,3.5)) +
  geom_abline(slope=1,intercept=0,lty=2,color="#357BA2FF",lwd=1.2) +
  annotate("text",label="X=Y",angle=52,x=3,y=2,size=3,color="#357BA2FF") +
  theme_pubclean() +
  xlab("Actual Avg. Max. Growth Rate") +
  ylab("Predicted Avg. Max. Growth Rate (1/Hr)") +
  ggtitle("Train 450bp/Test 450bp")

fast_df <- data.frame(Fast=mix_df$Fast,
                      pred.150.150=pred.150[,3],
                      pred.450.150=pred.150[,7])

pCf <- ggplot(fast_df,
       aes(x=Fast,y=exp(pred.150.150),group=Fast))  +
  geom_violin(fill="gray",color="gray") +
  geom_jitter(alpha=0.25,width=0.005) +
  geom_boxplot(width=0.03,fill=rgb(0,0,0,0),color="white",outliers=F) +
  theme_pubclean() +
  scale_y_log10() +
  xlab("Percent of Community with Doubling Times <1hr") +
  ylab("Predicted Avg. Min. Doubling Time (Hr)") +
  geom_hline(yintercept=1,lty=2,color="red")

cor.test(fast_df$Fast,fast_df$pred.150.150)

setwd("~/../../Documents/gRodon-assembly-free/figures")
png(file="mixtures_performance.png",width=9,height=8,units="in",res=600)
ggarrange(ggarrange(pCM,pCf,
                    ncol=2,
                    labels=c("(a)","(b)"),hjust=0),
          ggarrange(pCa,pCb,pCc,
                    ncol=3,
                    labels=c("(c)","(d)","(e)")),
          nrow=2)
dev.off()


# Synthetic Metagenomes ---------------------------------------------------


setwd("~/../../Documents/gRodon-assembly-free/data")
load("GrowthRates_Madin.rda")
load("Accession2Species_Madin.rda")
load("CodonStatistics_Madin.rda")
load("CodonStatistics_simmeta150.rda")
mix_names <- read.delim("simmeta_mixture_abundances.tsv",sep="\t",header=F)
names(mix_names) <- c("SIM","File","RelAb")

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

cu$Fast <- cu$d<1
mix_df <- sim_df %>%
  merge.easy(mix_names,.,key="SIM") %>%
  mutate(Accession=substr(File,1,13)) %>%
  merge.easy(.,cu,key="Accession") %>%
  mutate(mu=log(2)/d,
         mu.mmv2=log(2)/d.mmv2) %>%
  subset(select=c(MILC,ENCprime,B,SCUO,MCB,d.150,
                  d,d.mmv2,Fast,SIM,mu.mmv2,mu,RelAb)) %>%
  mutate(dW=log(d)*RelAb,
         FW=Fast*RelAb) %>%
  group_by(SIM) %>%
  summarise(MILC=mean(MILC,na.rm=T),
            ENCprime=mean(ENCprime,na.rm=T),
            B=mean(B,na.rm=T),
            SCUO=mean(SCUO,na.rm=T),
            MCB=mean(MCB,na.rm=T),
            d.150=mean(d.150,na.rm=T),
            d.mmv2=mean(d.mmv2,na.rm=T),
            d=exp(sum(dW)),
            Fast=sum(FW)) %>%
  mutate(Err.150=(log(d)-log(d.150))^2,
         Err.mmv2=(log(d)-log(d.mmv2))^2) %>%
  as.data.frame()

err_df <- mix_df %>%
  subset(select=c(Err.150,Err.mmv2)) %>%
  reshape2::melt()
pS1 <- ggplot(err_df,
       aes(x=variable,group=variable,y=value)) +
  geom_violin(fill="gray",color="gray") +
  geom_jitter(alpha=0.1,width=0.005) +
  geom_boxplot(width=0.25,fill=rgb(0,0,0,0),color="white",outliers=F) +
  theme_pubclean() +
  scale_y_log10() +
  scale_x_discrete(labels=c("Short-Read Mode","MMv2")) +
  xlab("") +
  ylab("Squared Prediction Error") +
  stat_pwc() +
  theme(axis.text.x=element_text(angle = 60, hjust = 1, vjust = 1))

pS2 <- ggplot(mix_df,
              aes(x=Fast,y=d.150))  +
  geom_point(alpha=0.1) +
  theme_pubclean() +
  scale_y_log10() +
  geom_smooth(color="white",fill="black") +
  xlab("Percent of Community with Doubling Times <1hr") +
  ylab("Predicted Avg. Min. Doubling Time (Hr)") +
  geom_hline(yintercept=1,lty=2,color="red")

cor.test(mix_df$Fast,mix_df$d.150)
wilcox.test(value~variable,data=err_df)

pS3 <- ggplot(mix_df,
              aes(x=log(2)/d,y=log(2)/d.150)) +
  geom_point(size=2,alpha=0.1,pch=21,fill="gray50") +
  scale_x_log10(limits=c(0.025,4.5)) +
  scale_y_log10(limits=c(0.025,4.5)) +
  geom_abline(slope=1,intercept=0,lty=2,color="#357BA2FF",lwd=1.2) +
  annotate("text",label="X=Y",angle=52,x=3,y=2,size=3.5,color="#357BA2FF") +
  theme_pubclean() +
  xlab("Actual Avg. Max. Growth Rate") +
  ylab("Predicted Avg. Max. Growth Rate (1/Hr)") +
  ggtitle("Short-Read Mode")


pS4 <- ggplot(mix_df,
              aes(x=log(2)/d,y=log(2)/d.mmv2)) +
  geom_point(size=2,alpha=0.1,pch=21,fill="gray50") +
  scale_x_log10(limits=c(0.025,4.5)) +
  scale_y_log10(limits=c(0.025,4.5)) +
  geom_abline(slope=1,intercept=0,lty=2,color="#357BA2FF",lwd=1.2) +
  annotate("text",label="X=Y",angle=52,x=3,y=2,size=3.5,color="#357BA2FF") +
  theme_pubclean() +
  xlab("Actual Avg. Max. Growth Rate") +
  ylab("Predicted Avg. Max. Growth Rate (1/Hr)") +
  ggtitle("Metagenome Mode v2")


setwd("~/../../Documents/gRodon-assembly-free/figures")
png(file="synthetic_metagenomes_performance.png",width=7,height=9,units="in",res=600)
ggarrange(ggarrange(pS1,
                    pS2,
                    ncol=2,
                    widths=c(5,10),
                    labels=c("(a)","(b)"),
                    hjust=0),
          ggarrange(pS3,
                    pS4,
                    ncol=2,
                    widths=c(1,1),
                    labels=c("(c)","(d)")),
          nrow=2)
dev.off()





