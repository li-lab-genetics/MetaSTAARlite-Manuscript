system("gsutil cp -r gs://fc-secure-797107a7-4402-4122-941c-9a486e0d633e/MetaSTAAR/AoU_Sumstats test/")
rm(list = ls())

coding_file_size_AoU <- 0
for(i in 1:381){
  coding_file_size_AoU <- coding_file_size_AoU + file.size(paste0("test/AoU_Sumstats/AoU_TC_Coding_sumstat_",i,".Rdata"))/(1024^3)
}

print(coding_file_size_AoU)

rm(list = ls())

for(j in 1:381){
  load(paste0("test/AoU_Sumstats/AoU_TC_Coding_sumstat_",j,".Rdata"))
  for(i in 1:length(coding_sumstat)){
    for(name in names(coding_sumstat[[i]])){
      coding_sumstat[[i]][[name]] <- subset(coding_sumstat[[i]][[name]],select = -c(CADD,LINSIGHT,FATHMM.XF,aPC.EpigeneticActive,aPC.EpigeneticRepressed,aPC.EpigeneticTranscription,aPC.Conservation,aPC.LocalDiversity,aPC.Mappability,aPC.TF,aPC.Protein))
    }
  }
  save(coding_sumstat,file=paste0("test/AoU_Sumstats/AoU_TC_Coding_sumstat_",j,".Rdata"),compress = "xz") 
}

rm(list = ls())

coding_file_size_AoU <- 0
for(i in 1:381){
  coding_file_size_AoU <- coding_file_size_AoU + file.size(paste0("test/AoU_Sumstats/AoU_TC_Coding_sumstat_",i,".Rdata"))/(1024^3)
}

print(coding_file_size_AoU)

rm(list = ls())

coding_file_size_AoU <- 0
for(i in 1:381){
  coding_file_size_AoU <- coding_file_size_AoU + file.size(paste0("test/AoU_Sumstats/AoU_TC_Coding_cov_",i,".Rdata"))/(1024^3)
}

print(coding_file_size_AoU)

rm(list = ls())

coding_time <- data.frame(job = 1:381)
for(i in 1:381){
  coding_time$time[i] <- unname(get(load(paste0("test/AoU_Sumstats/AoU_TC_Coding_time_",i,".Rdata"))))/3600
}

print(sum(coding_time$time))


































system("dx download -r UKB_PRS:MetaSTAAR/Results_470/")

rm(list = ls())

coding_file_size_470 <- 0
for(i in 1:381){
  coding_file_size_470 <- coding_file_size_470 + file.size(paste0("Results_470/UKB_TC_Coding_sumstat_",i,"_470.Rdata"))/(1024^3)
}

coding_file_size_470

rm(list = ls())

for(j in 1:381){
  load(paste0("Results_470/UKB_TC_Coding_sumstat_",j,"_470.Rdata"))
  for(i in 1:length(coding_sumstat)){
    for(name in names(coding_sumstat[[i]])){
      coding_sumstat[[i]][[name]] <- subset(coding_sumstat[[i]][[name]],select = -c(CADD,LINSIGHT,FATHMM.XF,aPC.EpigeneticActive,aPC.EpigeneticRepressed,aPC.EpigeneticTranscription,aPC.Conservation,aPC.LocalDiversity,aPC.Mappability,aPC.TF,aPC.Protein))
    }
  }
  save(coding_sumstat,file=paste0("Results_470/UKB_TC_Coding_sumstat_",j,"_470.Rdata"),compress = "xz") 
}

rm(list = ls())

coding_file_size_470 <- 0
for(i in 1:381){
  coding_file_size_470 <- coding_file_size_470 + file.size(paste0("Results_470/UKB_TC_Coding_sumstat_",i,"_470.Rdata"))/(1024^3)
}

coding_file_size_470

rm(list = ls())

coding_file_size_470 <- 0
for(i in 1:381){
  coding_file_size_470 <- coding_file_size_470 + file.size(paste0("Results_470/UKB_TC_Coding_cov_",i,"_470.Rdata"))/(1024^3)
}

coding_file_size_470

rm(list = ls())

coding_time <- data.frame(job = 1:381)
for(i in 1:381){
  coding_time$time[i] <- unname(get(load(paste0("Results_470/UKB_TC_Coding_time_",i,"_470.Rdata"))))/3600
}
coding_time$instance_type <- "V8"
coding_time$instance_type[c(8,9,28,55,58,86,108,112,113,132,138,145,163,204,216,255,263,268,282,332)] <- "V16"
coding_time$instance_type[c(112,113,332)] <- "V32"
coding_time$instance_type[c(380,381)] <- "V48"

sum(coding_time$time)
sum(c(coding_time$time[coding_time$instance_type == "V8"]*0.0528,coding_time$time[coding_time$instance_type == "V16"]*0.1056,coding_time$time[coding_time$instance_type == "V32"]*0.2112,coding_time$time[coding_time$instance_type == "V48"]*0.4224))





















































system("dx download -r UKB_PRS:MetaSTAAR/Results_31685/")
system("dx download -r UKB_PRS:MetaSTAAR/Results_63370/")
system("dx download -r UKB_PRS:MetaSTAAR/Results_95055/")

rm(list = ls())

coding_file_size_31685 <- 0
for(i in 1:381){
  coding_file_size_31685 <- coding_file_size_31685 + file.size(paste0("Results_31685/UKB_TC_Coding_sumstat_",i,"_31685.Rdata"))/(1024^3)
}

coding_file_size_63370 <- 0
for(i in 1:381){
  coding_file_size_63370 <- coding_file_size_63370 + file.size(paste0("Results_63370/UKB_TC_Coding_sumstat_",i,"_63370.Rdata"))/(1024^3)
}

coding_file_size_95055 <- 0
for(i in 1:381){
  coding_file_size_95055 <- coding_file_size_95055 + file.size(paste0("Results_95055/UKB_TC_Coding_sumstat_",i,"_95055.Rdata"))/(1024^3)
}

coding_file_size_31685
coding_file_size_63370
coding_file_size_95055

rm(list = ls())

for(j in 1:381){
  load(paste0("Results_31685/UKB_TC_Coding_sumstat_",j,"_31685.Rdata"))
  for(i in 1:length(coding_sumstat)){
    for(name in names(coding_sumstat[[i]])){
      coding_sumstat[[i]][[name]] <- subset(coding_sumstat[[i]][[name]],select = -c(CADD,LINSIGHT,FATHMM.XF,aPC.EpigeneticActive,aPC.EpigeneticRepressed,aPC.EpigeneticTranscription,aPC.Conservation,aPC.LocalDiversity,aPC.Mappability,aPC.TF,aPC.Protein))
    }
  }
  save(coding_sumstat,file=paste0("Results_31685/UKB_TC_Coding_sumstat_",j,"_31685.Rdata"),compress = "xz") 
}

rm(list = ls())

for(j in 1:381){
  load(paste0("Results_63370/UKB_TC_Coding_sumstat_",j,"_63370.Rdata"))
  for(i in 1:length(coding_sumstat)){
    for(name in names(coding_sumstat[[i]])){
      coding_sumstat[[i]][[name]] <- subset(coding_sumstat[[i]][[name]],select = -c(CADD,LINSIGHT,FATHMM.XF,aPC.EpigeneticActive,aPC.EpigeneticRepressed,aPC.EpigeneticTranscription,aPC.Conservation,aPC.LocalDiversity,aPC.Mappability,aPC.TF,aPC.Protein))
    }
  }
  save(coding_sumstat,file=paste0("Results_63370/UKB_TC_Coding_sumstat_",j,"_63370.Rdata"),compress = "xz") 
}

rm(list = ls())

for(j in 1:381){
  load(paste0("Results_95055/UKB_TC_Coding_sumstat_",j,"_95055.Rdata"))
  for(i in 1:length(coding_sumstat)){
    for(name in names(coding_sumstat[[i]])){
      coding_sumstat[[i]][[name]] <- subset(coding_sumstat[[i]][[name]],select = -c(CADD,LINSIGHT,FATHMM.XF,aPC.EpigeneticActive,aPC.EpigeneticRepressed,aPC.EpigeneticTranscription,aPC.Conservation,aPC.LocalDiversity,aPC.Mappability,aPC.TF,aPC.Protein))
    }
  }
  save(coding_sumstat,file=paste0("Results_95055/UKB_TC_Coding_sumstat_",j,"_95055.Rdata"),compress = "xz") 
}




rm(list = ls())

coding_file_size_31685 <- 0
for(i in 1:381){
  coding_file_size_31685 <- coding_file_size_31685 + file.size(paste0("Results_31685/UKB_TC_Coding_sumstat_",i,"_31685.Rdata"))/(1024^3)
}

coding_file_size_63370 <- 0
for(i in 1:381){
  coding_file_size_63370 <- coding_file_size_63370 + file.size(paste0("Results_63370/UKB_TC_Coding_sumstat_",i,"_63370.Rdata"))/(1024^3)
}

coding_file_size_95055 <- 0
for(i in 1:381){
  coding_file_size_95055 <- coding_file_size_95055 + file.size(paste0("Results_95055/UKB_TC_Coding_sumstat_",i,"_95055.Rdata"))/(1024^3)
}

coding_file_size_31685
coding_file_size_63370
coding_file_size_95055

rm(list = ls())

coding_file_size_31685 <- 0
for(i in 1:381){
  coding_file_size_31685 <- coding_file_size_31685 + file.size(paste0("Results_31685/UKB_TC_Coding_cov_",i,"_31685.Rdata"))/(1024^3)
}

coding_file_size_63370 <- 0
for(i in 1:381){
  coding_file_size_63370 <- coding_file_size_63370 + file.size(paste0("Results_63370/UKB_TC_Coding_cov_",i,"_63370.Rdata"))/(1024^3)
}

coding_file_size_95055 <- 0
for(i in 1:381){
  coding_file_size_95055 <- coding_file_size_95055 + file.size(paste0("Results_95055/UKB_TC_Coding_cov_",i,"_95055.Rdata"))/(1024^3)
}

coding_file_size_31685
coding_file_size_63370
coding_file_size_95055



rm(list = ls())

coding_time <- data.frame(job = 1:381)
for(i in 1:381){
  coding_time$time[i] <- unname(get(load(paste0("Results_31685/UKB_TC_Coding_time_",i,"_31685.Rdata"))))/3600
}
coding_time$instance_type <- "V2"
coding_time$instance_type[c(380,381)] <- "V8"

sum(coding_time$time)
sum(c(coding_time$time[coding_time$instance_type == "V2"]*0.0132,coding_time$time[coding_time$instance_type == "V8"]*0.0528))

rm(list = ls())

coding_time <- data.frame(job = 1:381)
for(i in 1:381){
  coding_time$time[i] <- unname(get(load(paste0("Results_63370/UKB_TC_Coding_time_",i,"_63370.Rdata"))))/3600
}
coding_time$instance_type <- "V2"
coding_time$instance_type[c(112,332)] <- "V4"
coding_time$instance_type[c(380,381)] <- "V8"

sum(coding_time$time)
sum(c(coding_time$time[coding_time$instance_type == "V2"]*0.0132,coding_time$time[coding_time$instance_type == "V4"]*0.0264,coding_time$time[coding_time$instance_type == "V8"]*0.0528))


rm(list = ls())

coding_time <- data.frame(job = 1:381)
for(i in 1:381){
  coding_time$time[i] <- unname(get(load(paste0("Results_95055/UKB_TC_Coding_time_",i,"_95055.Rdata"))))/3600
}
coding_time$instance_type <- "V2"
coding_time$instance_type[c(8,9,112,113,132,138,332)] <- "V4"
coding_time$instance_type[c(380,381)] <- "V8"

sum(coding_time$time)
sum(c(coding_time$time[coding_time$instance_type == "V2"]*0.0132,coding_time$time[coding_time$instance_type == "V4"]*0.0264,coding_time$time[coding_time$instance_type == "V8"]*0.0528))























rm(list = ls())

ncRNA_file_size_31685 <- 0
for(i in 1:223){
  ncRNA_file_size_31685 <- ncRNA_file_size_31685 + file.size(paste0("Results_31685/UKB_TC_ncRNA_sumstat_",i,"_31685.Rdata"))/(1024^3)
}

ncRNA_file_size_63370 <- 0
for(i in 1:223){
  ncRNA_file_size_63370 <- ncRNA_file_size_63370 + file.size(paste0("Results_63370/UKB_TC_ncRNA_sumstat_",i,"_63370.Rdata"))/(1024^3)
}

ncRNA_file_size_95055 <- 0
for(i in 1:223){
  ncRNA_file_size_95055 <- ncRNA_file_size_95055 + file.size(paste0("Results_95055/UKB_TC_ncRNA_sumstat_",i,"_95055.Rdata"))/(1024^3)
}

ncRNA_file_size_31685
ncRNA_file_size_63370
ncRNA_file_size_95055

rm(list = ls())

for(j in 1:223){
  load(paste0("Results_31685/UKB_TC_ncRNA_sumstat_",j,"_31685.Rdata"))
  for(i in 1:length(ncRNA_sumstat)){
    ncRNA_sumstat[[i]] <- subset(ncRNA_sumstat[[i]],select = -c(CADD,LINSIGHT,FATHMM.XF,aPC.EpigeneticActive,aPC.EpigeneticRepressed,aPC.EpigeneticTranscription,aPC.Conservation,aPC.LocalDiversity,aPC.Mappability,aPC.TF,aPC.Protein))
  }
  save(ncRNA_sumstat,file=paste0("Results_31685/UKB_TC_ncRNA_sumstat_",j,"_31685.Rdata"),compress = "xz") 
}

rm(list = ls())

for(j in 1:223){
  load(paste0("Results_63370/UKB_TC_ncRNA_sumstat_",j,"_63370.Rdata"))
  for(i in 1:length(ncRNA_sumstat)){
    ncRNA_sumstat[[i]] <- subset(ncRNA_sumstat[[i]],select = -c(CADD,LINSIGHT,FATHMM.XF,aPC.EpigeneticActive,aPC.EpigeneticRepressed,aPC.EpigeneticTranscription,aPC.Conservation,aPC.LocalDiversity,aPC.Mappability,aPC.TF,aPC.Protein))
  }
  save(ncRNA_sumstat,file=paste0("Results_63370/UKB_TC_ncRNA_sumstat_",j,"_63370.Rdata"),compress = "xz") 
}

rm(list = ls())

for(j in 1:223){
  load(paste0("Results_95055/UKB_TC_ncRNA_sumstat_",j,"_95055.Rdata"))
  for(i in 1:length(ncRNA_sumstat)){
    ncRNA_sumstat[[i]] <- subset(ncRNA_sumstat[[i]],select = -c(CADD,LINSIGHT,FATHMM.XF,aPC.EpigeneticActive,aPC.EpigeneticRepressed,aPC.EpigeneticTranscription,aPC.Conservation,aPC.LocalDiversity,aPC.Mappability,aPC.TF,aPC.Protein))
  }
  save(ncRNA_sumstat,file=paste0("Results_95055/UKB_TC_ncRNA_sumstat_",j,"_95055.Rdata"),compress = "xz") 
}




rm(list = ls())

ncRNA_file_size_31685 <- 0
for(i in 1:223){
  ncRNA_file_size_31685 <- ncRNA_file_size_31685 + file.size(paste0("Results_31685/UKB_TC_ncRNA_sumstat_",i,"_31685.Rdata"))/(1024^3)
}

ncRNA_file_size_63370 <- 0
for(i in 1:223){
  ncRNA_file_size_63370 <- ncRNA_file_size_63370 + file.size(paste0("Results_63370/UKB_TC_ncRNA_sumstat_",i,"_63370.Rdata"))/(1024^3)
}

ncRNA_file_size_95055 <- 0
for(i in 1:223){
  ncRNA_file_size_95055 <- ncRNA_file_size_95055 + file.size(paste0("Results_95055/UKB_TC_ncRNA_sumstat_",i,"_95055.Rdata"))/(1024^3)
}

ncRNA_file_size_31685
ncRNA_file_size_63370
ncRNA_file_size_95055

rm(list = ls())

ncRNA_file_size_31685 <- 0
for(i in 1:223){
  ncRNA_file_size_31685 <- ncRNA_file_size_31685 + file.size(paste0("Results_31685/UKB_TC_ncRNA_cov_",i,"_31685.Rdata"))/(1024^3)
}

ncRNA_file_size_63370 <- 0
for(i in 1:223){
  ncRNA_file_size_63370 <- ncRNA_file_size_63370 + file.size(paste0("Results_63370/UKB_TC_ncRNA_cov_",i,"_63370.Rdata"))/(1024^3)
}

ncRNA_file_size_95055 <- 0
for(i in 1:223){
  ncRNA_file_size_95055 <- ncRNA_file_size_95055 + file.size(paste0("Results_95055/UKB_TC_ncRNA_cov_",i,"_95055.Rdata"))/(1024^3)
}

ncRNA_file_size_31685
ncRNA_file_size_63370
ncRNA_file_size_95055



rm(list = ls())

ncRNA_time <- data.frame(job = 1:223)
for(i in 1:223){
  ncRNA_time$time[i] <- unname(get(load(paste0("Results_31685/UKB_TC_ncRNA_time_",i,"_31685.Rdata"))))/3600
}
ncRNA_time$instance_type <- "V2"
ncRNA_time$instance_type[c(223)] <- "V16"

sum(ncRNA_time$time)
sum(c(ncRNA_time$time[ncRNA_time$instance_type == "V2"]*0.0132,ncRNA_time$time[ncRNA_time$instance_type == "V16"]*0.1056))

rm(list = ls())

ncRNA_time <- data.frame(job = 1:223)
for(i in 1:223){
  ncRNA_time$time[i] <- unname(get(load(paste0("Results_63370/UKB_TC_ncRNA_time_",i,"_63370.Rdata"))))/3600
}
ncRNA_time$instance_type <- "V2"
ncRNA_time$instance_type[c(24,28,44,61,69,70,97,98,126,153,156,187,201,212,214)] <- "V4"
ncRNA_time$instance_type[c(223)] <- "V16"

sum(ncRNA_time$time)
sum(c(ncRNA_time$time[ncRNA_time$instance_type == "V2"]*0.0132,ncRNA_time$time[ncRNA_time$instance_type == "V4"]*0.0264,ncRNA_time$time[ncRNA_time$instance_type == "V16"]*0.1056))


rm(list = ls())

ncRNA_time <- data.frame(job = 1:223)
for(i in 1:223){
  ncRNA_time$time[i] <- unname(get(load(paste0("Results_95055/UKB_TC_ncRNA_time_",i,"_95055.Rdata"))))/3600
}
ncRNA_time$instance_type <- "V4"
ncRNA_time$instance_type[c(69,97,201,214)] <- "V8"
ncRNA_time$instance_type[c(223)] <- "V16"

sum(ncRNA_time$time)
sum(c(ncRNA_time$time[ncRNA_time$instance_type == "V2"]*0.0132,ncRNA_time$time[ncRNA_time$instance_type == "V4"]*0.0264,ncRNA_time$time[ncRNA_time$instance_type == "V8"]*0.0528,ncRNA_time$time[ncRNA_time$instance_type == "V16"]*0.1056))













rm(list = ls())

Noncoding_file_size_31685 <- 0
for(i in 1:387){
  Noncoding_file_size_31685 <- Noncoding_file_size_31685 + file.size(paste0("Results_31685/UKB_TC_Noncoding_sumstat_",i,"_31685.Rdata"))/(1024^3)
}

Noncoding_file_size_63370 <- 0
for(i in 1:387){
  Noncoding_file_size_63370 <- Noncoding_file_size_63370 + file.size(paste0("Results_63370/UKB_TC_Noncoding_sumstat_",i,"_63370.Rdata"))/(1024^3)
}

Noncoding_file_size_95055 <- 0
for(i in 1:387){
  Noncoding_file_size_95055 <- Noncoding_file_size_95055 + file.size(paste0("Results_95055/UKB_TC_Noncoding_sumstat_",i,"_95055.Rdata"))/(1024^3)
}

Noncoding_file_size_31685
Noncoding_file_size_63370
Noncoding_file_size_95055

rm(list = ls())

for(j in 1:387){
  load(paste0("Results_31685/UKB_TC_Noncoding_sumstat_",j,"_31685.Rdata"))
  for(i in 1:length(noncoding_sumstat)){
    for(name in names(noncoding_sumstat[[i]])){
      noncoding_sumstat[[i]][[name]] <- subset(noncoding_sumstat[[i]][[name]],select = -c(CADD,LINSIGHT,FATHMM.XF,aPC.EpigeneticActive,aPC.EpigeneticRepressed,aPC.EpigeneticTranscription,aPC.Conservation,aPC.LocalDiversity,aPC.Mappability,aPC.TF,aPC.Protein))
    }
  }
  save(noncoding_sumstat,file=paste0("Results_31685/UKB_TC_Noncoding_sumstat_",j,"_31685.Rdata"),compress = "xz") 
}

rm(list = ls())

for(j in 1:387){
  load(paste0("Results_63370/UKB_TC_Noncoding_sumstat_",j,"_63370.Rdata"))
  for(i in 1:length(noncoding_sumstat)){
    for(name in names(noncoding_sumstat[[i]])){
      noncoding_sumstat[[i]][[name]] <- subset(noncoding_sumstat[[i]][[name]],select = -c(CADD,LINSIGHT,FATHMM.XF,aPC.EpigeneticActive,aPC.EpigeneticRepressed,aPC.EpigeneticTranscription,aPC.Conservation,aPC.LocalDiversity,aPC.Mappability,aPC.TF,aPC.Protein))
    }
  }
  save(noncoding_sumstat,file=paste0("Results_63370/UKB_TC_Noncoding_sumstat_",j,"_63370.Rdata"),compress = "xz") 
}

rm(list = ls())

for(j in 1:387){
  load(paste0("Results_95055/UKB_TC_Noncoding_sumstat_",j,"_95055.Rdata"))
  for(i in 1:length(noncoding_sumstat)){
    for(name in names(noncoding_sumstat[[i]])){
      noncoding_sumstat[[i]][[name]] <- subset(noncoding_sumstat[[i]][[name]],select = -c(CADD,LINSIGHT,FATHMM.XF,aPC.EpigeneticActive,aPC.EpigeneticRepressed,aPC.EpigeneticTranscription,aPC.Conservation,aPC.LocalDiversity,aPC.Mappability,aPC.TF,aPC.Protein))
    }
  }
  save(noncoding_sumstat,file=paste0("Results_95055/UKB_TC_Noncoding_sumstat_",j,"_95055.Rdata"),compress = "xz") 
}




rm(list = ls())

Noncoding_file_size_31685 <- 0
for(i in 1:387){
  Noncoding_file_size_31685 <- Noncoding_file_size_31685 + file.size(paste0("Results_31685/UKB_TC_Noncoding_sumstat_",i,"_31685.Rdata"))/(1024^3)
}

Noncoding_file_size_63370 <- 0
for(i in 1:387){
  Noncoding_file_size_63370 <- Noncoding_file_size_63370 + file.size(paste0("Results_63370/UKB_TC_Noncoding_sumstat_",i,"_63370.Rdata"))/(1024^3)
}

Noncoding_file_size_95055 <- 0
for(i in 1:387){
  Noncoding_file_size_95055 <- Noncoding_file_size_95055 + file.size(paste0("Results_95055/UKB_TC_Noncoding_sumstat_",i,"_95055.Rdata"))/(1024^3)
}

Noncoding_file_size_31685
Noncoding_file_size_63370
Noncoding_file_size_95055

rm(list = ls())

Noncoding_file_size_31685 <- 0
for(i in 1:387){
  Noncoding_file_size_31685 <- Noncoding_file_size_31685 + file.size(paste0("Results_31685/UKB_TC_Noncoding_cov_",i,"_31685.Rdata"))/(1024^3)
}

Noncoding_file_size_63370 <- 0
for(i in 1:387){
  Noncoding_file_size_63370 <- Noncoding_file_size_63370 + file.size(paste0("Results_63370/UKB_TC_Noncoding_cov_",i,"_63370.Rdata"))/(1024^3)
}

Noncoding_file_size_95055 <- 0
for(i in 1:387){
  Noncoding_file_size_95055 <- Noncoding_file_size_95055 + file.size(paste0("Results_95055/UKB_TC_Noncoding_cov_",i,"_95055.Rdata"))/(1024^3)
}

Noncoding_file_size_31685
Noncoding_file_size_63370
Noncoding_file_size_95055



rm(list = ls())

Noncoding_time <- data.frame(job = 1:387)
for(i in 1:387){
  Noncoding_time$time[i] <- unname(get(load(paste0("Results_31685/UKB_TC_Noncoding_time_",i,"_31685.Rdata"))))/3600
}
Noncoding_time$instance_type <- "V2"
Noncoding_time$instance_type[c(2,11,12,13,29,33,41,44,46,47,53,54,58,59,64,67,75,85,168,176,191,195,208,319,324)] <- "V4"
Noncoding_time$instance_type[380:387] <- "V16"

sum(Noncoding_time$time)
sum(c(Noncoding_time$time[Noncoding_time$instance_type == "V2"]*0.0132,Noncoding_time$time[Noncoding_time$instance_type == "V4"]*0.0264,Noncoding_time$time[Noncoding_time$instance_type == "V16"]*0.1056))

rm(list = ls())

Noncoding_time <- data.frame(job = 1:387)
for(i in 1:387){
  Noncoding_time$time[i] <- unname(get(load(paste0("Results_63370/UKB_TC_Noncoding_time_",i,"_63370.Rdata"))))/3600
}
Noncoding_time$instance_type <- "V4"
Noncoding_time$instance_type[c(1,2,8,20,42,56,58,62,66,100,124,135,176,179,191,195,319,324)] <- "V8"
Noncoding_time$instance_type[380:387] <- "V16"

sum(Noncoding_time$time)
sum(c(Noncoding_time$time[Noncoding_time$instance_type == "V2"]*0.0132,Noncoding_time$time[Noncoding_time$instance_type == "V4"]*0.0264,Noncoding_time$time[Noncoding_time$instance_type == "V8"]*0.0528,Noncoding_time$time[Noncoding_time$instance_type == "V16"]*0.1056))


rm(list = ls())

Noncoding_time <- data.frame(job = 1:387)
for(i in 1:387){
  Noncoding_time$time[i] <- unname(get(load(paste0("Results_95055/UKB_TC_Noncoding_time_",i,"_95055.Rdata"))))/3600
}
Noncoding_time$instance_type <- "V4"
Noncoding_time$instance_type[c(2,3,11,12,20,26,29,31,34,40,41,42,44,45,46,47,48,53,54,56,57,58,59,61,62,63,64,67,69,71,74,75,76,78,81,85,89,90,95,103,104,105,106,107,109,110,114,115,116,117,125,127,130,131,133,138,139,145,146,148,151,160,162,165,166,168,170,173,174,175,176,177,178,179,182,183,184,187,191,193,195,200,202,208,211,224,225,229,230,242,243,245,248,249,251,254,255,256,261,262,263,266,269,270,273,274,275,279,286,293,294,295,296,305,319,321,324,325,357,360,362,363,371,377,378)] <- "V8"
Noncoding_time$instance_type[380:387] <- "V16"

sum(Noncoding_time$time)
sum(c(Noncoding_time$time[Noncoding_time$instance_type == "V2"]*0.0132,Noncoding_time$time[Noncoding_time$instance_type == "V4"]*0.0264,Noncoding_time$time[Noncoding_time$instance_type == "V8"]*0.0528),Noncoding_time$time[Noncoding_time$instance_type == "V16"]*0.1056)
































































# system("dx download -r UKB_PRS:MetaSTAAR/Results/")

rm(list = ls())

noncoding_file_size <- 0
for(i in 1:387){
  noncoding_file_size <- noncoding_file_size + file.size(paste0("Results/UKB_TC_Noncoding_sumstat_",i,".Rdata"))/(1024^3)
}

coding_file_size <- 0
for(i in 1:381){
  coding_file_size <- coding_file_size + file.size(paste0("Results/UKB_TC_Coding_sumstat_",i,".Rdata"))/(1024^3)
}

ncRNA_file_size <- 0
for(i in 1:223){
  ncRNA_file_size <- ncRNA_file_size + file.size(paste0("Results/UKB_TC_ncRNA_sumstat_",i,".Rdata"))/(1024^3)
}

noncoding_file_size
coding_file_size
ncRNA_file_size

rm(list = ls())

for(j in 1:381){
  load(paste0("Results/UKB_TC_Coding_sumstat_",j,".Rdata"))
  for(i in 1:length(coding_sumstat)){
    for(name in names(coding_sumstat[[i]])){
      coding_sumstat[[i]][[name]] <- subset(coding_sumstat[[i]][[name]],select = -c(CADD,LINSIGHT,FATHMM.XF,aPC.EpigeneticActive,aPC.EpigeneticRepressed,aPC.EpigeneticTranscription,aPC.Conservation,aPC.LocalDiversity,aPC.Mappability,aPC.TF,aPC.Protein))
    }
  }
  save(coding_sumstat,file=paste0("Results/UKB_TC_Coding_sumstat_",j,".Rdata"),compress = "xz") 
}

rm(list = ls())

for(j in 1:387){
  load(paste0("Results/UKB_TC_Noncoding_sumstat_",j,".Rdata"))
  for(i in 1:length(noncoding_sumstat)){
    for(name in names(noncoding_sumstat[[i]])){
      noncoding_sumstat[[i]][[name]] <- subset(noncoding_sumstat[[i]][[name]],select = -c(CADD,LINSIGHT,FATHMM.XF,aPC.EpigeneticActive,aPC.EpigeneticRepressed,aPC.EpigeneticTranscription,aPC.Conservation,aPC.LocalDiversity,aPC.Mappability,aPC.TF,aPC.Protein))
    }
  }
  save(noncoding_sumstat,file=paste0("Results/UKB_TC_Noncoding_sumstat_",j,".Rdata"),compress = "xz") 
}

rm(list = ls())

for(j in 1:223){
  load(paste0("Results/UKB_TC_ncRNA_sumstat_",j,".Rdata"))
  for(i in 1:length(ncRNA_sumstat)){
    ncRNA_sumstat[[i]] <- subset(ncRNA_sumstat[[i]],select = -c(CADD,LINSIGHT,FATHMM.XF,aPC.EpigeneticActive,aPC.EpigeneticRepressed,aPC.EpigeneticTranscription,aPC.Conservation,aPC.LocalDiversity,aPC.Mappability,aPC.TF,aPC.Protein))
  }
  save(ncRNA_sumstat,file=paste0("Results/UKB_TC_ncRNA_sumstat_",j,".Rdata"),compress = "xz") 
}


rm(list = ls())

noncoding_file_size <- 0
for(i in 1:387){
  noncoding_file_size <- noncoding_file_size + file.size(paste0("Results/UKB_TC_Noncoding_sumstat_",i,".Rdata"))/(1024^3)
}

coding_file_size <- 0
for(i in 1:381){
  coding_file_size <- coding_file_size + file.size(paste0("Results/UKB_TC_Coding_sumstat_",i,".Rdata"))/(1024^3)
}

ncRNA_file_size <- 0
for(i in 1:223){
  ncRNA_file_size <- ncRNA_file_size + file.size(paste0("Results/UKB_TC_ncRNA_sumstat_",i,".Rdata"))/(1024^3)
}

noncoding_file_size
coding_file_size
ncRNA_file_size



































rm(list = ls())

for(j in 1:381){
  load(paste0("Results/UKB_TC_Coding_sumstat_",j,".Rdata"))
  for(i in 1:length(coding_sumstat)){
    for(name in names(coding_sumstat[[i]])){
      coding_sumstat[[i]][[name]] <- coding_sumstat[[i]][[name]][coding_sumstat[[i]][[name]]$qc_label == "PASS" & coding_sumstat[[i]][[name]]$MAF < 0.05 & coding_sumstat[[i]][[name]]$MAF != 0, ]
      coding_sumstat[[i]][[name]] <- subset(coding_sumstat[[i]][[name]],select = -c(qc_label))
    }
  }
  save(coding_sumstat,file=paste0("Results/UKB_TC_Coding_sumstat_",j,".Rdata"),compress = "xz") 
}

rm(list = ls())

for(j in 1:387){
  load(paste0("Results/UKB_TC_Noncoding_sumstat_",j,".Rdata"))
  for(i in 1:length(noncoding_sumstat)){
    for(name in names(noncoding_sumstat[[i]])){
      noncoding_sumstat[[i]][[name]] <- noncoding_sumstat[[i]][[name]][noncoding_sumstat[[i]][[name]]$qc_label == "PASS" & noncoding_sumstat[[i]][[name]]$MAF < 0.05 & noncoding_sumstat[[i]][[name]]$MAF != 0, ]
      noncoding_sumstat[[i]][[name]] <- subset(noncoding_sumstat[[i]][[name]],select = -c(qc_label))
    }
  }
  save(noncoding_sumstat,file=paste0("Results/UKB_TC_Noncoding_sumstat_",j,".Rdata"),compress = "xz") 
}

rm(list = ls())

for(j in 1:223){
  load(paste0("Results/UKB_TC_ncRNA_sumstat_",j,".Rdata"))
  for(i in 1:length(ncRNA_sumstat)){
    ncRNA_sumstat[[i]] <- ncRNA_sumstat[[i]][ncRNA_sumstat[[i]]$qc_label == "PASS" & ncRNA_sumstat[[i]]$MAF < 0.05 & ncRNA_sumstat[[i]]$MAF != 0, ]
    ncRNA_sumstat[[i]] <- subset(ncRNA_sumstat[[i]],select = -c(qc_label))
  }
  save(ncRNA_sumstat,file=paste0("Results/UKB_TC_ncRNA_sumstat_",j,".Rdata"),compress = "xz") 
}


rm(list = ls())

noncoding_file_size <- 0
for(i in 1:387){
  noncoding_file_size <- noncoding_file_size + file.size(paste0("Results/UKB_TC_Noncoding_sumstat_",i,".Rdata"))/(1024^3)
}

coding_file_size <- 0
for(i in 1:381){
  coding_file_size <- coding_file_size + file.size(paste0("Results/UKB_TC_Coding_sumstat_",i,".Rdata"))/(1024^3)
}

ncRNA_file_size <- 0
for(i in 1:223){
  ncRNA_file_size <- ncRNA_file_size + file.size(paste0("Results/UKB_TC_ncRNA_sumstat_",i,".Rdata"))/(1024^3)
}

noncoding_file_size
coding_file_size
ncRNA_file_size

noncoding_file_size <- 0
for(i in 1:387){
  noncoding_file_size <- noncoding_file_size + file.size(paste0("Results/UKB_TC_Noncoding_cov_",i,".Rdata"))/(1024^3)
}

coding_file_size <- 0
for(i in 1:381){
  coding_file_size <- coding_file_size + file.size(paste0("Results/UKB_TC_Coding_cov_",i,".Rdata"))/(1024^3)
}

ncRNA_file_size <- 0
for(i in 1:223){
  ncRNA_file_size <- ncRNA_file_size + file.size(paste0("Results/UKB_TC_ncRNA_cov_",i,".Rdata"))/(1024^3)
}

noncoding_file_size
coding_file_size
ncRNA_file_size


rm(list = ls())

for(j in 1:381){
  load(paste0("Results/UKB_TC_Coding_sumstat_",j,".Rdata"))
  for(i in 1:length(coding_sumstat)){
    for(name in names(coding_sumstat[[i]])){
      coding_sumstat[[i]][[name]] <- subset(coding_sumstat[[i]][[name]],select = -c(CADD,LINSIGHT,FATHMM.XF,aPC.EpigeneticActive,aPC.EpigeneticRepressed,aPC.EpigeneticTranscription,aPC.Conservation,aPC.LocalDiversity,aPC.Mappability,aPC.TF,aPC.Protein))
    }
  }
  save(coding_sumstat,file=paste0("Results/UKB_TC_Coding_sumstat_",j,".Rdata"),compress = "xz") 
}

rm(list = ls())

for(j in 1:387){
  load(paste0("Results/UKB_TC_Noncoding_sumstat_",j,".Rdata"))
  for(i in 1:length(noncoding_sumstat)){
    for(name in names(noncoding_sumstat[[i]])){
      noncoding_sumstat[[i]][[name]] <- subset(noncoding_sumstat[[i]][[name]],select = -c(CADD,LINSIGHT,FATHMM.XF,aPC.EpigeneticActive,aPC.EpigeneticRepressed,aPC.EpigeneticTranscription,aPC.Conservation,aPC.LocalDiversity,aPC.Mappability,aPC.TF,aPC.Protein))
    }
  }
  save(noncoding_sumstat,file=paste0("Results/UKB_TC_Noncoding_sumstat_",j,".Rdata"),compress = "xz") 
}

rm(list = ls())

for(j in 1:223){
  load(paste0("Results/UKB_TC_ncRNA_sumstat_",j,".Rdata"))
  for(i in 1:length(ncRNA_sumstat)){
    ncRNA_sumstat[[i]] <- subset(ncRNA_sumstat[[i]],select = -c(CADD,LINSIGHT,FATHMM.XF,aPC.EpigeneticActive,aPC.EpigeneticRepressed,aPC.EpigeneticTranscription,aPC.Conservation,aPC.LocalDiversity,aPC.Mappability,aPC.TF,aPC.Protein))
  }
  save(ncRNA_sumstat,file=paste0("Results/UKB_TC_ncRNA_sumstat_",j,".Rdata"),compress = "xz") 
}

rm(list = ls())

noncoding_file_size <- 0
for(i in 1:387){
  noncoding_file_size <- noncoding_file_size + file.size(paste0("Results/UKB_TC_Noncoding_sumstat_",i,".Rdata"))/(1024^3)
}

coding_file_size <- 0
for(i in 1:381){
  coding_file_size <- coding_file_size + file.size(paste0("Results/UKB_TC_Coding_sumstat_",i,".Rdata"))/(1024^3)
}

ncRNA_file_size <- 0
for(i in 1:223){
  ncRNA_file_size <- ncRNA_file_size + file.size(paste0("Results/UKB_TC_ncRNA_sumstat_",i,".Rdata"))/(1024^3)
}

noncoding_file_size
coding_file_size
ncRNA_file_size


rm(list = ls())

noncoding_time <- data.frame(job = 1:387)
for(i in 1:387){
  noncoding_time$time[i] <- unname(get(load(paste0("Results/UKB_TC_Noncoding_time_",i,".Rdata"))))/3600
}
noncoding_time$instance_type <- "V8"
noncoding_time$instance_type[c(2,11,12,22,29,41,45,46,53,54,57,58,59,64,67,75,78,85,89,95,109,110,117,127,130,133,138,139,145,148,151,162,165,166,168,173,175,176,178,182,183,191,193,195,200,202,208,211,229,230,242,243,248,249,251,254,255,256,262,263,266,269,270,273,275,279,293,294,295,296,319,324,325,357,360,362,363,371)] <- "V16"
noncoding_time$instance_type[c(380,382,384,386,387)] <- "V16"
noncoding_time$instance_type[c(381,383,385)] <- "V32"

sum(noncoding_time$time)
sum(c(noncoding_time$time[noncoding_time$instance_type == "V8"]*0.0528,noncoding_time$time[noncoding_time$instance_type == "V16"]*0.1056,noncoding_time$time[noncoding_time$instance_type == "V32"]*0.2112))

rm(list = ls())

coding_time <- data.frame(job = 1:381)
for(i in 1:379){
  coding_time$time[i] <- unname(get(load(paste0("Results/UKB_TC_Coding_time_",i,".Rdata"))))/3600
}
for(i in 380:381){
  coding_time$time[i] <- unname(get(load(paste0("Results/UKB_TC_Coding_time_",i,".Rdata"))))[3]/3600
}
coding_time$instance_type <- "V4"
coding_time$instance_type[c(8,112,113,332)] <- "V8"
coding_time$instance_type[c(380,381)] <- "V16"

sum(coding_time$time)
sum(c(coding_time$time[coding_time$instance_type == "V4"]*0.0264,coding_time$time[coding_time$instance_type == "V8"]*0.0528,coding_time$time[coding_time$instance_type == "V16"]*0.1056))


rm(list = ls())

ncRNA_time <- data.frame(job = 1:223)
for(i in 1:223){
  ncRNA_time$time[i] <- unname(get(load(paste0("Results/UKB_TC_Coding_time_",i,".Rdata"))))/3600
}
ncRNA_time$instance_type <- "V4"
ncRNA_time$instance_type[c(1,19,24,26,28,39,44,46,52,54,59,61,65,70,73,98,101,113,114,121,126,142,153,156,159,161,187,212)] <- "V8"
ncRNA_time$instance_type[c(69,97,201,214)] <- "V16"
ncRNA_time$instance_type[c(223)] <- "V32"

sum(ncRNA_time$time)
sum(c(ncRNA_time$time[ncRNA_time$instance_type == "V4"]*0.0264,ncRNA_time$time[ncRNA_time$instance_type == "V8"]*0.0528,ncRNA_time$time[ncRNA_time$instance_type == "V16"]*0.1056,ncRNA_time$time[ncRNA_time$instance_type == "V32"]*0.2112))
