rm(list=ls())
gc()

time <- system.time({
  ## load required packages
  library(gdsfmt)
  library(SeqArray)
  library(SeqVarTools)
  library(STAAR)
  library(STAARpipeline)
  library(readr)
  library(dplyr)
  library(stringr)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  
  library(devtools)
  devtools::install_github("xihaoli/MetaSTAAR",ref="main")
  library(MetaSTAAR)
  library(expm)
  library(Matrix)
  
  # dx run app-swiss-army-knife -iin=UKB_PRS:JW/Software/r_with_plink.tar.gz -iin=UKB_PRS:MetaSTAAR/Scripts/Meta_ncRNA_63370_LongMasks.R -iin=UKB_PRS:MetaSTAAR/Scripts/Meta_ncRNA_63370_LongMasks.sh -icmd="bash Meta_ncRNA_63370_LongMasks.sh" -y --destination UKB_PRS:MetaSTAAR/Results_63370/ --priority low --extra-args '{"executionPolicy":{"maxSpotTries":5,"spotOnly":true}}' --instance-type mem3_ssd1_v2_x16

  
  MetaSTAARlite_worker_sumstat <- function(genotype,obj_nullmodel,variant_info,qc_label=NULL,
                                           annotation_phred=NULL,for_individual_analysis=FALSE){
    
    if(!inherits(genotype, "matrix") && !inherits(genotype, "Matrix")){
      stop("genotype is not a matrix!")
    }
    
    if(inherits(genotype, "sparseMatrix")){
      genotype <- as.matrix(genotype)
    }
    
    if(dim(genotype)[2] != dim(variant_info)[1]){
      stop(paste0("Dimensions don't match for genotype and variant_info!"))
    }
    
    if(!is.null(qc_label) && dim(variant_info)[1] != length(qc_label)){
      stop(paste0("Dimensions don't match for variant_info and qc_label!"))
    }
    
    if(!is.null(annotation_phred) && dim(annotation_phred)[1] != dim(variant_info)[1]){
      stop(paste0("Dimensions don't match for annotation_phred and variant_info!"))
    }
    
    N <- dim(genotype)[1]
    alt_AC <- as.integer(colSums(2 - genotype))
    MAC <- as.integer(pmin(alt_AC, 2 * N - alt_AC))
    genotype <- matrix_flip(genotype)
    MAF <- genotype$MAF
    variant_label <- as.vector(MAF>0)
    Geno <- genotype$Geno[,variant_label,drop=FALSE]
    Geno <- as(Geno,"dgCMatrix")
    rm(genotype)
    gc()
    
    if(obj_nullmodel$relatedness){
      if(!obj_nullmodel$sparse_kins){
        stop(paste0("Please use a sparse genetic relatedness matrix when fitting the null model!"))
      }
    }else{
      if(obj_nullmodel$family[1] == "binomial"){
        obj_nullmodel$Sigma_i <- Diagonal(x = obj_nullmodel$weights)
      }else if(obj_nullmodel$family[1] == "gaussian"){
        obj_nullmodel$Sigma_i <- Diagonal(length(obj_nullmodel$y)) / summary(obj_nullmodel)$dispersion
      }
      obj_nullmodel$Sigma_iX <- obj_nullmodel$Sigma_i %*% model.matrix(obj_nullmodel)
      obj_nullmodel$cov <- solve(t(model.matrix(obj_nullmodel)) %*% obj_nullmodel$Sigma_i %*% model.matrix(obj_nullmodel))
      obj_nullmodel$scaled.residuals <- (obj_nullmodel$y - obj_nullmodel$fitted.values) / summary(obj_nullmodel)$dispersion
    }
    
    GTSinvX_cov <- as(t(Geno) %*% obj_nullmodel$Sigma_iX,"matrix") %*% sqrtm(obj_nullmodel$cov)
    #V <- diag(t(Geno) %*% obj_nullmodel$Sigma_i %*% Geno - GTSinvX_cov %*% t(GTSinvX_cov))
    V <- colSums(Geno * (obj_nullmodel$Sigma_i %*% Geno)) - rowSums(GTSinvX_cov^2) # faster for large number of variants
    U <- as.vector(t(Geno) %*% obj_nullmodel$scaled.residuals)
    rm(Geno)
    gc()
    
    if (!is.null(qc_label)){
      sumstat <- data.frame(variant_info[,c("chr","pos","ref","alt")],qc_label,alt_AC,MAC,MAF,N)
    }else{
      sumstat <- data.frame(variant_info[,c("chr","pos","ref","alt")],alt_AC,MAC,MAF,N)
    }
    if (!is.null(annotation_phred)){
      sumstat <- cbind(sumstat,annotation_phred)
    }
    sumstat <- sumstat[variant_label,]
    if(!for_individual_analysis){
      sumstat <- cbind(sumstat,U,V,GTSinvX_cov)
    }else{
      sumstat <- cbind(sumstat,U,V)
    }
    
    return(sumstat)
  }
  
  MetaSTAARlite_worker_cov <- function(genotype,obj_nullmodel,cov_maf_cutoff=0.05,
                                       qc_label=NULL,signif.digits=3){
    
    if(!inherits(genotype, "matrix") && !inherits(genotype, "Matrix")){
      stop("genotype is not a matrix!")
    }
    
    if(cov_maf_cutoff < 0 | cov_maf_cutoff > 0.5){
      stop("cov_maf_cutoff should be a number between 0 and 0.5!")
    }
    
    if(cov_maf_cutoff == 0.5){
      cov_maf_cutoff <- 0.5 + 1e-16
    }
    
    if(inherits(genotype, "sparseMatrix")){
      MAF <- colMeans(genotype)/2
      if (!is.null(qc_label)){
        RV_label <- as.vector((MAF<cov_maf_cutoff)&(MAF>0)&(qc_label=="PASS"))
      }else{
        RV_label <- as.vector((MAF<cov_maf_cutoff)&(MAF>0))
      }
      Geno_rare <- genotype[,RV_label,drop=FALSE]
      Geno_rare <- as(Geno_rare,"dgCMatrix")
    }else{
      genotype <- matrix_flip(genotype)
      MAF <- genotype$MAF
      if (!is.null(qc_label)){
        RV_label <- as.vector((MAF<cov_maf_cutoff)&(MAF>0)&(qc_label=="PASS"))
      }else{
        RV_label <- as.vector((MAF<cov_maf_cutoff)&(MAF>0))
      }
      Geno_rare <- genotype$Geno[,RV_label,drop=FALSE]
      Geno_rare <- as(Geno_rare,"dgCMatrix")
    }
    
    rm(genotype,MAF)
    gc()
    
    if(obj_nullmodel$relatedness){
      if(!obj_nullmodel$sparse_kins){
        stop(paste0("Please use a sparse genetic relatedness matrix when fitting the null model!"))
      }
    }else{
      if(obj_nullmodel$family[1] == "binomial"){
        obj_nullmodel$Sigma_i <- Diagonal(x = obj_nullmodel$weights)
      }else if(obj_nullmodel$family[1] == "gaussian"){
        obj_nullmodel$Sigma_i <- Diagonal(length(obj_nullmodel$y)) / summary(obj_nullmodel)$dispersion
      }
      obj_nullmodel$Sigma_iX <- obj_nullmodel$Sigma_i %*% model.matrix(obj_nullmodel)
      obj_nullmodel$cov <- solve(t(model.matrix(obj_nullmodel)) %*% obj_nullmodel$Sigma_i %*% model.matrix(obj_nullmodel))
      obj_nullmodel$scaled.residuals <- (obj_nullmodel$y - obj_nullmodel$fitted.values) / summary(obj_nullmodel)$dispersion
    }
    
    GTSinvG_rare <- t(obj_nullmodel$Sigma_i %*% Geno_rare) %*% Geno_rare
    rm(Geno_rare)
    gc()
    
    GTSinvG_rare <- as(GTSinvG_rare,"TsparseMatrix")
    remove_ind <- (GTSinvG_rare@j < GTSinvG_rare@i) # Be careful about multi-allelic issue
    GTSinvG_rare@i <- GTSinvG_rare@i[!remove_ind]
    GTSinvG_rare@j <- GTSinvG_rare@j[!remove_ind]
    GTSinvG_rare@x <- GTSinvG_rare@x[!remove_ind]
    rm(remove_ind)
    gc()
    GTSinvG_rare <- as(GTSinvG_rare,"CsparseMatrix")
    
    ### Create a version of GTSinvG_rare with rounded significant digits
    if(!is.null(signif.digits)){
      GTSinvG_rare <- signif(GTSinvG_rare, digits = signif.digits)
    }
    
    return(GTSinvG_rare)
  }
  
  ncRNA_MetaSTAARlite_worker <- function(chr,gene_name,genofile,obj_nullmodel,known_loci=NULL,
                                         cov_maf_cutoff=0.05,signif.digits=NULL,
                                         QC_label="annotation/filter",check_qc_label=FALSE,variant_type=c("SNV","Indel","variant"),
                                         Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                         Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                                         silent=FALSE){
    
    ## evaluate choices
    variant_type <- match.arg(variant_type)
    
    phenotype.id <- as.character(obj_nullmodel$id_include)
    
    filter <- seqGetData(genofile, QC_label)
    if(variant_type=="variant")
    {
      if(check_qc_label)
      {
        SNVlist <- TRUE
      }else
      {
        SNVlist <- filter == "PASS"
      }
    }
    
    if(variant_type=="SNV")
    {
      if(check_qc_label)
      {
        SNVlist <- isSNV(genofile)
      }else
      {
        SNVlist <- (filter == "PASS") & isSNV(genofile)
      }
    }
    
    if(variant_type=="Indel")
    {
      if(check_qc_label)
      {
        SNVlist <- !isSNV(genofile)
      }else
      {
        SNVlist <- (filter == "PASS") & (!isSNV(genofile))
      }
    }
    
    variant.id <- seqGetData(genofile, "variant.id")
    
    if(!is.null(known_loci))
    {
      position <- as.integer(seqGetData(genofile, "position"))
      allele <- seqGetData(genofile, "allele")
      
      loc_SNV <- c()
      for (i in 1:dim(known_loci)[1]){
        loc_SNV <- c(loc_SNV,which((position==known_loci$POS[i])&(allele==paste0(known_loci$REF[i],",",known_loci$ALT[i]))))
      }
      
      seqSetFilter(genofile,variant.id=variant.id[loc_SNV],sample.id=phenotype.id)
      
      G_SNV <- seqGetData(genofile, "$dosage")
      
      if (!is.null(G_SNV)){
        id.SNV <- seqGetData(genofile,"sample.id")
        id.SNV.match <- rep(0,length(id.SNV))
        
        for(i in 1:length(id.SNV))
        {
          id.SNV.match[i] <- which.max(id.SNV==phenotype.id[i])
        }
        
        G_SNV <- G_SNV[id.SNV.match,,drop=FALSE]
        G_SNV_MAF <- apply(G_SNV, 2, function(x){sum(x[!is.na(x)])/2/sum(!is.na(x))})
        for (i in 1:length(G_SNV_MAF)){
          if (G_SNV_MAF[i] <= 0.5){
            G_SNV[is.na(G_SNV[,i]),i] <- 0
          }
          else{
            G_SNV[is.na(G_SNV[,i]),i] <- 2
          }
        }
        
        pos_adj <- as.integer(seqGetData(genofile, "position"))
        ref_adj <- as.character(seqGetData(genofile, "$ref"))
        alt_adj <- as.character(seqGetData(genofile, "$alt"))
        variant_adj_info <- data.frame(chr,pos_adj,ref_adj,alt_adj)
        colnames(variant_adj_info) <- c("chr","pos","ref","alt")
        variant_adj_info
        
        seqResetFilter(genofile)
      }
    }
    
    rm(filter)
    gc()
    
    ## ncRNA SNVs
    GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
    is.in <- ((GENCODE.Category=="ncRNA_exonic")|(GENCODE.Category=="ncRNA_exonic;splicing")|(GENCODE.Category=="ncRNA_splicing"))&(SNVlist)
    
    variant.id.ncRNA <- variant.id[is.in]
    
    rm(GENCODE.Category)
    gc()
    
    seqSetFilter(genofile,variant.id=variant.id.ncRNA,sample.id=phenotype.id)
    
    rm(variant.id.ncRNA)
    gc()
    
    GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
    GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[;]")
    Gene <- as.character(sapply(GENCODE.Info.split,function(z) gsub("\\(.*\\)","",z[1])))
    
    Gene_list_1 <- as.character(sapply(strsplit(Gene,','),'[',1))
    Gene_list_2 <- as.character(sapply(strsplit(Gene,','),'[',2))
    Gene_list_3 <- as.character(sapply(strsplit(Gene,','),'[',3))
    
    rm(GENCODE.Info)
    gc()
    
    rm(GENCODE.Info.split)
    gc()
    
    variant.id.ncRNA <- seqGetData(genofile, "variant.id")
    
    seqResetFilter(genofile)
    
    ### Gene
    is.in <- union(which(Gene_list_1==gene_name),which(Gene_list_2==gene_name))
    is.in <- union(is.in,which(Gene_list_3==gene_name))
    
    variant.is.in <- variant.id.ncRNA[is.in]
    
    seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)
    
    pos <- as.integer(seqGetData(genofile, "position"))
    ref <- as.character(seqGetData(genofile, "$ref"))
    alt <- as.character(seqGetData(genofile, "$alt"))
    if(check_qc_label){
      qc_label <- as.character(seqGetData(genofile, QC_label))
    }else{
      qc_label <- NULL
    }
    
    ## genotype id
    id.genotype <- seqGetData(genofile,"sample.id")
    # id.genotype.match <- rep(0,length(id.genotype))
    
    id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
    phenotype.id.merge <- data.frame(phenotype.id)
    phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
    id.genotype.match <- phenotype.id.merge$index
    
    ## Genotype
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,,drop=FALSE]
    
    summary_stat <- NULL
    GTSinvG_rare <- NULL
    cov_cond <- NULL
    
    if(!is.null(Geno))
    {
      # Summary statistics
      genotype <- matrix_impute(Geno)
      variant_info <- data.frame(chr,pos,ref,alt,row.names=NULL)
      
      ## Annotation
      Anno.Int.PHRED.sub <- NULL
      Anno.Int.PHRED.sub.name <- NULL
      
      if(variant_type=="SNV")
      {
        if(Use_annotation_weights)
        {
          for(k in 1:length(Annotation_name))
          {
            if(Annotation_name[k]%in%Annotation_name_catalog$name)
            {
              Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
              Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
              
              if(Annotation_name[k]=="CADD")
              {
                Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
              }
              Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
            }
          }
          
          Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
          colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
        }
      }
      
      try(summary_stat <- MetaSTAARlite_worker_sumstat(genotype,obj_nullmodel,variant_info,qc_label,
                                                       Anno.Int.PHRED.sub),silent=silent)
      
      # Covariance matrices
      genotype <- matrix_flip(Geno)$Geno
      genotype <- as(genotype,"dgCMatrix")
      
      rm(Geno)
      gc()
      
      try(GTSinvG_rare <- MetaSTAARlite_worker_cov(genotype,obj_nullmodel,cov_maf_cutoff,
                                                   qc_label,signif.digits),silent=silent)
      
      # Covariance matrices for conditional analysis
      if(!is.null(known_loci))
      {
        try(cov_cond <- MetaSTAARlite_worker_cov_cond(genotype,G_SNV,obj_nullmodel,variant_info,variant_adj_info),silent=silent)
      }
    }
    
    seqResetFilter(genofile)
    
    if(!is.null(known_loci))
    {
      return(list(summary_stat=summary_stat,
                  GTSinvG_rare=GTSinvG_rare,
                  cov_cond=cov_cond))
    }else
    {
      return(list(summary_stat=summary_stat,
                  GTSinvG_rare=GTSinvG_rare))
    }
  }
  
  ###########################################################
  #           User Input
  ###########################################################
  ## Null model
  obj_nullmodel <- get(load("obj.STAAR.UKB.TC.20240530.63370.Rdata"))
  
  ## QC_label
  QC_label <- "annotation/info/QC_label2"
  geno_missing_imputation <- "mean"
  variant_type <- "SNV"	
  
  ## Annotation_dir
  Annotation_dir <- "annotation/info/FunctionalAnnotation"
  ## Annotation channel
  Annotation_name_catalog <- read.csv("Annotation_name_catalog.csv")
  ## Use_annotation_weights
  Use_annotation_weights <- TRUE
  ## Annotation name
  Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                       "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")
  
  ## output file name
  output_file_name <- "UKB_TC_ncRNA"
  
  agds_dir <- paste0("ukb.200k.wgs.chr",1:22,".pass.annotated.gds")
  
  
  ###########################################################
  #           Main Function 
  ###########################################################
  arrayid <- c(117,218,220,220,221,156,219)
  sub_seq_id <- c(53,19,208,274,311,41,103)
  
  region_spec <- data.frame(arrayid,sub_seq_id)
  
  gene_num_in_array <- 100
  group.num.allchr <- ceiling(table(ncRNA_gene[,1])/gene_num_in_array)
  sum(group.num.allchr)
  
  ncRNA_sumstat <- list()
  ncRNA_cov <- list()
  for(kk in 1:dim(region_spec)[1])
  {
    arrayid <- region_spec$arrayid[kk]
    sub_seq_id <- region_spec$sub_seq_id[kk]
    
    chr <- which.max(arrayid <= cumsum(group.num.allchr))
    ncRNA_gene_chr <- ncRNA_gene[ncRNA_gene[,1]==chr,]
    
    ## aGDS file
    agds.path <- agds_dir[chr]
    genofile <- seqOpen(agds.path)
    
    gene_name <- ncRNA_gene_chr[sub_seq_id,2]
    results_temp <- ncRNA_MetaSTAARlite_worker(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                               cov_maf_cutoff=0.05,signif.digits=NULL,
                                               QC_label=QC_label,check_qc_label=TRUE,variant_type=variant_type,
                                               Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                               Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
    ncRNA_sumstat[[gene_name]] <- results_temp$summary_stat
    ncRNA_cov[[gene_name]] <- results_temp$GTSinvG_rare
    
    seqClose(genofile)
  }
  
  system("rm obj.STAAR.UKB.TC.20240226.Rdata")
  for(i in 1:length(agds_dir)){
    system(paste0("rm ",agds_dir[i])) 
  }
  system("rm Annotation_name_catalog.csv")
  
  save(ncRNA_sumstat,file=paste0(output_file_name,"_sumstat_",223,"_63370.Rdata"),compress = "xz")
  save(ncRNA_cov,file=paste0(output_file_name,"_cov_",223,"_63370.Rdata"),compress = "xz")
})[3]

save(time,file = paste0(output_file_name,"_time_",223,"_63370.Rdata"))