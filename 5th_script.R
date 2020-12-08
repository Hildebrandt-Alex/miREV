#--------------------------------------------------------#
# 10. Apply RefGene Algorithms and create a ranked table ####
#--------------------------------------------------------#
# ___prepare readcount tables for algorithms ####
#--------------------------------------------#
# replace all 0, inf, NaN with NA
counts_list <- list()
for (i in 1:length(count_list_normalized)) {
  counts_list_Na <- count_list_normalized[[i]]
  name <- names(count_list_normalized[i])
  for (j in 1:length(count_list_normalized[[i]])) {
    counts_list_Na[[j]][counts_list_Na[[j]]==0]         <- NA
    counts_list_Na[[j]][counts_list_Na[[j]]=="NaN"]     <- NA
    counts_list_Na[[j]][counts_list_Na[[j]]==Inf]       <- NA
  }
  counts_list[[name]] <- counts_list_Na
}

# transpose data for NormFinder and geNorm Algorithms
counts_list_trans <- list()
for (i in 1:length(counts_list)) {
  counts_list_t <- counts_list[[i]]
  name <- names(counts_list[i])
  for (j in 1:length(counts_list[[i]])) {
    counts_list_t[[j]]<- t(counts_list_t[[j]])
  }
  counts_list_trans[[name]] <- counts_list_t
}

#--------------------------------------#
# ___bestKeeper with calculation of CV ####
#--------------------------------------#
ref_gene_CV_list_complete <- list()
ref_gene_CV_complete <- list()
ref_gene_CV_list  <- list()
for (i in 1:length(counts_list_trans)) {
  ref_gene_CV <- counts_list_trans[[i]]
  name <- names(counts_list_trans[i])
  for (j in 1:length(counts_list_trans[[i]])) {
    
    ref_gene_CV[[j]]          <- apply(ref_gene_CV[[j]], 2, raster::cv, na.rm = TRUE)
    ref_gene_CV[[j]]          <- as.data.frame(ref_gene_CV[[j]])
    ref_gene_CV[[j]][,"gene"] <- rownames(ref_gene_CV[[j]])
    ref_gene_CV[[j]]          <- ref_gene_CV[[j]][order(ref_gene_CV[[j]][1]),]
    ref_gene_CV[[j]]$gene     <- make.names(ref_gene_CV[[j]]$gene)
    rownames(ref_gene_CV[[j]]) <- NULL
    rownames(ref_gene_CV[[j]]) <- ref_gene_CV[[j]]$gene
    # copy complete list
    norm_name                  <- names(ref_gene_CV[j])
    ref_gene_CV_complete[[norm_name]]  <- ref_gene_CV[[j]] 
    # apply treshhold for candidate genes (50% above under best stability value)
    ref_gene_CV[[j]]           <- ref_gene_CV[[j]][which(ref_gene_CV[[j]]$`ref_gene_CV[[j]]` <= min(ref_gene_CV[[j]]$`ref_gene_CV[[j]]`)*treshold_CV),]
    
  }
  ref_gene_CV_list[[name]] <- ref_gene_CV
  ref_gene_CV_list_complete[[name]]  <- ref_gene_CV_complete
}

#----------------------------------#
# ___geNorm with meanM calculation ####
#----------------------------------#
ref_gene_geNorm_list_complete <- list()
ref_gene_geNorm_complete  <- list()
ref_gene_geNorm_list <- list()
for (i in 1:length(counts_list_trans)) {
  ref_gene_geNorm <- counts_list_trans[[i]]
  name <- names(counts_list_trans[i])
  for (j in 1:length(counts_list_trans[[i]])) {
    gene_names                            <- dimnames(ref_gene_geNorm[[j]])
    gene_names                            <- gene_names[[2]]
    ref_gene_geNorm[[j]]                  <- NormqPCR::selectHKs(ref_gene_geNorm[[j]], method = "geNorm", log = FALSE, Symbols = gene_names, trace = FALSE, na.rm =TRUE) 
    ref_gene_geNorm[[j]]$meanM[length(ref_gene_geNorm[[j]]$meanM)+1] <- NA
    ref_gene_geNorm[[j]]                <- data.frame("gene" = ref_gene_geNorm[[j]]$ranking,
                                                      "meanM" = rev(ref_gene_geNorm[[j]]$meanM))
    row.names(ref_gene_geNorm[[j]])     <- NULL
    ref_gene_geNorm[[j]]$gene           <- make.names(ref_gene_geNorm[[j]]$gene)
    row.names(ref_gene_geNorm[[j]])     <- ref_gene_geNorm[[j]]$gene
    ref_gene_geNorm[[j]][1,2]           <- ref_gene_geNorm[[j]][2,2]
    # copy complete list
    norm_name                           <- names(ref_gene_geNorm[j])
    ref_gene_geNorm_complete[[norm_name]]  <- ref_gene_geNorm[[j]] 
    
    # apply treshhold for candidate genes (50% above under best stability value)
    ref_gene_geNorm[[j]]           <- ref_gene_geNorm[[j]][which(ref_gene_geNorm[[j]]$meanM <= min(ref_gene_geNorm[[j]]$meanM)*treshold_genorm),]
    
  }
  ref_gene_geNorm_list[[name]] <- ref_gene_geNorm
  ref_gene_geNorm_list_complete[[name]]  <- ref_gene_geNorm_complete
}
#------------------------------------#
# ___normFinder with roh calculation ####
#------------------------------------#
groups_normF_list <- list()
for (i in 1:length(counts_list_trans)) {
  groups_normF                          <- as.factor(split_list[[i]]$Disease)  
  name <- names(counts_list_trans[i])
  # split patients in two groups, if not assigned by split_list ("0" and "1" in split_list$Diesease )
  if (nlevels(groups_normF) == 1) {
    v <- rep(c("group1", "group2"),1000)
    groups_normF <- v[1:length(groups_normF)]
  }else {NULL}
  groups_normF_list[[name]] <- groups_normF
}

ref_gene_normF_list_complete <- list()
ref_gene_normF_complete <- list()
ref_gene_normF_list <- list()
for (i in 1:length(counts_list_trans)) {
  ref_gene_normF <- counts_list_trans[[i]]
  name <- names(counts_list_trans[i])
  for (j in 1:length(counts_list_trans[[i]])) {
    gene_names                            <- dimnames(ref_gene_normF[[j]])
    gene_names                            <- gene_names[[2]]
    groups                                <- unlist(groups_normF_list[i])
    ref_gene_normF[[j]]                  <- NormqPCR::selectHKs(ref_gene_normF[[j]], method = "NormFinder", log = FALSE, Symbols = gene_names, trace = FALSE, group = groups, minNrHKs = length(gene_names)-1 , na.rm =TRUE) 
    ref_gene_normF[[j]]                  <- data.frame("gene" = ref_gene_normF[[j]]$ranking,
                                                       "stab_measure_rho" = ref_gene_normF[[j]]$rho)
    row.names(ref_gene_normF[[j]]) <- make.names(ref_gene_normF[[j]]$gene)
    ref_gene_normF[[j]]            <- ref_gene_normF[[j]][order(ref_gene_normF[[j]]$stab_measure_rho, decreasing = TRUE),]
    # copy complete list
    norm_name                             <- names(ref_gene_normF[j])
    ref_gene_normF_complete[[norm_name]]  <- ref_gene_normF[[j]]
    ref_gene_normF_complete[[j]]$gene     <- make.names(as.character(ref_gene_normF_complete[[j]]$gene))
    
    # apply treshhold for candidate genes (50% above under best stability value)
    ref_gene_normF[[j]]           <- ref_gene_normF[[j]][which(ref_gene_normF[[j]]$stab_measure_rho >= max(ref_gene_normF[[j]]$stab_measure_rho)*treshold_NormF),]
    ref_gene_normF[[j]]$gene      <- make.names(as.character(ref_gene_normF[[j]]$gene))
  }
  ref_gene_normF_list[[name]] <- ref_gene_normF
  ref_gene_normF_list_complete[[name]]  <- ref_gene_normF_complete
}
