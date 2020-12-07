### ---------------------------------------------------------------------------------------------------------###
### Automated finding of appropriate reference miRNAs in blood EVs readcount data sets (Alex Hildebrandt)
###----------------------------------------------------------------------------------------------------------###

# Description of the script ####

### tasks of the script:
### 
### 1. load necessary libraries 
### 2. load readcount table and metadata of interest - one readcount file for all samples and one metadata file for all samples 
### 3. VALIDATION: check if readcount table and metadata have same order - lines with code to setcolorder are included in this part, and need to be executed if necessary
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
### 4. PARAMETERS: define all parameters which are necessary to run the script
#                  1: Filtering samples (% of mapped miRNA in ration to all readcounts)                    (-> 7% in first calculation)
#                  2: Filtering genes   (% of obligatory appearance of gene across all samples)            (-> 95% in first calculation)
#                  3: split combinations for readcount table
#                  4: treshhold CV      (% of genes under best ranked, which should be taken into account) (-> 5% above best stability value in first calculation)
#                  5: treshhold geNomr  (% of genes under best ranked, which should be taken into account) (-> 5% above best stability value in first calculation)
#                  6: treshhold NormF   (% of genes under best ranked, which should be taken into account) (-> 5% under best stability value in first calculation)
#                  7: treshhold sample  (number of samples, which have to appear in each list - lists with less samples will be removed -> necessary for NormF: 10 samples are recommended)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
### 5. FILTERING : samples of bad quality are removed - according to "Filtering samples" - PARAMETER
### 6. VALIDATION: check if readcount table and metadata have same order after filtering
### 7. TASK      : Produce splitted readcounds in a list according to "split combination"- PARAMETER                         --> count_list()
### 8. FILTERING : remove genes, which are under the treshholdd of: "Filtering genes" - PARAMETER of the samples             --> count_list_filtered()  
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
### 9. TASK      : NORMALIZATION: six different normalization Methods are used for each list:                                --> count_list_normalized()
#                                 --> therefore a selfmade normalization function is written (filename: FUNCTION_1_norm_fun ; you can call it with "norm_function(count_matrix_raw, metdata)")  
#                  1:TC         - total count,
#                  2:Med        - median,
#                  3:Q          - quantile,
#                  4:UQ         - upper quartile,
#                  5:TMM        - trimmed mean of M-values,
#                  6:DESeq      - median of expression ratios
#                  7:raw_counts - raw counts are as well in the list 
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
### 10. TASK     : RefGeneFinder Algorithms are applied to each normalized list of each combination list
#                  1: BestKeeper - CV is calculated for each gene                                                   --> ref_gene_CV_list()             - with "treshhold CV"-PARAMETER
#                                                                                                                   --> ref_gene_CV_list_complete()    - without "treshhold CV"-PARAMETER
#                  2: NormFinder - Average expression stability value rho is calculated for each gene               --> ref_gene_normF_list()          - with "treshhold normF"-PARAMETER
#                                                                                                                   --> ref_gene_normF_list()          - without "treshhold normF"-PARAMETER  
#                  3: geNorm     - Average expression stability value (M) is calculated for each gene               --> ref_gene_geNorm_list()         - with "treshhold geNorm"-PARAMETER                                
#                                                                                                                   --> ref_gene_geNorm_list()         - without "treshhold geNorm"-PARAMETER      
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
### 11. TASK     : EVALUATION Tasks:
#                  1: TASK       : create gene lists for each RefGeneFinder Algorithms with and without stability values:
#                                                 --> CV_list()       -> without stability values
#                                                 --> CV_list2()      -> with    stability values
#                                                 --> geNorm_list()   -> without stability values
#                                                 --> geNorm_list2()  -> with    stability values
#                                                 --> normF_list()    -> without stability values
#                                                 --> normF_list2()   -> with    stability values
#                  2: TASK          :merge lits:                                                                     --> RefGene_list()  -> without stability values
#                                                                                                                    --> RefGene_list2() -> with    stability values
#                  3: PARAMETERS: SET COMBO of interest to be anlysed via Overlap analysis; Currently the following parameters are selectable:
#                                         - RefGene Algorithm (all or some of them)
#                                         - Study (all datasets also possible)
#                                         - kind of disease (which disease are you interested in)
#                                         - kind of EV-preparation type (which method was used)
#                                         - normalization method (one or some of them, which method is used by the user of the app)
#                  4: TASK      : create shiny list with COMBO of interest: therefore run only lines which you have changed (if you are interested in a disease just run first line, which is necessary and line with disease parameter)
#                  5: VALIDATION: check names of shiny list - should include choosen names of choosen parameters
#                  6: TASK      : create 3x18 data frame with stability values  (refGeneAlg x normalization methods)
#                  7: TASK      : create 1x6 data frame  (normalization methods irrespective of RefGeneAlg)
#                  
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#


# 1.load libraries ####
library(data.table)
# for normalization function
#BiocManager::install(c("edgeR", "preprocessCore", "raster"))
library(edgeR)   
library(graphics)
library(stats)
library(preprocessCore) 
library(DESeq2)
library(raster) 
library(matrixStats)
library(data.table)
# for RefGene Algorithms and extraction & ranking of genes
#BiocManager::install(c("ctrlGene", "NormqPCR", "rgl"))
library(ctrlGene)
library(NormqPCR)
library(ggplot2)
library(raster)
library(rgl)
library(plot3Drgl)
library(data.table)####
# for Overlap analysis
#BiocManager::install(c("veccompare", "rowr", "plot3Drgl"))
library(veccompare)
library(rowr)
library(tidyr)
library(VennDiagram)
library(gridExtra)
library(grid)
library(plyr)
library(stringr)
library(gimme)

#--------------------------#
# save/load workspace      #
#--------------------------#
save.image(file = "WORKSPACE_654_samples.RData")
load("Worksace654_samples_EV_TRACK_variables.RData")
#----------------------------------------#
# 2. load readcount table and metadata   #####
#----------------------------------------#
# set working directory
setwd("~/hard_drive/usrs_folder/Alex_1984/1_LRZ_Rstudio_projects_sync_folder/1_LRZ_Rstudio_projects_sync_folder/1_miRNA_Ref_Gene/")
# load files
load(file = "./data/metaDB_654_EV_track.Rda")                        # -> metaDB_both_batches
load(file = "./data/parent_miRNA_1mm_both_batches_ordered.Rda")              # -> parent_both_batches

# replace by loaded data to use for the pipeline
metaDB_both_batches <- as.data.frame(metaDB_both_batches_EV_Track)
parent_both_batches <- parent_both_batches
# setcolorder metaDB and countmatrix
gene_names <- rownames(parent_both_batches)
colnames(parent_both_batches) <- make.names(colnames(parent_both_batches))
metaDB_both_batches$sample.name <- make.names(metaDB_both_batches$sample.name)
parent_both_batches <- as.data.table(parent_both_batches)
setcolorder(parent_both_batches, as.character(metaDB_both_batches$sample.name))
parent_both_batches <- as.data.frame(parent_both_batches)
rownames(parent_both_batches) <- gene_names

#----------------------------------------------#
# 3.check if both sheets are in the same order ####
#----------------------------------------------#
table(colnames(parent_both_batches)==metaDB_both_batches$sample.name)    # --> must be TRUE in all observations of metaDB   


#------------------------#
# 4. define Parameters   ####
#------------------------#
# Filter sample: define % of mapping filter (samples to remove) and gene filter (gene_filter <- 0.05 => keep genes which have counts in 95% of all samples)
mapped_miRNA <- 7
# Filter genes:  (inverse percentage have to be assigned) define how many samples at least should have the gene included (here 95% of all samples need to have the gene, otherwise it will be removed)
gene_filter  <- 0.05
# treshold for candidate genes (50% under above best stability value)
treshold_genorm <- 1.5
treshold_NormF  <- 0.5
treshold_CV     <- 1.5
# treshold to remove lists less than sample_tresh
sample_tresh    <- 10

#------------------------#
# 5. filtering samples   #####
#------------------------#
# ___1.1 filtering samples: remove samples with less mapped miRNAs of total library size than 7%  
table(metaDB_both_batches$`%_mapped_miRNA_to_lib.size` >mapped_miRNA)
keep.samples <- metaDB_both_batches$`%_mapped_miRNA_to_lib.size` >mapped_miRNA

parent_both_batches_1_filtered <- parent_both_batches[,keep.samples]
metaDB_both_batches_1_filtered <- metaDB_both_batches[keep.samples,]

#-----------------------------------------------#
# 6. check if both sheets are in the same order ####
#-----------------------------------------------#
table(colnames(parent_both_batches_1_filtered)==metaDB_both_batches_1_filtered$sample.name)    # must be TRUE in all observations of metaDB   



#----------------------------------------------------------------#
# 7. splitting counts and metaDB according to pheno-combinations #####
#----------------------------------------------------------------#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! use combn() to determine unique combinations 
#___1st level of kombination (with one variable) ####
# _______split disease #####
#first make a list of all possible comnination in one factor
number_levels <- length(levels(as.factor(metaDB_both_batches_1_filtered$kind_of_disease)))
temp_1 <- list()
for (i in 1:number_levels) {
  temp_1[[i]]  <- combn(levels(as.factor(metaDB_both_batches_1_filtered$kind_of_disease)),i)
}
#second split/merge count tables in list
split_list_disease <- list()
for (i in 1:length(temp_1)) {
  dd_sub <- temp_1[[i]]
  dd_sub <- as.data.frame(dd_sub)
  length_dd <- length(dd_sub)
  
  for (j in 1:length_dd) {
    name <- paste(temp_1[[i]][,j], collapse = ".")
    split_list_disease[[name]] <- metaDB_both_batches_1_filtered[metaDB_both_batches_1_filtered$kind_of_disease %in% as.character(temp_1[[i]][,j]),] 
  }
}
# _______split EV_isolation method #####-------------------------------------####
#first make a list of all possible comnination in one factor
number_levels <- length(levels(as.factor(metaDB_both_batches_1_filtered$EV_Isolation_Method_EV_track)))
temp_1 <- list()
for (i in 1:number_levels) {
  temp_1[[i]]  <- combn(levels(as.factor(metaDB_both_batches_1_filtered$EV_Isolation_Method_EV_track)),i)
}
#second split/merge count tables in list
split_list_ev_iso <- list()
for (i in 1:length(temp_1)) {
  dd_sub <- temp_1[[i]]
  dd_sub <- as.data.frame(dd_sub)
  length_dd <- length(dd_sub)
  
  for (j in 1:length_dd) {
    name <- paste(temp_1[[i]][,j], collapse = ".")
    split_list_ev_iso[[name]] <- metaDB_both_batches_1_filtered[metaDB_both_batches_1_filtered$EV_Isolation_Method_EV_track %in% as.character(temp_1[[i]][,j]),] 
  }
}
# _______split species #####----------------------------------------####
#first make a list of all possible comnination in one factor
number_levels <- length(levels(as.factor(metaDB_both_batches_1_filtered$species)))
temp_1 <- list()
for (i in 1:number_levels) {
  temp_1[[i]]  <- combn(levels(as.factor(metaDB_both_batches_1_filtered$species)),i)
}
#second split/merge count tables in list
split_list_species <- list()
for (i in 1:length(temp_1)) {
  dd_sub <- temp_1[[i]]
  dd_sub <- as.data.frame(dd_sub)
  length_dd <- length(dd_sub)
  
  for (j in 1:length_dd) {
    name <- paste(temp_1[[i]][,j], collapse = ".")
    split_list_species[[name]] <- metaDB_both_batches_1_filtered[metaDB_both_batches_1_filtered$species %in% as.character(temp_1[[i]][,j]),] 
  }
}
# _______split biofluids #####----------------------------------------####
#first make a list of all possible comnination in one factor
number_levels <- length(levels(as.factor(metaDB_both_batches_1_filtered$biofluids)))
temp_1 <- list()
for (i in 1:number_levels) {
  temp_1[[i]]  <- combn(levels(as.factor(metaDB_both_batches_1_filtered$biofluids)),i)
}
#second split/merge count tables in list
split_list_biofluids <- list()
for (i in 1:length(temp_1)) {
  dd_sub <- temp_1[[i]]
  dd_sub <- as.data.frame(dd_sub)
  length_dd <- length(dd_sub)
  
  for (j in 1:length_dd) {
    name <- paste(temp_1[[i]][,j], collapse = ".")
    split_list_biofluids[[name]] <- metaDB_both_batches_1_filtered[metaDB_both_batches_1_filtered$biofluids %in% as.character(temp_1[[i]][,j]),] 
  }
}
#---------------------------------------------------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------#
split_species <- split_list_species
split_biofluids <- split_list_biofluids
split_ev_iso <- split_list_ev_iso
split_disease <- split_list_disease

# ___2nd level of kombination (with two variables) ####
list_kombos <- list()
for (i in 1:length(split_species)) {
   #species & biofluids all kombinations with appearance 
 for (j in 1:length(split_biofluids)) {
    name_ele2 <- paste(names(split_species[i]), names(split_biofluids[j]), sep = ".")
    list_kombos[[name_ele2]] <- fintersect(as.data.table(split_species[[i]]),as.data.table(split_biofluids[[j]]))
  }
 # species $ iso all kombinations with appearance
 for (j in 1:length(split_ev_iso)) {
   name_ele3 <- paste(names(split_species[i]), names(split_ev_iso[j]), sep = ".")
   list_kombos[[name_ele3]] <- fintersect(as.data.table(split_species[[i]]),as.data.table(split_ev_iso[[j]]))
 }
 # species & disease
 for (j in 1:length(split_disease)) {
   name_ele4 <- paste(names(split_species[i]), names(split_disease[j]), sep = ".")
   list_kombos[[name_ele4]] <- fintersect(as.data.table(split_species[[i]]),as.data.table(split_disease[[j]]))
 }
 for (j in 1:length(split_biofluids)) {
   name_ele5 <- paste(names(split_species[i]), names(split_biofluids[j]), sep = ".")
   list_kombos[[name_ele5]] <- fintersect(as.data.table(split_species[[i]]),as.data.table(split_biofluids[[j]]))
 }
}
#------------------------------------------------------------------------------------#
list_kombos2 <- list()
for (i in 1:length(split_biofluids)) {
  #biofluids & isolation all kombinations with appearance 
  for (j in 1:length(split_ev_iso)) {
    name_ele2 <- paste(names(split_biofluids[i]), names(split_ev_iso[j]), sep = ".")
    list_kombos2[[name_ele2]] <- fintersect(as.data.table(split_biofluids[[i]]),as.data.table(split_ev_iso[[j]]))
  }
  #biofluids & disease all kombinations with appearance 
  for (j in 1:length(split_disease)) {
    name_ele3 <- paste(names(split_biofluids[i]), names(split_disease[j]), sep = ".")
    list_kombos2[[name_ele3]] <- fintersect(as.data.table(split_biofluids[[i]]),as.data.table(split_disease[[j]]))
  }
}
#--------------------------------------------------------------------------------------#
list_kombos3 <- list()
for (i in 1:length(split_ev_iso)) {
 #isolation & disease all kombinations with appearance 
  for (j in 1:length(split_disease)) {
    name_ele2 <- paste(names(split_ev_iso[i]), names(split_disease[j]), sep = ".")
    list_kombos3[[name_ele2]] <- fintersect(as.data.table(split_ev_iso[[i]]),as.data.table(split_disease[[j]]))
  }
}
#-------------------------------------------------------------------------------------------------#
# ___3nd level of kombination (with three variables) ####
list_kombos4 <- list()
for (i in 1:length(split_species)) {
  #species & biofluids & isolation all kombinations with appearance 
  for (j in 1:length(split_biofluids)) {
    for (k in 1:length(split_ev_iso)) {
      name_ele2 <- paste(names(split_species[i]), names(split_biofluids[j]), names(split_ev_iso[k]), sep = ".")
      temp_2                    <- fintersect(as.data.table(split_species[[i]]),as.data.table(split_biofluids[[j]]))
      list_kombos4[[name_ele2]] <- fintersect(as.data.table(temp_2),as.data.table(split_ev_iso[[k]]))
    }
  }
  #species & biofluids & disease all kombinations with appearance 
  for (j in 1:length(split_biofluids)) {
    for (k in 1:length(split_disease)) {
      name_ele2 <- paste(names(split_species[i]), names(split_biofluids[j]), names(split_disease[k]), sep = ".")
      temp_2                    <- fintersect(as.data.table(split_species[[i]]),as.data.table(split_biofluids[[j]]))
      list_kombos4[[name_ele2]] <- fintersect(as.data.table(temp_2),as.data.table(split_disease[[k]]))
    }
  }
  #species & isolation & disease all kombinations with appearance 
  for (j in 1:length(split_ev_iso)) {
    for (k in 1:length(split_disease)) {
      name_ele2 <- paste(names(split_species[i]), names(split_ev_iso[j]), names(split_disease[k]), sep = ".")
      temp_2                    <- fintersect(as.data.table(split_species[[i]]),as.data.table(split_ev_iso[[j]]))
      list_kombos4[[name_ele2]] <- fintersect(as.data.table(temp_2),as.data.table(split_disease[[k]]))
    }
  }
}
#---------------------------------------------------------------------------------------------------------------------------------#
# ___4th level of kombination (with four variables) ####
list_kombos5 <- list()
for (i in 1:length(split_species)) {
  #species & biofluids & isolation all kombinations with appearance 
  for (j in 1:length(split_biofluids)) {
    for (k in 1:length(split_ev_iso)) {
       for (l in 1:length(split_disease)) {
         name_ele2 <- paste(names(split_species[i]), names(split_biofluids[j]), names(split_ev_iso[k]),names(split_disease[l]), sep = ".")
         temp_2                    <- fintersect(as.data.table(split_species[[i]]),as.data.table(split_biofluids[[j]]))
         temp_3                    <- fintersect(as.data.table(temp_2),as.data.table(split_ev_iso[[k]]))
         list_kombos5[[name_ele2]] <- fintersect(as.data.table(temp_3),as.data.table(split_disease[[l]]))
      }
    }
  }
}

#--------------------------------------------------------------------------------------------#
# ___7.2 concatenate lists in one list ####                     
split_list <- append(split_list_species, split_list_ev_iso)
split_list <- append(split_list, split_list_biofluids)
split_list <- append(split_list, split_list_disease)
split_list <- append(split_list, list_kombos)
split_list <- append(split_list, list_kombos2)
split_list <- append(split_list, list_kombos3)
split_list <- append(split_list, list_kombos4)
split_list <- append(split_list, list_kombos5)

# ___7.3 filter split list: remove lists with less than 10 samples ####
remove_obj  <- list()
for (i in 1:length(split_list)) {
  name <- names(split_list[i])
  remove_obj[i] <- length(split_list[[i]]$sample.name)>=sample_tresh 
  
}
split_list <- split_list[unlist(remove_obj)]

#split  list innto 10 sublists
split_list_1 <- split_list[1:3000]
split_list_2 <- split_list[3001:6000]
split_list_3 <- split_list[6001:9000]
split_list_4 <- split_list[9001:12000]
split_list_5 <- split_list[12001:15000]
split_list_6 <- split_list[15001:18000]
split_list_7 <- split_list[18001:21000]
split_list_8 <- split_list[21001:24000]
split_list_9 <- split_list[24001:27000]
split_list_10 <- split_list[27001:29591]

split_list <- split_list_2
# ___7.4 split readcount table according to split_list ####                   
    count_list <- list()
    for (i in 1:length(split_list)) {
      
      name <- names(split_list[i])
      count_list[[name]] <- parent_both_batches_1_filtered[, unlist(split_list[[i]][,2])]  
    }
    
    
    #--------------------#
    # 8. filtering genes ####
    #--------------------#
    # ___3.1 fitering genes: Remove rows which assign more "NA's" than "Filtering genes" - PARAMETER  
    count_list_filtered <- list()
    for (i in 1:length(count_list)) {
      # convert 0 to NA
      count_list[[i]][count_list[[i]]==0] <- NA
      name <- names(count_list[i])
      count_list_filtered[[name]] <- count_list[[i]][-which(rowMeans(is.na(count_list[[i]]))>gene_filter),]
      #convert NA to 0
      count_list_filtered[[i]][is.na(count_list_filtered[[i]])] <- 0
    }
    # check dimension of each count table
    for (i in 1:length(count_list)) {b <- dim(count_list_filtered[[i]])
    print(b)
    }
    #-----------------------------------------------------------#
    # 9. Normalization of each lists readcount table seperately ####
    #-----------------------------------------------------------#
    # load normalization function
    source("./R_code/FUNCTION_1_norm_fun.R")
    # run function
    count_list_normalized <- list()
    
    for (i in 1:length(count_list_filtered)) {
      
      countmatrix <- count_list_filtered[[i]]
      metadata   <- split_list[[i]]
      name       <- names(count_list[i])
      
      count_list_normalized[[name]] <- norm_function(count_matrix_raw = countmatrix,
                                                     metadata = metadata)
    }
    
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
    #--------------------------------------------------------------------------------------------------------------------------------------------#
    ###----------------###
    ############ 11. EVALUATION ##############
    ###----------------###
    #--------------------------------------------------------------------------------------------------------------------------------------------#
    
    #--------------------------------------------------------------------------------------------------------------------------------------------#
    # ____1.create gene lists of each RefGeneAlgorithm (each list has the genes of .... combination multiplied with 7 normalization methods -> ...x7=??? elements in the list) ####
    #-------------------------------------------------------------------------------------------------------------------------------------------------------------#
    # CV list - create a vector list with all genes from all tables
    CV_list    <- list()
    for (i in 1:length(ref_gene_CV_list)) {
      temp_list <- ref_gene_CV_list[[i]]
      for (j in 1:length(ref_gene_CV_list[[j]])) {
        name <- paste("CV",names(ref_gene_CV_list[i]), names(temp_list[j]),sep="__") 
        CV_list[[name]] <- temp_list[[j]]$gene
      }
    }
    # CV list - create a vector list with all genes and stability values
    CV_list2    <- list()
    for (i in 1:length(ref_gene_CV_list)) {
      temp_list <- ref_gene_CV_list[[i]]
      for (j in 1:length(ref_gene_CV_list[[j]])) {
        name <- paste("CV",names(ref_gene_CV_list[i]), names(temp_list[j]),sep="__") 
        CV_list2[[name]] <- temp_list[[j]]
      }
    }
    # geNorm list - create a vector list with all genes from all tables
    geNorm_list    <- list()
    for (i in 1:length(ref_gene_geNorm_list)) {
      temp_list <- ref_gene_geNorm_list[[i]]
      for (j in 1:length(ref_gene_geNorm_list[[j]])) {
        name <- paste("geNorm",names(ref_gene_geNorm_list[i]), names(temp_list[j]),sep="__") 
        geNorm_list[[name]] <- temp_list[[j]]$gene
      }
    }
    # geNorm list - create a vector list with all genes and stability values
    geNorm_list2    <- list()
    for (i in 1:length(ref_gene_geNorm_list)) {
      temp_list <- ref_gene_geNorm_list[[i]]
      for (j in 1:length(ref_gene_geNorm_list[[j]])) {
        name <- paste("geNorm",names(ref_gene_geNorm_list[i]), names(temp_list[j]),sep="__") 
        geNorm_list2[[name]] <- temp_list[[j]]
      }
    }
    # normF list - create a vector list with all genes from all tables
    normF_list    <- list()
    for (i in 1:length(ref_gene_normF_list)) {
      temp_list <- ref_gene_normF_list[[i]]
      for (j in 1:length(ref_gene_normF_list[[j]])) {
        name <- paste("normF",names(ref_gene_normF_list[i]), names(temp_list[j]),sep="__") 
        normF_list[[name]] <- temp_list[[j]]$gene
      }
    }
    # normF list - create a vector list with all genes and stability values
    normF_list2    <- list()
    for (i in 1:length(ref_gene_normF_list)) {
      temp_list <- ref_gene_normF_list[[i]]
      for (j in 1:length(ref_gene_normF_list[[j]])) {
        name <- paste("normF",names(ref_gene_normF_list[i]), names(temp_list[j]),sep="__") 
        normF_list2[[name]] <- temp_list[[j]]
      }
    }
    
    #-----------------#
    # ____2.merge lists ####  
    #-----------------#
    RefGene_list2 <- append(CV_list, geNorm_list)
    RefGene_list2 <- append(RefGene_list, normF_list)
    
    RefGene_list2_2 <- append(CV_list2, geNorm_list2)
    RefGene_list2_2 <- append(RefGene_list2, normF_list2)
    
    save(RefGene_list2_2, RefGene_list2, file = "./data/RefGene2.Rda")
    #----------------------------#
# ____3.SET COMBO of interest ####
#----------------------------#
# if you didnt choose one of the option per parameter, all of them are included
# the inclusion of options will decrease the number of lists which will be taken into account for overlap analysis 

# choose RefGeneFinder Algorithm 
refGene_method      <- as.character("CV|geNorm|normF")            # "*CV|geNorm|normF*"       "CV"
# choose study; all_datasets also possible
table(metaDB_both_batches_1_filtered$Study)                       # "all_datasets" also available
study               <- as.character("Tiez_cancer")
# choose kind of disease
table(metaDB_both_batches_1_filtered$kind_of_disease)
kind_of_disease     <- as.character("heart_disease")
# choose kind of EV-preparation type
table(metaDB_both_batches_1_filtered$EV_preparation_type)
EV_isolation_Method <- as.character("*Precipitation")
# choose normalization method
names(count_list_normalized[[1]])
normalization_method<- as.character("TC_matrix|UQ_matrix|Med_matrix|DESeq_matrix|TMM_matrix|Quantile_matrix")        #"TC_matrix|UQ_matrix|Med_matrix|DESeq_matrix|TMM_matrix|Quantile_matrix"

#--------------------------------------------------------------------#
# ____4.make a shiny list / final list to table genes for overlap analysis ####
#--------------------------------------------------------------------#
# with RefGene_list()
list_shiny <- RefGene_list[grep(refGene_method,names(RefGene_list))]              # refGeneAlgorithm
list_shiny <- list_shiny[grep(study, names(list_shiny))]                          # study                 !!!!!! 
list_shiny <- list_shiny[grep(kind_of_disease, names(list_shiny))]                # kind_of_disease
list_shiny <- list_shiny[grep(EV_isolation_Method, names(list_shiny))]            # EV_preparation
list_shiny <- list_shiny[grep(normalization_method, names(list_shiny))]           # normalization method  !!!!!!
# with RefGene_list2()
list_shiny2 <- RefGene_list2[grep(refGene_method,names(RefGene_list2))]              # refGeneAlgorithm
list_shiny2 <- list_shiny2[grep(study, names(list_shiny2))]                          # study                 !!!!!! 
list_shiny2 <- list_shiny2[grep(kind_of_disease, names(list_shiny2))]                # kind_of_disease
list_shiny2 <- list_shiny2[grep(EV_isolation_Method, names(list_shiny2))]            # EV_preparation
list_shiny2 <- list_shiny2[grep(normalization_method, names(list_shiny2))]           # normalization method  !!!!!!

#-------------------------------------------------------------------------------------#
# ____5.check if shiny list have the right choosen parameters in the names of the lists ####
#-------------------------------------------------------------------------------------#
names(list_shiny)
names(list_shiny2)

#----------------------------------------------------#
# ____6. create 3x18 data frame with stability values  (refGeneAlg x normalization methods) ####
#--------------------------------------------------------# 
#TC genes extracted
overlap_normMeth_TC_lists2    <- RefGene_list2
#UQ genes extracted
overlap_normMeth_UQ_lists2     <- RefGene_list2
#Med genes extracted
overlap_normMeth_Med_lists2     <- RefGene_list2
#DESeq genes extracted
overlap_normMeth_DESeq_lists2     <- RefGene_list2
#TMM genes extracted
overlap_normMeth_TMM_lists2     <- RefGene_list2
#Q genes extracted
overlap_normMeth_Q_lists2     <- RefGene_list2

#join TC normalization
overlap_normMeth_TC_join <- merge(overlap_normMeth_TC_lists2$CV__all_datasets__TC_matrix,
                                  overlap_normMeth_TC_lists2$geNorm__all_datasets__TC_matrix, all = TRUE)
overlap_normMeth_TC_join <- merge(overlap_normMeth_TC_join,
                                  overlap_normMeth_TC_lists2$normF__all_datasets__TC_matrix, all = TRUE)
names(overlap_normMeth_TC_join) <- c("gene", "TC_CV", "TC_geNorm", "TC_normF")
# join UQ normalization
overlap_normMeth_join    <- merge(overlap_normMeth_TC_join,
                                  overlap_normMeth_UQ_lists2$CV__all_datasets__UQ_matrix, all = TRUE)
overlap_normMeth_join    <- merge(overlap_normMeth_join,
                                  overlap_normMeth_UQ_lists2$geNorm__all_datasets__UQ_matrix, all = TRUE)
overlap_normMeth_join    <- merge(overlap_normMeth_join,
                                  overlap_normMeth_UQ_lists2$normF__all_datasets__UQ_matrix, all = TRUE)
names(overlap_normMeth_join) <- c("gene", "TC_CV", "TC_geNorm", "TC_normF","UQ_CV", "UQ_geNorm", "UQ_normF")
# join med_normalization 
overlap_normMeth_join    <- merge(overlap_normMeth_join,
                                  overlap_normMeth_Med_lists2$CV__all_datasets__Med_matrix, all = TRUE)
overlap_normMeth_join    <- merge(overlap_normMeth_join,
                                  overlap_normMeth_Med_lists2$geNorm__all_datasets__Med_matrix, all = TRUE)
overlap_normMeth_join    <- merge(overlap_normMeth_join,
                                  overlap_normMeth_Med_lists2$normF__all_datasets__Med_matrix, all = TRUE)
names(overlap_normMeth_join) <- c("gene", "TC_CV", "TC_geNorm", "TC_normF","UQ_CV", "UQ_geNorm", "UQ_normF",";ed_CV", "Med_geNorm", "Med_normF")
# join DESeq normalization
overlap_normMeth_join    <- merge(overlap_normMeth_join,
                                  overlap_normMeth_DESeq_lists2$CV__all_datasets__DESeq_matrix, all = TRUE)
overlap_normMeth_join    <- merge(overlap_normMeth_join,
                                  overlap_normMeth_DESeq_lists2$geNorm__all_datasets__DESeq_matrix, all = TRUE)
overlap_normMeth_join    <- merge(overlap_normMeth_join,
                                  overlap_normMeth_DESeq_lists2$normF__all_datasets__DESeq_matrix, all = TRUE)
names(overlap_normMeth_join) <- c("gene", "TC_CV", "TC_geNorm", "TC_normF","UQ_CV", "UQ_geNorm", "UQ_normF","Med_CV", "Med_geNorm", "Med_normF","DESeq_CV", "DESeq_geNorm", "DESeq_normF")
# join TMM normalization
overlap_normMeth_join    <- merge(overlap_normMeth_join,
                                  overlap_normMeth_TMM_lists2$CV__all_datasets__TMM_matrix, all = TRUE)
overlap_normMeth_join    <- merge(overlap_normMeth_join,
                                  overlap_normMeth_TMM_lists2$geNorm__all_datasets__TMM_matrix, all = TRUE)
overlap_normMeth_join    <- merge(overlap_normMeth_join,
                                  overlap_normMeth_TMM_lists2$normF__all_datasets__TMM_matrix, all = TRUE)
names(overlap_normMeth_join) <- c("gene", "TC_CV", "TC_geNorm", "TC_normF","UQ_CV", "UQ_geNorm", "UQ_normF","Med_CV", "Med_geNorm", "Med_normF","DESeq_CV", "DESeq_geNorm", "DESeq_normF","TMM_CV", "TMM_geNorm", "TMM_normF")
# join Q normalization
overlap_normMeth_join    <- merge(overlap_normMeth_join,
                                  overlap_normMeth_Q_lists2$CV__all_datasets__Quantile_matrix, all = TRUE)
overlap_normMeth_join    <- merge(overlap_normMeth_join,
                                  overlap_normMeth_Q_lists2$geNorm__all_datasets__Quantile_matrix, all = TRUE)
overlap_normMeth_join    <- merge(overlap_normMeth_join,
                                  overlap_normMeth_Q_lists2$normF__all_datasets__Quantile_matrix, all = TRUE)
names(overlap_normMeth_join) <- c("gene", "TC_CV", "TC_geNorm", "TC_normF","UQ_CV", "UQ_geNorm", "UQ_normF","Med_CV", "Med_geNorm", "Med_normF","DESeq_CV", "DESeq_geNorm", "DESeq_normF","TMM_CV", "TMM_geNorm", "TMM_normF","Q_CV", "Q_geNorm", "Q_normF")
# count appearance of genes
counts_true <- as.data.frame(apply(overlap_normMeth_join[,2:19],1,function(x)(x != 0)))
counts_true <- t(counts_true)
overlap_normMeth_join[["total"]] <- as.vector(apply(counts_true, 1, table))

#----------------------------------------------------#
# ____7. create 1x6 data frame  (normalization methods irrespective of RefGeneAlg) ####
#--------------------------------------------------------# 
# FILTER for gene number in barplt
gene_number_plot_norm_methods <- 25
# generate 1 x 6 matrix
norm_meth_df <- data.frame("genes"  = overlap_normMeth_join$gene,
                           "TC_values" =paste(overlap_normMeth_join$TC_CV, overlap_normMeth_join$TC_geNorm, overlap_normMeth_join$TC_normF, sep = ","),
                           "UQ_values" =paste(overlap_normMeth_join$UQ_CV, overlap_normMeth_join$UQ_geNorm, overlap_normMeth_join$UQ_normF, sep = ","),
                           "Med_values" =paste(overlap_normMeth_join$Med_CV, overlap_normMeth_join$Med_geNorm, overlap_normMeth_join$Med_normF, sep = ","),
                           "DESeq_values" =paste(overlap_normMeth_join$DESeq_CV, overlap_normMeth_join$DESeq_geNorm, overlap_normMeth_join$DESeq_normF, sep = ","),
                           "TMM_values" =paste(overlap_normMeth_join$TMM_CV, overlap_normMeth_join$TMM_geNorm, overlap_normMeth_join$TMM_normF, sep = ","),
                           "Q_values" =paste(overlap_normMeth_join$Q_CV, overlap_normMeth_join$Q_geNorm, overlap_normMeth_join$Q_normF, sep = ",")
                           
)
# rename appearance of genes
norm_meth_df[] <- lapply(norm_meth_df, as.character)
for (i in 2:length(norm_meth_df)) {
  norm_meth_df[[i]] <-  mapvalues(norm_meth_df[[i]], from = "NA,NA,NA", to = "not appeared")
}

for (i in 2:length(norm_meth_df)) {
  norm_meth_df[[i]] <-  gsub(c(".*,.*"), c( "appeared"), norm_meth_df[[i]], perl = TRUE)
  
}
gene_number <- rowSums(norm_meth_df[-1] == "appeared")
norm_meth_df["total"] <- gene_number
# order according to frequency
norm_meth_df2 <- norm_meth_df[order(norm_meth_df$total, decreasing = TRUE),]
#plot results
p <- ggplot(data=norm_meth_df2[1:gene_number_plot_norm_methods,], aes(x=reorder(genes, -total), y=total))+
  ggtitle(paste("Frequency of best ranked genes \n across six normalization methods in", study))+
  geom_bar(stat="identity", fill="steelblue")+
  theme_minimal()+
  labs(y="Frequency", x="Genes")+
  theme(axis.text.x = element_text(angle = 90))+
  geom_text(aes(label = total), vjust=2.2, position = "identity", size = 3.5)

p

#----------------------------------------------------#
# ____8. create 1x3 data frame  (RefGeneAlg irrespective of normalization methods ) ####
#--------------------------------------------------------# 








#--------------------------------------------------------------------------------------------------------------------------------------------#
###----------------###
############ 12. Validation ##############
###----------------###
#--------------------------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------#
# ___1.analysis: overlap of RefGeneAlgorithms by choosing one dataset irrespective of normalization method ####
#                           -> choose dataset through shiny list + one refGeneAlg  + all 6 normalization methods
#                           -> output is a venn diagram and a list with combinations of overlap
#--------------------------------------------------------# 
#bestKeeper genes extracted
overlap_refGeneAlg_CV_lists     <- RefGene_list[grep("CV",names(RefGene_list))]
overlap_refGeneAlg_CV_lists     <- unlist(overlap_refGeneAlg_CV_lists)
overlap_refGeneAlg_CV_lists     <- as.character(overlap_refGeneAlg_CV_lists)
#geNorm genes extracted
overlap_refGeneAlg_geNorm_lists <- RefGene_list[grep("geNorm",names(RefGene_list))]
overlap_refGeneAlg_geNorm_lists <- unlist(overlap_refGeneAlg_geNorm_lists)
overlap_refGeneAlg_geNorm_lists <- as.character(overlap_refGeneAlg_geNorm_lists)
#normF genes extracted
overlap_refGeneAlg_normF_lists  <- RefGene_list[grep("normF",names(RefGene_list))]
overlap_refGeneAlg_normF_lists  <- unlist(overlap_refGeneAlg_normF_lists)
overlap_refGeneAlg_normF_lists  <- as.character(overlap_refGeneAlg_normF_lists)

# make a vectorlist wiht genes only once included
refGeneAlg_list <- list("bestKeeper" = as.vector(unique(overlap_refGeneAlg_CV_lists)),
                        "geNorm"= as.vector(unique(overlap_refGeneAlg_geNorm_lists)),
                        "normFinder" = as.vector(unique(overlap_refGeneAlg_normF_lists))
)
# veccompare is performed in server site of shiny to plot venn diagramm from vectorlists refGeneAlg_list


#----------------------------------------------------#
# ___2.analysis: overlap of normalization methods by choosing one dataset irrespective of RefGeneAlg ####
#                           1. choose dataset through shiny list + all refGeneAlg  + one of 6 normalization methods
#--------------------------------------------------------# 
#TC genes extracted
overlap_normMeth_TC_lists     <- RefGene_list[grep("TC_matrix",names(RefGene_list))]
overlap_normMeth_TC_lists     <- unlist(overlap_normMeth_TC_lists)
overlap_normMeth_TC_lists     <- as.character(overlap_normMeth_TC_lists)
#UQ genes extracted
overlap_normMeth_UQ_lists     <- RefGene_list[grep("UQ_matrix",names(RefGene_list))]
overlap_normMeth_UQ_lists     <- unlist(overlap_normMeth_UQ_lists)
overlap_normMeth_UQ_lists     <- as.character(overlap_normMeth_UQ_lists)
#Med genes extracted
overlap_normMeth_Med_lists     <- RefGene_list[grep("Med_matrix",names(RefGene_list))]
overlap_normMeth_Med_lists     <- unlist(overlap_normMeth_Med_lists)
overlap_normMeth_Med_lists     <- as.character(overlap_normMeth_Med_lists)
#DESeq genes extracted
overlap_normMeth_DESeq_lists     <- RefGene_list[grep("DESeq_matrix",names(RefGene_list))]
overlap_normMeth_DESeq_lists     <- unlist(overlap_normMeth_DESeq_lists)
overlap_normMeth_DESeq_lists     <- as.character(overlap_normMeth_DESeq_lists)
#TMM genes extracted
overlap_normMeth_TMM_lists     <- RefGene_list[grep("TMM_matrix",names(RefGene_list))]
overlap_normMeth_TMM_lists     <- unlist(overlap_normMeth_TMM_lists)
overlap_normMeth_TMM_lists     <- as.character(overlap_normMeth_TMM_lists)
#Q genes extracted
overlap_normMeth_Q_lists     <- RefGene_list[grep("Quantile_matrix",names(RefGene_list))]
overlap_normMeth_Q_lists     <- unlist(overlap_normMeth_Q_lists)
overlap_normMeth_Q_lists     <- as.character(overlap_normMeth_Q_lists)
# make a vectorlist wiht genes only once included
norm_meth_overlap_res <- list("TC" = as.vector(unique(overlap_normMeth_TC_lists)),
                        "Med"= as.vector(unique(overlap_normMeth_Med_lists)),
                        "FQ" = as.vector(unique(overlap_normMeth_Q_lists)),
                        "UQ" = as.vector(unique(overlap_normMeth_UQ_lists)),
                        "TMM" = as.vector(unique(overlap_normMeth_TMM_lists)),
                        "MoR" = as.vector(unique(overlap_normMeth_DESeq_lists))
)
# overlap through veccompare package
temp_over <- veccompare::compare.vectors(norm_meth_overlap_res)
# show which genes are overlapped in all normalization methods
temp_over[32][[1]]$overlap_of_elements
# show elements of overlap - than you can choose with the index, which genes overlap across different normalization methods
for (i in 1:length(temp_over)) {
  print(temp_over[[i]]$elements_involved)
}

overlap_normaalization <- list()
for (i in 1: length(temp_over)) {
  name <- paste("Overlap of" ,paste(temp_over[i][[1]]$elements_involved, collapse = " & "), collapse = " " )
  overlap_normaalization[[name]] <- as.factor(temp_over[i][[1]]$overlap_of_elements)
}

overlap_norm_meth <- as.data.frame(stri_list2matrix(overlap_normaalization, fill = NA_character_))
names(gg) <- names(overlap_normaalization)

#-------------------------------------------#
# ____Plot Frequency of genes in shiny list ####     
#-------------------------------------------#
# table result of list_shiny result
overlap <- data.frame(unlist(list_shiny))
overlap_result <- data.frame(table(overlap))
overlap_result <- overlap_result[order(overlap_result$Freq, decreasing = TRUE),]
names(overlap_result)[names(overlap_result) == "overlap"] <- paste("overlap_result_of",length(list_shiny), "combo lists and",length(overlap_result[,1]),"unique genes", sep = "__")

# plot best 10 genes with frequency
number_genes <- 15
df <- data.frame(supp = c(rep("total", number_genes), rep("appeared", number_genes)),
                 freq  = c( overlap_result[1:number_genes,2], rep(length(list_shiny), number_genes)-overlap_result[1:number_genes,2]),
                 gene = rep(overlap_result$`overlap_result_of__18__combo lists and__43__unique genes`[1:number_genes],2))

frequen <- c(overlap_result$Freq[1:number_genes], rep("",number_genes))
p <- ggplot(df, aes(x = reorder(gene, freq), y = freq))+
  geom_col(aes(fill = supp), width = 0.7)+
  ggtitle(paste("Frequency of best ranked genes in", length(list_shiny),"lists with", length(overlap_result[,1]),"unique genes" ))+
  labs(y="Number of lists", x="Best ranked genes")+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90))+
  geom_text(aes(label = frequen), vjust=3, position = "identity", size = 2.5) 
p

#----------------------------------------------------#
# ____Example evaluation: Precipitation & disease  ####     
#-----------------------------------------------------#
# precipitation & colorectal cancer
preci_colorec_canc_list <- list_shiny
preci_colorec_canc_gene <- unlist(preci_colorec_canc_list)
preci_colorec_canc_gene <- as.character(unique(preci_colorec_canc_gene))

# precipitation & none
preci_colorec_none_list <- list_shiny
preci_colorec_none_gene <- unlist(preci_colorec_none_list)
preci_colorec_none_gene <- as.character(unique(preci_colorec_none_gene))

# precipitation & heart disease
preci_colorec_heart_d_list <- list_shiny
preci_colorec_heart_d_gene <- unlist(preci_colorec_heart_d_list)
preci_colorec_heart_d_gene <- as.character(unique(preci_colorec_heart_d_gene))

# make a vectorlist 
refGeneAlg_list <- list("colorectal cancer" = as.vector(preci_colorec_canc_gene),
                        "none"= as.vector(preci_colorec_none_gene),
                        "heart disease" = as.vector(preci_colorec_heart_d_gene))
# make a data frame: bind all lists together and fill empty rows with "NA"
df      <- cbind.fill(preci_colorec_canc_gene,
                      preci_colorec_none_gene,
                      preci_colorec_heart_d_gene, fill = NA)
names(df) <- c("colorectal cancer", "none", "heart disease")
# veccompare
ra <- veccompare::compare.vectors(refGeneAlg_list)
# make venn diagram
h <- table(c(levels(df[,1]),levels(df[,2])), exclude = NA)   
length(which(h==2))
i <- table(c(levels(df[,2]),levels(df[,3])), exclude = NA)   
length(which(i==2))
j <- table(c(levels(df[,1]),levels(df[,3])), exclude = NA)   
length(which(j==2))
g <- table(c(levels(df[,1]),levels(df[,2]), levels(df[,3])), exclude = NA)
length(which(g == 3))
names_3 <- unique(names(which(g==3)))
names_2 <- unique(c(names(which(h==2)),names(which(i==2)), names(which(j==2))))

png(filename = "./precipitation_kind_of_disease.png", res = 300, width = 25, height = 13, units = "cm")
k <- draw.triple.venn(area1 = length(unique(df$`colorectal cancer`))-1,
                      area2 =  length(unique(df$none))-1,
                      area3 = length(unique(df$`heart disease`))-1,
                      n12 =  length(which(h == 2)),
                      n23 = length(which(i == 2)),
                      n13 = length(which(j == 2)),
                      n123 = length(which(g == 3)),
                      category = c("colorectal cancer", "none", "heart disease"),
                      lty = 3,
                      fill = c("skyblue", "pink1", "mediumorchid") ,
                      cex=2, cat.cex=1.2, cat.pos = 1,
                      cat.fontfamily = rep("serif", 3),
                      scaled = TRUE)
#grid.arrange(gTree(children=k), top="Title", bottom="Raw_Count matrix")
grid.arrange(gTree(children=k), top="Overlap of candidate miRNAs from precipitated samples \n and different kind of disease")
dev.off()
# overlap combinations data frame
precipitation_disease <- list("overlap all" = names(which(g==3)),
                              "overlap colorec_none"    = names(which(j==2)),
                              "overlap colorec_heart_disease"    = names(which(h==2)),
                              "overlap none_heart_disease"    = names(which(i==2)))

#----------------------------------------------------#
# ____Example evaluation: Precipitation & study  ####     
#-----------------------------------------------------#
# precipitation & Tiez_cancer
preci_tiez_list <- list_shiny
preci_tiez_gene <- unlist(preci_tiez_list)
preci_tiez_gene <- as.character(unique(preci_tiez_gene))

# precipitation & Pfaffl_St_HCH
preci_HCH_list <- list_shiny
preci_HCH_gene <- unlist(preci_HCH_list)
preci_HCH_gene <- as.character(unique(preci_HCH_gene))

# precipitation & Pfaffl_an_depr
preci_depr_list <- list_shiny
preci_depr_gene <- unlist(preci_depr_list)
preci_depr_gene <- as.character(unique(preci_depr_gene))

# make a vectorlist 
refGeneAlg_list <- list("tiez_cancer" = as.vector(preci_tiez_gene),
                        "pfaffl_hch"= as.vector(preci_HCH_gene),
                        "pfaffl_depr" = as.vector(preci_depr_gene))
# make a data frame: bind all lists together and fill empty rows with "NA"
df      <- cbind.fill(preci_tiez_gene,
                      preci_HCH_gene,
                      preci_depr_gene, fill = NA)
names(df) <- c("tiez_cancer", "pfaffl_hch", "pfaffl_depr")
# veccompare
ra <- veccompare::compare.vectors(refGeneAlg_list)
# make venn diagram
h <- table(c(levels(df[,1]),levels(df[,2])), exclude = NA)   
length(which(h==2))
i <- table(c(levels(df[,2]),levels(df[,3])), exclude = NA)   
length(which(i==2))
j <- table(c(levels(df[,1]),levels(df[,3])), exclude = NA)   
length(which(j==2))
g <- table(c(levels(df[,1]),levels(df[,2]), levels(df[,3])), exclude = NA)
length(which(g == 3))
names_3 <- unique(names(which(g==3)))
names_2 <- unique(c(names(which(h==2)),names(which(i==2)), names(which(j==2))))

png(filename = "./precipitation_study.png", res = 300, width = 25, height = 13, units = "cm")
k <- draw.triple.venn(area1 = length(unique(df$tiez_cancer))-1,
                      area2 =  length(unique(df$pfaffl_hch))-1,
                      area3 = length(unique(df$pfaffl_depr))-1,
                      n12 =  length(which(h == 2)),
                      n23 = length(which(i == 2)),
                      n13 = length(which(j == 2)),
                      n123 = length(which(g == 3)),
                      category = c("tiez_cacner", "pfaffl_hch", "pfaffl_depr"),
                      lty = 3,
                      fill = c("skyblue", "pink1", "mediumorchid") ,
                      cex=2, cat.cex=1.2, cat.pos = 1,
                      cat.fontfamily = rep("serif", 3),
                      scaled = TRUE)
#grid.arrange(gTree(children=k), top="Title", bottom="Raw_Count matrix")
grid.arrange(gTree(children=k), top="Overlap of candidate miRNAs from precipitated samples \n and different studies")
dev.off()
# overlap combinations data frame
precipitation_study <- list("overlap all" = names(which(g==3)),
                            "overlap tiez_hch"    = names(which(j==2)),
                            "overlap tiez_depr"    = names(which(h==2)),
                            "overlap hch_depr"    = names(which(i==2)))





