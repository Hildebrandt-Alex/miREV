#----------------------------------------------------------------#
# 6. splitting counts and metaDB according to pheno-combinations #####
#----------------------------------------------------------------#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! use combn() to determine unique combinations 
#___1st level of kombination (with one variable) ####
# _______split disease #####
#first make a list of all possible comnination in one factor
number_levels <- length(levels(as.factor(meta_Data_1_filtered$kind_of_disease)))
temp_1 <- list()
for (i in 1:number_levels) {
  temp_1[[i]]  <- combn(levels(as.factor(meta_Data_1_filtered$kind_of_disease)),i)
}
#second split/merge count tables in list
split_list_disease <- list()
for (i in 1:length(temp_1)) {
  dd_sub <- temp_1[[i]]
  dd_sub <- as.data.frame(dd_sub)
  length_dd <- length(dd_sub)
  
  for (j in 1:length_dd) {
    name <- paste(temp_1[[i]][,j], collapse = ".")
    split_list_disease[[name]] <- meta_Data_1_filtered[meta_Data_1_filtered$kind_of_disease %in% as.character(temp_1[[i]][,j]),] 
  }
}
# _______split EV_isolation method #####-------------------------------------####
#first make a list of all possible comnination in one factor
number_levels <- length(levels(as.factor(meta_Data_1_filtered$EV_Isolation_Method_EV_track)))
temp_1 <- list()
for (i in 1:number_levels) {
  temp_1[[i]]  <- combn(levels(as.factor(meta_Data_1_filtered$EV_Isolation_Method_EV_track)),i)
}
#second split/merge count tables in list
split_list_ev_iso <- list()
for (i in 1:length(temp_1)) {
  dd_sub <- temp_1[[i]]
  dd_sub <- as.data.frame(dd_sub)
  length_dd <- length(dd_sub)
  
  for (j in 1:length_dd) {
    name <- paste(temp_1[[i]][,j], collapse = ".")
    split_list_ev_iso[[name]] <- meta_Data_1_filtered[meta_Data_1_filtered$EV_Isolation_Method_EV_track %in% as.character(temp_1[[i]][,j]),] 
  }
}
# _______split species #####----------------------------------------####
#first make a list of all possible comnination in one factor
number_levels <- length(levels(as.factor(meta_Data_1_filtered$species)))
temp_1 <- list()
for (i in 1:number_levels) {
  temp_1[[i]]  <- combn(levels(as.factor(meta_Data_1_filtered$species)),i)
}
#second split/merge count tables in list
split_list_species <- list()
for (i in 1:length(temp_1)) {
  dd_sub <- temp_1[[i]]
  dd_sub <- as.data.frame(dd_sub)
  length_dd <- length(dd_sub)
  
  for (j in 1:length_dd) {
    name <- paste(temp_1[[i]][,j], collapse = ".")
    split_list_species[[name]] <- meta_Data_1_filtered[meta_Data_1_filtered$species %in% as.character(temp_1[[i]][,j]),] 
  }
}
# _______split biofluids #####----------------------------------------####
#first make a list of all possible comnination in one factor
number_levels <- length(levels(as.factor(meta_Data_1_filtered$biofluids)))
temp_1 <- list()
for (i in 1:number_levels) {
  temp_1[[i]]  <- combn(levels(as.factor(meta_Data_1_filtered$biofluids)),i)
}
#second split/merge count tables in list
split_list_biofluids <- list()
for (i in 1:length(temp_1)) {
  dd_sub <- temp_1[[i]]
  dd_sub <- as.data.frame(dd_sub)
  length_dd <- length(dd_sub)
  
  for (j in 1:length_dd) {
    name <- paste(temp_1[[i]][,j], collapse = ".")
    split_list_biofluids[[name]] <- meta_Data_1_filtered[meta_Data_1_filtered$biofluids %in% as.character(temp_1[[i]][,j]),] 
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
