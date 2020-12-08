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
