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

