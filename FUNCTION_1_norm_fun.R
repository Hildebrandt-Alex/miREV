#----------------------------#
### normalization function ####
#---------------------------#


norm_function          <- function(count_matrix_raw,
                                   metadata
                                   )
{
  #---------------------------------------------------#
  # 1. create empty matrix sheet for each Normalization method ############
  #-------------------------------------------------------#
  RawCountMatrix          <-  matrix(NA, nrow(count_matrix_raw), ncol(count_matrix_raw))
  TotalCountMatrix        <-  matrix(NA, nrow(count_matrix_raw), ncol(count_matrix_raw))
  UpperQuantileMatrix     <-  matrix(NA, nrow(count_matrix_raw), ncol(count_matrix_raw))
  MedianMatrix            <-  matrix(NA, nrow(count_matrix_raw), ncol(count_matrix_raw))
  DESeqMatrix             <-  matrix(NA, nrow(count_matrix_raw), ncol(count_matrix_raw))
  TrimmedMeanMatrix       <-  matrix(NA, nrow(count_matrix_raw), ncol(count_matrix_raw))
  QuantileMatrix          <-  matrix(NA, nrow(count_matrix_raw), ncol(count_matrix_raw))
  
  rownames_matrix <- rownames(count_matrix_raw)
  colnames_matrix <- colnames(count_matrix_raw)
  #----------------------------------------------------------------------------------------------------#
  raw_df <- count_matrix_raw
  
  #---------------------------------------#
  # 2. apply normalization on the input matrix #########
  #------------------------------------------#
  ## quantile normalization     preprocessCore
  QuantileMatrix <- preprocessCore::normalize.quantiles(as.matrix(count_matrix_raw))
  rownames(QuantileMatrix) <- rownames(count_matrix_raw)
  colnames(QuantileMatrix) <- colnames(count_matrix_raw)
  #------------------------------------------------------------------------------------------------------------------------#
  ## TMM normalization  edgeR
  TrimmedMeanMatrix  <- edgeR::calcNormFactors(as.matrix(count_matrix_raw), method="TMM")
  TrimmedMeanMatrix  <- colSums(count_matrix_raw)*TrimmedMeanMatrix / mean(colSums(count_matrix_raw)*TrimmedMeanMatrix)
  TrimmedMeanMatrix  <- scale(count_matrix_raw, center = FALSE, scale = TrimmedMeanMatrix)
  #----------------------------------------------------------------------------------------------------------------------#
  ##Total count normalization 
  TotalCountMatrix  <- colSums(count_matrix_raw)/mean(colSums(count_matrix_raw))
  TotalCountMatrix  <- scale(count_matrix_raw, center = FALSE, scale = TotalCountMatrix)
  #-----------------------------------------------------------------------------------------------------------------------#
  ## Median normalizaiton   
  MedianMatrix  <- apply(count_matrix_raw, 2, median)
  MedianMatrix  <- MedianMatrix/mean(MedianMatrix)
  MedianMatrix  <- scale(count_matrix_raw, center = FALSE, scale = MedianMatrix)
  #-----------------------------------------------------------------------------------------------------------------------#
  ## DESeq2 normalization - Median ratio
  DESeqMatrix <-  DESeq2::DESeqDataSetFromMatrix(countData = count_matrix_raw, colData = metadata, design = ~ 1)
  dds <- DESeq2::estimateSizeFactors(DESeqMatrix)
  DESeqMatrix <- DESeq2::counts(dds, normalized=TRUE)
  #-----------------------------------------------------------------------------------------------------------------------#
  ## UQ normalizaiton    edgeR
  UpperQuantileMatrix  <- apply(count_matrix_raw,2, quantile, 0.75)
  UpperQuantileMatrix  <- UpperQuantileMatrix/mean(UpperQuantileMatrix)
  UpperQuantileMatrix  <- scale(count_matrix_raw, center = FALSE, scale = UpperQuantileMatrix)
  
  #--------------------------------------------------------------------#
  # 3. output of the function is list with different Normalized Matrix's ####
  #-----------------------------------------------------------------------#
  list(Raw_matrix       = raw_df,
       TC_matrix        = TotalCountMatrix,
       UQ_matrix        = UpperQuantileMatrix,
       Med_matrix       = MedianMatrix,
       DESeq_matrix     = DESeqMatrix,
       TMM_matrix       = TrimmedMeanMatrix,
       Quantile_matrix  = QuantileMatrix
  )
}