
#----------------------------------------#
# 1. load readcount table and metadata   #####
#----------------------------------------#
# load files
load(file = "./miREV/metadata.Rda")                        
load(file = "./miREV/read_count.Rda")

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
