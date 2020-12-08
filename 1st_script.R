### ---------------------------------------------------------------------------------------------------------###
### Automated finding of appropriate reference miRNAs in blood EVs readcount data sets (Alex Hildebrandt)
###----------------------------------------------------------------------------------------------------------###

# Description of the script ####

### tasks of the script:
### 
### 1. load necessary libraries 
### 2. load readcount table and metadata of interest - one readcount file for all samples and one metadata file for all samples 
### 3. define PARAMETERS: define all parameters which are necessary to run the script
#                  1: Filtering samples (% of mapped miRNA in ration to all readcounts)                    (-> 7% in first calculation)
#                  2: Filtering genes   (% of obligatory appearance of gene across all samples)            (-> 95% in first calculation)
#                  3: split combinations for readcount table
#                  4: treshhold CV      (% of genes under best ranked, which should be taken into account) (-> 5% above best stability value in first calculation)
#                  5: treshhold geNomr  (% of genes under best ranked, which should be taken into account) (-> 5% above best stability value in first calculation)
#                  6: treshhold NormF   (% of genes under best ranked, which should be taken into account) (-> 5% under best stability value in first calculation)
#                  7: treshhold sample  (number of samples, which have to appear in each list - lists with less samples will be removed -> necessary for NormF: 10 samples are recommended)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
### 4. FILTERING : samples of bad quality are removed - according to "Filtering samples" - PARAMETER
### 5. VALIDATION: check if readcount table and metadata have same order after filtering




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

#1. load libraries ####

library(data.table)
#----------------------------------------#
# 2. load readcount table and meta_Data   #####
#----------------------------------------#
# load files
load(file = "./miREV/meta_Data.Rda")                        
load(file = "./miREV/read_count.Rda")

#------------------------#
# 3. define Parameters   ####
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
# 4. filtering samples   #####
#------------------------#
# ___1.1 filtering samples: remove samples with less mapped miRNAs of total library size than 7%  
table(meta_Data$`%_mapped_miRNA_to_lib.size` >mapped_miRNA)
keep.samples <- meta_Data$`%_mapped_miRNA_to_lib.size` >mapped_miRNA

read_count_Data_1_filtered <- read_count_Data[,keep.samples]
meta_Data_1_filtered <- meta_Data[keep.samples,]

#-----------------------------------------------#
# 5. check if both sheets are in the same order ####
#-----------------------------------------------#
table(colnames(read_count_Data_1_filtered)==meta_Data_1_filtered$sample.name)    
