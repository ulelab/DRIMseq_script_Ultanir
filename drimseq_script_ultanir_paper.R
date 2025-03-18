#load packages
library(data.table)
library(DRIMSeq)
library(BiocParallel)
library(tictoc)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(RColorBrewer)
library(stageR)

#set working directory to location of count file directory, metadata file and poly(A) site atlas file
setwd("/Users/llywelyngriffith/Documents/AZ_postdoc/Simeon/Nobby_APA_analysis/common_atlas/mm10/")

#define PAS thresholds using PAS atlas
all_pas.dt = fread('merged_polya.filteredunique.annotated.bed', col.names = c("seqnames", "start", "end", "id", "score", "strand", "ensg", "hgnc", "region"))
pas.dt = all_pas.dt[region == "UTR3"]
pas.dt[, total_score := sum(score), by = ensg]
pas.dt[, percent_5 := 0.05 * total_score]
pas.dt[, percent_1 := 0.01 * total_score]
pas.dt[, multi := .N, by = ensg]

# Load and merge all count files
all.count.files = list.files("counts", full.names = TRUE, pattern = ".bed")
count.list = lapply(1:length(all.count.files), function(i) fread(all.count.files[[i]], select = c(4, 7, 8, 10), col.names = c("PAS", "ensg", "hgnc", gsub("_r1.polya_trimmed_window_merged_polya.bed", "", basename(all.count.files)[i]))))
count.dt = Reduce(function(x, y) merge(x = x, y = y, by = c("PAS", "ensg", "hgnc")), count.list)


# Remove intergenic PAS and those genes with only 1 PAS
count.dt = count.dt[ensg != "intergenic"]
count.dt = count.dt[PAS %in% pas.dt[multi > 1]$id]

#make columns called gene_id and feature_id for drimseq
counts.df = as.data.frame(count.dt[, `:=` (gene_id = paste0(hgnc, "_", ensg),
                                           feature_id = PAS)])

counts.df=counts.df[,-c(1:3)]

#make metadata dataframe
quantseq_metadata = read.table("quantseq_metadata.txt",
                               header = TRUE, as.is = TRUE)
quantseq_samples = data.frame(sample_id = quantseq_metadata$SampleName,
                              group = quantseq_metadata$Condition)

#### drimseq WT vs KO
quantseq_samples_pairwise = quantseq_samples[1:10,]
quantseq_samples_pairwise$group = factor(quantseq_samples_pairwise$group, levels = c("WT_brain","KO_brain"))
counts.df_pairwise=counts.df[,c(11:22)]

#create drimseq object
d_pairwise = dmDSdata(counts = counts.df_pairwise, samples = quantseq_samples_pairwise)
head(counts(d_pairwise), 3)
plotData(d_pairwise)

# Require gene expression in at least 75% of all samples, and PAS expression in 75% of samples for either DMSO or inhibition (whichever has fewest)
d_pairwise = dmFilter(d_pairwise, 
                      min_samps_gene_expr = floor(nrow(quantseq_samples_pairwise) * 0.75), 
                      min_samps_feature_expr = floor(min(table(quantseq_samples_pairwise$group)) * 0.75), 
                      min_gene_expr = 10, 
                      min_feature_expr = 5)

#drimseq design
design_full = model.matrix(~ group, data = DRIMSeq::samples(d_pairwise))
set.seed(42)
d_pairwise = dmPrecision(d_pairwise, design = design_full, verbose = 1)
d_pairwise = dmFit(d_pairwise, design = design_full, verbose = 1)
plotPrecision(d_pairwise)
head(coefficients(d_pairwise))

#run drimseq test
d_pairwise = dmTest(d_pairwise, coef = "groupKO_brain", verbose = 1)
View(results(d_pairwise))
View(results(d_pairwise,level="feature"))

#visualise results
res_pairwise_genes = results(d_pairwise)
res_pairwise_pA_sites = results(d_pairwise,level="feature")
res_pairwise_genes = res_pairwise_genes[order(res_pairwise_genes$adj_pvalue, decreasing = FALSE), ]
res_pairwise_pA_sites = res_pairwise_pA_sites[order(res_pairwise_pA_sites$adj_pvalue, decreasing = FALSE), ]
top_res_pairwise_genes_gene_id = res_pairwise_genes$gene_id[1]
plotProportions(d_pairwise, gene_id = top_res_pairwise_genes_gene_id , group_variable = "group")

#change NA values to 1s for two stage test
no.na = function(x) ifelse(is.na(x), 1, x)
res_pairwise_genes$pvalue = no.na(res_pairwise_genes$pvalue)
res_pairwise_pA_sites$pvalue = no.na(res_pairwise_pA_sites$pvalue)

##two stage test
#assign gene-level pvalues to the screening stage
pScreen_pairwise = res_pairwise_genes$pvalue
names(pScreen_pairwise) = res_pairwise_genes$gene_id

#assign transcript level p values to the confirmation stage
pConfirmation_pairwise = matrix(res_pairwise_pA_sites$pvalue, ncol = 1)
rownames(pConfirmation_pairwise) = res_pairwise_pA_sites$feature_id

#create the gene-transcript mapping
tx2gene_pairwise = res_pairwise_pA_sites[,c("feature_id","gene_id")]

#create the stageRTx object and perform the stage-wise analysis
stageR0bj_pairwise = stageRTx(pScreen = pScreen_pairwise, 
                              pConfirmation = pConfirmation_pairwise, 
                              pScreenAdjusted = FALSE, 
                              tx2gene = tx2gene_pairwise)

stageR0bj_pairwise  = stageWiseAdjustment(object = stageR0bj_pairwise, method = "dtu", alpha = 0.05)

getSignificantGenes(stageR0bj_pairwise)
getSignificantTx(stageR0bj_pairwise)

two_stage_padj_pairwise =getAdjustedPValues(stageR0bj_pairwise, order=FALSE, onlySignificantGenes = FALSE)

#change column names
colnames(two_stage_padj_pairwise)[1]='gene_id'
colnames(two_stage_padj_pairwise)[2]='feature_id'
colnames(two_stage_padj_pairwise)[3]='twostep_gene_padj'
colnames(two_stage_padj_pairwise)[4]='twostep_transcript_padj'

#add column of differences in usage between each condition
proportions_pairwise=proportions(d_pairwise)
proportions_pairwise$change_in_usage = (rowMeans(proportions_pairwise[10]) - rowMeans(proportions_pairwise[5]))
res_pairwise_pA_sites_with_change_in_usage = left_join(res_pairwise_pA_sites,proportions_pairwise, by = c("gene_id","feature_id"))

#combine two-test table with change in usage table
two_step_WT_vs_KO_brain_pA_sites_with_change_in_usage=inner_join(two_stage_padj_pairwise,res_pairwise_pA_sites_with_change_in_usage,by=c("gene_id","feature_id"))
sig_two_step_WT_vs_KO_brain_pA_sites_with_change_in_usage=two_step_WT_vs_KO_brain_pA_sites_with_change_in_usage %>% filter(twostep_transcript_padj <= 0.05 & abs(change_in_usage)>=0.1)

#process PAS atlas table
colnames(all_pas.dt)[1:4] <- c("chr", "start", "end", "feature_id")
colnames(all_pas.dt)[6] <- "strand"
polyA_coords <- all_pas.dt[, c(1:4, 6)]

#combine poly A table with coordinate table
two_step_WT_vs_KO_brain_pA_sites_with_change_in_usage_with_coords = inner_join(two_step_WT_vs_KO_brain_pA_sites_with_change_in_usage,polyA_coords, by = c("feature_id"))

#write to file
write_csv(two_step_WT_vs_KO_brain_pA_sites_with_change_in_usage_with_coords,"WT_vs_KO_brain_pA_sites_drimseq_results.csv")

#### drimseq WT vs HOM
quantseq_samples_pairwise = quantseq_samples[11:20,]
quantseq_samples_pairwise$group = factor(quantseq_samples_pairwise$group, levels = c("ELAVL1_WT","ELAVL1_HOM"))
counts.df_pairwise=counts.df[,c(1:10,21:22)]

#create drimseq object
d_pairwise = dmDSdata(counts = counts.df_pairwise, samples = quantseq_samples_pairwise)
head(counts(d_pairwise), 3)
plotData(d_pairwise)

# Require gene expression in at least 75% of all samples, and PAS expression in 75% of samples for either DMSO or inhibition (whichever has fewest)
d_pairwise = dmFilter(d_pairwise, 
                      min_samps_gene_expr = floor(nrow(quantseq_samples_pairwise) * 0.75), 
                      min_samps_feature_expr = floor(min(table(quantseq_samples_pairwise$group)) * 0.75), 
                      min_gene_expr = 10, 
                      min_feature_expr = 5)

#drimseq design
design_full = model.matrix(~ group, data = DRIMSeq::samples(d_pairwise))
set.seed(42)
d_pairwise = dmPrecision(d_pairwise, design = design_full, verbose = 1)
d_pairwise = dmFit(d_pairwise, design = design_full, verbose = 1)
plotPrecision(d_pairwise)
head(coefficients(d_pairwise))

#run drimseq test
d_pairwise = dmTest(d_pairwise, coef = "groupELAVL1_HOM", verbose = 1)

#visualise results
res_pairwise_genes = results(d_pairwise)
res_pairwise_pA_sites = results(d_pairwise,level="feature")
res_pairwise_genes = res_pairwise_genes[order(res_pairwise_genes$adj_pvalue, decreasing = FALSE), ]
res_pairwise_pA_sites = res_pairwise_pA_sites[order(res_pairwise_pA_sites$adj_pvalue, decreasing = FALSE), ]
top_res_pairwise_genes_gene_id = res_pairwise_genes$gene_id[1]
plotProportions(d_pairwise, gene_id = top_res_pairwise_genes_gene_id , group_variable = "group")

#change NA values to 1s for two stage test
no.na = function(x) ifelse(is.na(x), 1, x)
res_pairwise_genes$pvalue = no.na(res_pairwise_genes$pvalue)
res_pairwise_pA_sites$pvalue = no.na(res_pairwise_pA_sites$pvalue)

##two stage test
#assign gene-level pvalues to the screening stage
pScreen_pairwise = res_pairwise_genes$pvalue
names(pScreen_pairwise) = res_pairwise_genes$gene_id

#assign transcript level pvalues to the confirmation stage
pConfirmation_pairwise = matrix(res_pairwise_pA_sites$pvalue, ncol = 1)
rownames(pConfirmation_pairwise) = res_pairwise_pA_sites$feature_id

#create the gene-transcript mapping
tx2gene_pairwise = res_pairwise_pA_sites[,c("feature_id","gene_id")]

#create the stageRTx object and perform the stage-wise analysis
stageR0bj_pairwise = stageRTx(pScreen = pScreen_pairwise, 
                              pConfirmation = pConfirmation_pairwise, 
                              pScreenAdjusted = FALSE, 
                              tx2gene = tx2gene_pairwise)

stageR0bj_pairwise  = stageWiseAdjustment(object = stageR0bj_pairwise, method = "dtu", alpha = 0.05)

getSignificantGenes(stageR0bj_pairwise)
getSignificantTx(stageR0bj_pairwise)

two_stage_padj_pairwise =getAdjustedPValues(stageR0bj_pairwise, order=FALSE, onlySignificantGenes = FALSE)

colnames(two_stage_padj_pairwise)[1]='gene_id'
colnames(two_stage_padj_pairwise)[2]='feature_id'
colnames(two_stage_padj_pairwise)[3]='twostep_gene_padj'
colnames(two_stage_padj_pairwise)[4]='twostep_transcript_padj'

#add column of differences in usage between each condition
proportions_pairwise=proportions(d_pairwise)
proportions_pairwise$change_in_usage = (rowMeans(proportions_pairwise[10]) - rowMeans(proportions_pairwise[5]))
res_pairwise_pA_sites_with_change_in_usage = left_join(res_pairwise_pA_sites,proportions_pairwise, by=c("gene_id","feature_id"))

#combine two-test table with change in usage table table
two_step_WT_vs_HOM_pA_sites_with_change_in_usage=inner_join(two_stage_padj_pairwise,res_pairwise_pA_sites_with_change_in_usage,by=c("gene_id","feature_id"))
sig_two_step_WT_vs_HOM_pA_sites_with_change_in_usage=two_step_WT_vs_HOM_pA_sites_with_change_in_usage %>% filter(twostep_transcript_padj <= 0.05 & abs(change_in_usage)>=0.1)

#combine poly A table with coordinate table
two_step_WT_vs_HOM_pA_sites_with_change_in_usage_with_coords = inner_join(two_step_WT_vs_HOM_pA_sites_with_change_in_usage,polyA_coords, by = "feature_id")

#write to file
write_csv(two_step_WT_vs_HOM_pA_sites_with_change_in_usage_with_coords,"WT_vs_HOM_brain_pA_sites_drimseq_results.csv")
