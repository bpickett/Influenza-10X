library(cellrangerRkit)

#set pwd to parent folder for matrix file first **Note this has to be the parent directory for the "outs" folder.
##**Note** directory structure MUST be: "/sample_name/outs/filtered_gene_bc_matrices/outs", with the last "outs" containing the barcodes.tsv, genes.tsv, and matrix.mtx files.
cellranger_pipestance_path <- getwd()#"Desktop/filtered_gene_bc_matrices/Galgal_WF10_1"
gbm <- load_cellranger_matrix(cellranger_pipestance_path)
#analysis_results = load_cellranger_analysis_results(cellranger_pipestance_path)

###extract tSNE plot data for preliminary visualization
##tsne_proj=analysis_results$tsne
##visualize_umi_counts(gbm,tsne_proj[c("TSNE.1","TSNE.2")],limits=c(3,4),marker_size = 0.05)

#keep non-zero genes, normalize, and log10-transform the expression levels
use_genes=get_nonzero_genes(gbm)
gbm_bcnorm=normalize_barcode_sums_to_median(gbm[use_genes,])
gbm_log=log_gene_bc_matrix(gbm_bcnorm,base=10)

print(dim(gbm_log)) #gives dimensions of gene expression matrix (# non-zero genes=first value, # cells=second value)
##genes <- c("WF10PB2","WF10PB1","WF10PA","WF10HA","WF10NP","WF10NA","WF10M","WF10NS")
##visualize_gene_markers(gbm_log,genes,tsne_proj[c("TSNE.1","TSNE.2")],limits=c(0,1.5))

##first <- exprs(gbm_log)
##second <- fData(gbm_log) #data frame of genes
##third <- pData(gbm_log) #data frame of cell barcodes

###make a densely-populated csv file from the 10X matrix (MEX) output
subset_by_gene_id <- gbm_log[c("WF10PB2","WF10PB1","WF10PA","WF10HA","WF10NP","WF10NA","WF10M","WF10NS")]
print(dim(subset_by_gene_id))
first_subset <- as.matrix(exprs(subset_by_gene_id))
out <- t(first_subset)

tempString <- c("DF1_18","_","Normalized_Expression_Values_flu.csv")
outfile <- paste(tempString, collapse="")
write.csv(file=outfile,out, row.names=TRUE)
