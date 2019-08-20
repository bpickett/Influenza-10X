library(Seurat)
library(dplyr)
library(Matrix)
library(homologene)
#code adapted from https://davetang.org/muse/2017/08/01/getting-started-seurat/

# change this to your working directory
setwd("~/Downloads/10X_files/")
data_dir <- "MDCK_0.06_outs"
#from_taxid_string <- "9031" #note that 9031 is chicken
from_taxid_string <- "9615" #note that 9615 is dog
pbmc.data <- Read10X(data.dir = data_dir)

# see https://stat.ethz.ch/R-manual/R-devel/library/Matrix/html/dgTMatrix-class.html
class(pbmc.data)
dim(pbmc.data)

# check out the first six genes and cells
pbmc.data[1:6, 1:6]
# summary of total expression per single cell
summary(colSums(pbmc.data))

# check how many genes have at least one transcript in each cell
at_least_one <- apply(pbmc.data, 2, function(x) sum(x>0))
hist(at_least_one, breaks = 100,
     main = "Distribution of detected genes",
     xlab = "Genes with at least one tag")

hist(colSums(pbmc.data),
     breaks = 100, main = "Expression sum per cell",
     xlab = "Sum expression")

# manually check the number of genes detected in three or more cells
# a lot of genes are not detected in 3 or more cells
tmp <- apply(pbmc.data, 1, function(x) sum(x>0))
table(tmp>=3)

# all cells have at least 200 detected genes
keep <- tmp>=3
tmp <- pbmc.data[keep,]
at_least_one <- apply(tmp, 2, function(x) sum(x>0))
summary(at_least_one)

pbmc <- CreateSeuratObject(raw.data = pbmc.data,
                           min.cells = 3,
                           min.genes = 200,
                           project = "10X_PBMC")

# see ?seurat for more information on the class
class(pbmc)
# same numbers as above 
pbmc

# mitochondria genes conveniently start with MT
mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
length(mito.genes)
if(length(mito.genes)>0){
  percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ]) / Matrix::colSums(pbmc@raw.data)
  # add some more meta data
  pbmc <- AddMetaData(object = pbmc,
                      metadata = percent.mito,
                      col.name = "percent.mito")
  # plot number of genes, UMIs, and % mitochondria
  VlnPlot(object = pbmc,
          features.plot = c("nGene", "nUMI", "percent.mito"),
          nCol = 3)
  GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito", pch.use = '.')
  }
# check out the meta data
head(pbmc@meta.data)

par(mfrow = c(1, 2))
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene", pch.use = '.')

# perform the filtering using FilterCells()
#pbmc <- FilterCells(object = pbmc,
#                    subset.names = c("nGene", "percent.mito"),
#                    low.thresholds = c(200, -Inf),
#                    high.thresholds = c(2500, 0.05))

# cells are filtered out; numbers consistent with above
pbmc
hist(colSums(pbmc@data),
     breaks = 100,
     main = "Total expression before normalization",
     xlab = "Sum of expression")

# currently there is only log-normalisation
pbmc <- NormalizeData(object = pbmc,
                      normalization.method = "LogNormalize",
                      scale.factor = 1e4)
hist(colSums(pbmc@data),
     breaks = 100,
     main = "Total expression after normalization",
     xlab = "Sum of expression")
pbmc_normalized <- as.data.frame(as.matrix(pbmc@data))

# the variable genes slot is empty before the analysis
#pbmc@var.genes
# refer to ?FindVariableGenes
pbmc <- FindVariableGenes(object = pbmc,
                          mean.function = ExpMean,
                          dispersion.function = LogVMR,
                          x.low.cutoff = 0.0125,
                          x.high.cutoff = 3,
                          y.cutoff = 0.5)
## vector of variable genes
#head(pbmc@var.genes)
## number of variable genes
#length(pbmc@var.genes)
## mean and variance of genes are stored pbmc@hvg.info
#head(pbmc@hvg.info)
## slot is empty before running ScaleData()
#pbmc@scale.data
## build linear model using nUMI and percent.mito
pbmc <- ScaleData(object = pbmc,
                  vars.to.regress = "nUMI")
#prepare to write scaled results to file
pbmc_scaled <- as.data.frame(as.matrix(pbmc@scale.data))
#pbmc@scale.data

pbmc <- RunPCA(object = pbmc,
               pc.genes = pbmc@var.genes,
               do.print = TRUE,
               pcs.print = 1:5,
               genes.print = 10)

#label data frames with genes that were included in PCA
pca_list <- pbmc@calc.params[["RunPCA"]][["pc.genes"]]
pca_list1 <- as.data.frame(pca_list)
pca_list1$is_pca <- rep(1,nrow(pca_list1)) # make new column
colnames(pca_list1) <- c("Genes","is_PCA")

pbmc_normalized$Genes <- rownames(pbmc_normalized)
pbmc_normalized = merge(pbmc_normalized, pca_list1, by.x="Genes", by.y="Genes", all = TRUE)
gene_list <- as.vector(pbmc_normalized$Genes)
temp <- homologene(gene_list, from_taxid_string, 9606)#convert gene symbols to human equivalent (where possible)
pbmc_normalized = merge(pbmc_normalized, temp, by.x="Genes", by.y=from_taxid_string, all = TRUE)

pbmc_scaled$Genes <- rownames(pbmc_scaled)
pbmc_scaled = merge(pbmc_scaled, pca_list1, by.x="Genes", by.y="Genes", all = TRUE)
gene_list <- as.vector(pbmc_scaled$Genes)
temp <- homologene(gene_list, from_taxid_string, 9606)#convert gene symbols to human equivalent (where possible)
pbmc_scaled = merge(pbmc_scaled, temp, by.x="Genes", by.y=from_taxid_string, all = TRUE)

write.csv(pbmc_normalized,file = paste0(data_dir,"_normalized.csv"))
write.csv(pbmc_scaled,file = paste0(data_dir,"_scaled.csv"))

#PrintPCAParams(pbmc)
# visualise top genes associated with principal components
#VizPCA(object = pbmc, pcs.use = 1:2)
#pbmc <- FindClusters(object = pbmc,
#                     reduction.type = "pca",
#                     dims.use = 1:10,
#                     resolution = 0.6,
#                     print.output = 1,
#                     save.SNN = TRUE)
#
#Output signficant genes to file
#output_file <- paste0("sig_genes","-",data_dir,".csv")
#df1 <- as.data.frame(pbmc@data)
#write.csv(pbmc@data,sep = ",",file = output_file)


# find all markers of cluster 1
#cluster1.markers <- FindMarkers(object = pbmc,
#                                test.use = "t",
#                                min.cells.feature = 3, 
#                                min.cells.group = 3,
#                                return.thresh = 0.25
#                                )
#cluster1.markers <- FindMarkers(object = pbmc,
#                                ident.1 = 1,
#                                min.pct = 0.25)
#head(cluster1.markers)
#
## find all markers distinguishing cluster 5 from clusters 0 and 3
#cluster5.markers <- FindMarkers(object = pbmc,
#                                ident.1 = 5,
#                                ident.2 = c(0,3),
#                                min.pct = 0.25)
#head(cluster5.markers)
#pbmc.markers <- FindAllMarkers(object = pbmc,
#                               only.pos = TRUE,
#                               min.pct = 0.25,
#                               thresh.use = 0.25)
#head(pbmc.markers)
#pbmc.markers %>% group_by(cluster) %>% top_n(2, avg_diff)
#
#table(pbmc@ident)
#my_t <- FindMarkers(object = pbmc,
#                    ident.1 = 1,
#                    thresh.use = 0.25,
#                    test.use = "t",
#                    only.pos = TRUE)
#dim(my_t)
#
#VlnPlot(object = pbmc, features.plot = c("TNPO3", "SMO"))
#FeaturePlot(object = pbmc,
#            features.plot = c("TNPO3", "SMO", "K123", "ENSGALG00000038834", "ADAMTS20", "PTPRB", "H2A-VII.6", "PLBD1", "HSP90B1"),
#            cols.use = c("grey", "blue"),
#            reduction.use = "tsne")
#head(pbmc.markers)
## stash cluster identities for later
#pbmc <- StashIdent(object = pbmc, save.name = paste0("ClusterNames_","data_dir"))
#
#                   #