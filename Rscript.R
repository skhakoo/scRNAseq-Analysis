### PACKAGES ###
library(SoupX)
library(Seurat)
library(dplyr)
library(SeuratObject)
library(Matrix)
library(fields)
library(KernSmooth)
library(parallel)
library(ROCR)
library(EnhancedVolcano)
library(DoubletFinder)
library(ggplot2)
library(SingleR)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationHub)
library(tidyverse)
library(msigdbr)
library(clusterProfiler)
library(fgsea)



### READING DATA AND QC ###

#read in cellrangercount data specifying the pathway to the feature matrices
sample <- Read10X(data.dir = "~/Documents/Sophie/cellranger_counts/run_count70/outs/filtered_feature_bc_matrix")

#convert into a Seurat object for downstream analysis
sample <- CreateSeuratObject(counts = sample, project="chond_met") 

## quality control on seurat object ##
sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^MT-") #add a column with the percent of mitochondrial genes in each cell
sample <- subset(sample, subset = nFeature_RNA > 300) #remove cells with fewer than 300 genes 
sample <- subset(sample, subset = nFeature_RNA < 3000)#remove cells with more than 3000 genes 
sample <- subset(sample, subset = percent.mt < 10)#remove cells with more than 10% mitochondrial genes 

## removing doublets using doublet finder ##
data.filt <- NormalizeData(sample) #normalize data  
data.filt <- FindVariableFeatures(data.filt, selection.method = "vst", nfeatures = 2000) #find the top 2000 highly variable genes 
data.filt <- ScaleData(data.filt) #scale sata 
data.filt <- RunPCA(data.filt) #runPCA
data.filt <- RunUMAP(data.filt, dims = 1:10) #run UMAP with the first 10PCs
nExp <- round(ncol(data.filt)*0.04) #simulate the expected number of doublets (is 0.8% per 1000 cells)
data.filt <- doubletFinder_v3(data.filt, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10) #run doublet finder
DF.name = colnames(data.filt@meta.data)[grepl("DF.classification", colnames(data.filt@meta.data))]
filt.sample = data.filt[, data.filt@meta.data[, DF.name] =="Singlet"] #only keep cells which are defined as singlets
save(filt.sample, file = "~/Documents/Sophie/filtered_data/sample")


### NORMALIZATION ###
## SCTransform and normalization ## 

#calculate cell cycle scores to mitigate effects of cell cycle heterogeneity
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
sample <- CellCycleScoring(filt.sample, s.features = s.genes, g2m.features = g2m.genes)

#SCTransform replaces NormalizeData, ScaleData, and FindVariableFeatures
#correct for cell cycle scores and %MT genes by regressing them out 
sample <- SCTransform(sample, method = "glmGamPoi", vars.to.regress = c("percent.mt", 
                                                                        "S.Score", "G2M.Score"), verbose = F)

### CLUSTERING AND VISUALIZATION ###
## clustering ##
sample <- RunPCA(sample, verbose = F, dims = 1:20)
sample <- RunUMAP(sample, dims = 1:20)
sample <- FindNeighbors(sample, dims = 1:20, verbose = F)
sample <- FindClusters(sample, resolution = 0.5, verbose = F)
save(sample, file = "~/Documents/Sophie/analysed_data/GSM4601029")

## visualization ##

#open a jpeg file for the plot
library(ggplot2)
jpeg("sample76_cluster.jpeg", units="in", width=11, height=6, res=400) #open a jpeg file for the plot 
DimPlot(sample76, label=T) + labs(title='Cluster Graph GSM4601029') #generate the cluster plot for the sample
#saves the plot in the current working directory
dev.off()


jpeg("sample76_CAPN6_expression.jpeg", units="in", width=11, height=6, res=400) 
#generating a feature plot to visualise CAPN6 expression within the different clusters
FeaturePlot(sample76, features = "CAPN6", max.cutoff = 1) + labs(title='CAPN6 expression GSM4601029')
dev.off()



### DIFFERENTIAL GENE EXPRESSION ###

#Finding the markers that are differentially expressed in high CAPN6 expressing clusters (ident.1)
DEGs <- FindMarkers(sample, ident.1 = c(0, 2, 7), ident.2 = c(1, 3, 4, 5, 6, 7, 9), #idents.2 are the low CAPN6 expressing clusters 
                    min.pct = 0.25, logfc.threshold = 0.25) #min.pct specifies a minimum threshold of 25% of cells expressing the gene in either of the two groups of cells
#logfc.threshold removes genes with a minimum logfc of 0.25
DEGs %>%
  filter(p_val_adj < 0.05)
#saving the differentially expressed genes for further analysis 
save(DEGs, file = "~/Documents/Sophie/DEG/DEG_sample.RData")


## visualization of DEGs ##
df <- CAPN6_sample45_markers

#creating key-value pairs for 'high', 'mid', and 'low' epxression by fold-change
keyvals <- ifelse(df$avg_log2FC < -0.5, 'red',
                  ifelse(df$avg_log2FC >0.5, 'blue',
                         'grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals=='blue'] <- 'high'
names(keyvals)[keyvals=='grey'] <- 'mid'
names(keyvals)[keyvals=='red'] <- 'low'

library(EnhancedVolcano)
#opening a jpeg for the volcano plot 
jpeg("sample45_DEGs.jpeg", units="in", width=15, height=10, res=400)
EnhancedVolcano(df, 
                lab = rownames(df),
                x = 'avg_log2FC', 
                y = 'p_val', 
                selectLab = c('PTN','LUM','COL1A1','TIMP1','COL1A2',
                              'TNNC2','S100A13','FGFBP2','TPM2','CCN2',
                              'HLA-DRA','CD74','HLA-DRB1','C1QA','HLA-DPA1',
                              'C1QB','TMSB4X','CST3','HLA-DPB1','CCL3'),
                title = 'DGE in CAPN6 expressing clusters GSM4600998',
                drawConnectors = T, pointSize = 2, arrowheads = F, #drawConnectors adds connectors to fit more labels on the graph without arrowheads
                labSize = 6, colCustom = keyvals)  #colour graph based on custom key-value pairs
dev.off()     


### GO AND KEGG ENRICHMENT ###
CAPN6_sample70_markers$GEA <- "GSM4601023"
CAPN6_sample75_markers$GEA <- "GSM4601028"
CAPN6_sample76_markers$GEA <- "GSM4601029"
CAPN6_sample52_markers$GEA <- "GSM4601005"

rm(df)
chond_met<- rbind(CAPN6_sample45_markers, CAPN6_sample48_markers, CAPN6_sample50_markers, CAPN6_sample52_markers)
df <- chond_met %>% 
  filter(p_val_adj < 0.05,
         avg_log2FC >= 1) 


#turning the index genes into a column in the DEG matrix
list <- bitr(df$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") #generate a list of the DEGs and corresponding entrezID based on human annotation
genelist <- list$ENTREZID #extract the list of entrez ids to generate gene list for GO and KEGG

## GO enrichment ##

#return the enriched GO categories based on the list of entrez ids for all ontologies with a pvalue and qvalue cutoff of 0.05
GO_prim_all <- enrichGO(gene = genelist, OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, qvalueCutoff=0.05, ont="BP")
save(GO_prim_all, file = "~/Documents/Sophie/chondroblastic_primary/chond_prim_all2_GO.RData")

#create a dotplot for the enriched gene ontologies from the sample and save as a jpeg in the current directory
jpeg("chond_prim_all2_GO.jpeg", units="in", width=8.5, height=10, res=400)
dotplot(GO_prim_all, showCategory = 15) + labs(title = "Chondroblastic Primary GO")
dev.off()

## KEGG enrichment ##

#return the KEGG enrichment pathways based on list of entrez ids with a p and qvalue cutoff of 0.05
KEGG_met_all <- enrichKEGG(gene = genelist, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
save(KEGG_met_all, file = "~/Documents/Sophie/chondroblastic_primary/chond_met_all2_KEGG.RData")
#create a dotplot for the enriched KEGG ontologies from the sample and save as jpeg in the current directory
jpeg("chond_met_all2_KEGG.jpeg", units="in", width=8, height=9, res=400)
dotplot(KEGG_met_all, showCategory = 10) + labs(title = "Chondroblastic Metastasis KEGG")
dev.off()


### GENE SET ENRICHMENT ANALYSIS ###

## #HALLMARK ENRICHMENT ##

#Get gene set database
df <- DEGs

#turning the indexed genes into a column in the DEG matrix
df <- cbind(gene = rownames(df), df)
rownames(df) <- 1:nrow(df)

#obtain the 50 hallmark gene sets from the molecular signatures database 
H <- msigdbr(species = "Homo sapiens", category = "H")

## format the database ##

#creating a list of 50 hallmark gene sets and the corresponding gene symbols 
H.genes.ls <- H %>%
  select(gs_name, gene_symbol) %>% #extract the gene set name and gene symbol
  group_by(gs_name) %>% #group by gene set name (hallmark gene set)
  summarise(all.genes = list(unique(gene_symbol))) %>% #create a list of every unique gene symbols for each gene set name
  deframe() #

## Format for GSEA ##

df <- df %>%
  arrange(desc(avg_log2FC)) #rank the fold change in descending order
df.vec <- df$avg_log2FC #create a vector of the ranked logfc
names(df.vec) <- df$gene #create a vector with the gene symbol and logfc (named numeric vector)

## Set score type ##
min(df.vec)
max(df.vec)
scoreType <- "std"

## Run GSEA ##
gsea.H <- fgseaSimple(pathways = H.genes.ls, #list of 50 hallmark genes 
                      stats = df.vec, #list of gene symbols and logfc
                      scoreType = scoreType,
                      nperm = 100000) #the number of permeatations

## Plot gsea results ##
jpeg("sample_GSEA.jpeg", units="in", width=7, height=11, res=400)
gsea.H %>%
  filter(padj < 0.05) %>% #defining significant enriched gene sets as less than 0.05
  #Beautify descriptions by removing _and HALLMARK
  mutate(pathway = gsub("HALLMARK_", "", pathway),
         pathway = gsub("_", " ", pathway)) %>%
  
  ggplot(aes(x=reorder(pathway, NES), #Reorder gene sets by NES values
             y=NES)) +
  geom_col() +
  theme_classic() +
  #Flip x and y so long lables can be read
  coord_flip() +
  #fix labels
  labs(y="Normalized enrichment score (NES)",
       x="Gene set",
       title = "Hallmark GSEA GSM01029")
dev.off()


## CELL TYPE ENRICHMENT ANALYSIS ##
df <- DEGs

#obtain the cell type signature gene sets from the molecular signatures database 
C <- msigdbr(species = "Homo sapiens", category = "C8")

## Format GSEA database ##

#creating a list of the cell type signature gene sets and the corresponding gene symbols 
C.genes.ls <- C %>%
  select(gs_name, gene_symbol) %>% #selecting the gene set name and gene symbol
  group_by(gs_name) %>% #grouping these by gene set name 
  summarise(all.genes = list(unique(gene_symbol))) %>% #creating a list of unique gene symbols for each gene set
  deframe()

## format for gsea ##
df <- df %>%
  arrange(desc(avg_log2FC)) #rank the fold change in descending order
df.vec <- df$avg_log2FC #create a vector of the ranked logfc
names(df.vec) <- df$gene #create a vector with the gene symbol and logfc (named numeric vector)
df.vec

###set score type ###
min(df.vec)
max(df.vec)
scoreType <- "std"

### Run GSEA ###
gsea.C8 <- fgseaSimple(pathways = C.genes.ls, #list of cell type signature gene sets
                       stats = df.vec, #list of gene symbols and logfc
                       scoreType = scoreType,
                       nperm = 1000) #the number of permeatations

#subset vector to only include neutrophil siganture gene sets 
gsea.C8_subset <- gsea.C8[grep("NEUTROPHIL", gsea.C8$pathway), ]

## plot gsea results ##
jpeg("sample_neutrophil_GSEA.jpeg", units="in", width=7, height=4, res=400)
gsea.C8_subset %>%
  filter(padj < 0.05) %>% #defining significant enriched gene sets as less than 0.05
  ggplot(aes(x=reorder(pathway, NES), #Reorder gene sets by NES values
             y= NES)) +
  geom_col() +
  theme_classic() +
  #Flip x and y so long lables can be read
  coord_flip() +
  #fix labels
  labs(y="Normalized enrichment score (NES)",
       x="Gene set",
       title = "Neutrophil GSEA GSM4601029") +
  ylim(-3.5, 3.5)
dev.off()