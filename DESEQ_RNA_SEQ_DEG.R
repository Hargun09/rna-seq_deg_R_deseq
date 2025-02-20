

## Install all packages
     
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

install.packages("ggplot2")
install.packages("gplots")
install.packages("readr")
install.packages("pheatmap")



if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EnhancedVolcano")


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("apeglm")
       
library(readr)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)


ed <- read_csv("data_1.csv",col_types = cols(C1 = col_integer(), 
                                             C2 = col_integer(), 
                                             C3 = col_integer(), 
                                             C4 = col_integer(), 
                                             C5 = col_integer(),
                                             C6 = col_integer(), 
                                             C7 = col_integer(), 
                                             C8 = col_integer(), 
                                             C9 = col_integer(), 
                                             C10 = col_integer(),
                                             C11 = col_integer(),
                                             C12 = col_integer(),
                                             C13 = col_integer(),
                                             C14 = col_integer(),
                                             C15 = col_integer(),
                                             C16 = col_integer(),
                                             C17 = col_integer(),
                                             C18 = col_integer(),
                                             D1 = col_integer(),
                                             D2 = col_integer(),
                                             D3 = col_integer(),
                                             D4 = col_integer(),
                                             D5 = col_integer(),
                                             D6 = col_integer(),
                                             D7 = col_integer(),
                                             D8 = col_integer(),
                                             D9 = col_integer(),
                                             D10 = col_integer(),
                                             D11 = col_integer(),
                                             D12 = col_integer(),
                                             D13 = col_integer(),
                                             D14 = col_integer(),
                                             D15 = col_integer(),
                                             D16 = col_integer(),
                                             D17 = col_integer(),
                                             D18 = col_integer()))





counts <- data.frame(ed)
row.names(counts) <- counts$Geneid
counts
final<- counts [, -1]
head(final)
write.csv(final ,'data_primary_indexed_final.csv')


  
condition <- factor(c("C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","D","D","D","D","D","D","D","D","D","D","D","D","D","D","D","D","D","D"))
       
coldata <- data.frame(row.names = colnames(final), condition)
        
any(is.na(final))

Counts <- final[complete.cases(final), ]`
   


dds <- DESeqDataSetFromMatrix(Counts,coldata,~condition)

dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "C", "D"))

head(res)

summary(res)
         
sigs <- na.omit(res)
sigs <- sigs[sigs$padj < 0.05,]
 
vsd <- vst(dds)
plotPCA(vsd)
plotPCA(vsd, ntop=500, "condition") #An alternate type of PCA plot for top 500 most variable genes.
   
down_genes <- sigs[sigs$log2FoldChange < 0,]
up_genes <- sigs[sigs$log2FoldChange > 0,]
       
write.csv(sigs, "sigs.csv")
write.csv(down_genes, "down_genes.csv")
write.csv(up_genes, "up_genes.csv")
      

plotDispEsts(dds)
    

library(pheatmap)

de_genes <- rownames(sigs)

counts_de <- Counts[de_genes,]

pheatmap(counts_de, 
         scale = "row", 
         clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation", 
         clustering_method = "complete",
         show_rownames = FALSE,
         show_colnames = FALSE)

EnhancedVolcano(res,
                lab = as.character(rownames(res)),  
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-5,5),  
                ylim = c(0, 7))

plotMA(res)


pca <- prcomp(t(assay(vsd))) # Performing PCA on data transformed previously by vst()

pca_data <- pca$x[, 1:3] #Extracting PCA into three parts

conditions <- coldata$condition
colors <- ifelse(conditions == "C", "pink", "darkgreen")

open3d()

plot3d(pca_data, size=2, type="s", box=F,
       xlab="PC1", ylab="PC2", zlab="PC3",
       col=colors)


bbox3d(color=c("#483D8B", "black"), emission="#483D8B",
       specular="#483D8B", shininess=1, alpha=0.8, nticks=3) #defining colours of the bounding box

legend3d("topright", legend = levels(conditions), pch = 16, col = c("pink", "darkgreen"))




