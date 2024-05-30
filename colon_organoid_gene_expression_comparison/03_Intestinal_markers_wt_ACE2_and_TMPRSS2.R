library(ComplexHeatmap) # Heatmap
library(circlize)

##### Create a heatmap of log2FC of ACE2, TMPRSS2, and intestinal marker genes (Supplementary Figure 4A)
colon_gene_df <- read.delim("~/240213_Colon_expression/240213_colon_gene_expression.txt", header=TRUE)
colnames(colon_gene_df) <- c("Group","ACE2","CYP2C19","SLC51A","SLC51B","CA2","CDX2","TMPRSS2","VIL1","ABCC2")
colon_gene_df$Group <- c("Colon organoid","Caco-2","Adult small intestine","Fetal small intestine")
rownames(colon_gene_df) <- colon_gene_df$Group
colon_gene_df[,2:10] <- log2(colon_gene_df[,2:10])
colon_gene_df[,2:10] <- apply(colon_gene_df[,2:10], 2, scale)
colon_gene_df$Group <- NULL
colon_gene_df <- colon_gene_df[,sort(c("ACE2","CYP2C19","SLC51A","SLC51B","CA2","CDX2","TMPRSS2","VIL1","ABCC2"))]

# pdf("colon_gene_heatmap.pdf", width=5.7, height=6)
colon_gene_heatmap <- Heatmap(t(colon_gene_df), show_row_names=TRUE, row_names_side="left", row_title_rot=90, row_names_gp=gpar(fontsize=14, fontface="italic"), row_gap=unit(2.5, "mm"), cluster_rows=FALSE, row_title_gp=gpar(fontface="bold", fontsize=18), clustering_method_rows="ward.D2", 
                              show_column_names=TRUE, column_gap=unit(2.5, "mm"), column_title_gp=gpar(fontface="bold", fontsize=18), cluster_columns=TRUE, column_names_gp=gpar(fontsize=14), clustering_method_columns="ward.D2", 
                              col=colorRamp2(c(-2, 0, 2), c("blue","white","red")), na_col="lightgrey", heatmap_legend_param=list(title="z score of log2 fold change", title_gp=gpar(fontface="bold", fontsize=16), labels_gp=gpar(fontsize=16)))
colon_gene_heatmap
# dev.off()


