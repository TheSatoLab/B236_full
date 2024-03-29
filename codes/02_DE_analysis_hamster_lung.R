library(zFPKM)           # zFPKM analysis
library(ComplexHeatmap)  # Heatmap
library(DESeq2)          # DESeq2
library(PCAtools)        # PCA
library(clusterProfiler) # Enrichment analysis
library(rstatix)         # Hypothesis testing

# Miscellaneous
library(readr)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(circlize)
library(data.table)
library(RColorBrewer)
library(cowplot)

### Change working directory
dir <- "~/Desktop/BANAL_project/230404_Hamster_lung/"
setwd(dir)

######################################## Analyze read counts with total viral reads ########################################
########## Read count table and record grouping information
count.df <- read.table(file.path(dir, "B236_hamster_counts.txt"), header=TRUE)
rownames(count.df) <- count.df$Geneid

samples <- read.delim(file.path(dir, "B236_hamster_samples.tsv"), header=TRUE)

samples <- samples %>% mutate(Virus = ifelse(Virus=="Control", "Saline",
                                      ifelse(Virus=="WK-521", "SC2",
                                      ifelse(Virus=="WK-521 dFCS", "SC2dFCS",
                                      ifelse(Virus=="BANAL-236", "B236", NA)))))
samples <- samples %>% mutate(Grouping = paste(Virus, Tissue, Day, sep="_"))

sample_dict <- unique(samples[,c("Sample_bam","Sample")])

passed_samples <- samples %>% filter(RIN > 7)
passed_samples.v <- (samples %>% filter(RIN > 7))$Sample

colnames(count.df) <- c(colnames(count.df)[1:6], sample_dict$Sample[match(colnames(count.df[7:ncol(count.df)]), sample_dict$Sample_bam)])

count.df_deseq <- count.df
count.df_deseq[1:6] <- NULL
count.df_deseq <- count.df_deseq[,passed_samples.v]

##### Since the genome length of B236 (29,844 bp) and that of SC2 (29,903 bp) are similar, it's unnecessary to normalize counts with virus genome length. #####

total_count_per_each_exp <- colSums(count.df_deseq)
total_viral_count.df <- count.df_deseq[(nrow(count.df_deseq)-1):nrow(count.df_deseq),]
norm_viral_count.mx <- sweep(total_viral_count.df*1000000,2,total_count_per_each_exp,FUN="/")
norm_viral_count.df <- as.data.frame.table(as.matrix(norm_viral_count.mx))
colnames(norm_viral_count.df) <- c("Genome","Grouping","norm_count")

norm_viral_count.df$Genome <- str_replace(norm_viral_count.df$Genome, "WK-521", "SARS-CoV-2")
norm_viral_count.df$Genome <- str_replace(norm_viral_count.df$Genome, "MZ937003.2", "B236")
norm_viral_count.df$Virus <- str_split(str_split(norm_viral_count.df$Grouping, "dpi_", simplify=T)[,2], "_", simplify=T)[,1]
norm_viral_count.df$Virus <- str_replace(norm_viral_count.df$Virus, "WK", "SC2")
norm_viral_count.df$Virus <- str_replace(norm_viral_count.df$Virus, "delFCS", "SC2dFCS")
norm_viral_count.df$Virus <- str_replace(norm_viral_count.df$Virus, "BANAL", "B236")
norm_viral_count.df$Tissue <- str_split_i(norm_viral_count.df$Grouping, "_", -1)
norm_viral_count.df$Tissue <- str_replace(norm_viral_count.df$Tissue, "h", "Hilum")
norm_viral_count.df$Tissue <- str_replace(norm_viral_count.df$Tissue, "p", "Periphery")
norm_viral_count.df$Day <- str_split_i(norm_viral_count.df$Grouping, "_", -4)
norm_viral_count.df$log_norm_count <- log10(norm_viral_count.df$norm_count)

norm_viral_count.df <- norm_viral_count.df %>% filter(Virus!="Saline") %>% filter((Genome=="SARS-CoV-2" & Virus=="SC2") | (Genome=="SARS-CoV-2" & Virus=="SC2dFCS") | (Genome=="B236" & Virus=="B236"))
norm_viral_count.df$Virus <- factor(norm_viral_count.df$Virus, levels=c("SC2","SC2dFCS","B236"))
norm_viral_count.df$Genome <- factor(norm_viral_count.df$Genome, levels=c("SARS-CoV-2","B236"))
norm_viral_count.df$Tissue <- factor(norm_viral_count.df$Tissue, levels=c("Hilum","Periphery"))
norm_viral_count.df$Day <- factor(norm_viral_count.df$Day, levels=c("2dpi","5dpi"))
norm_viral_count.df$Grouping <- apply(norm_viral_count.df[,c("Virus","Tissue","Day")], 1, paste, collapse="_")



########## Calculate the normalized counts of SC2 and B236 in hamster lung hila and peripheries (Figure 3F)
norm_viral_count.df %>% group_by(Tissue, Virus, Day) %>% summarise(median=median(log_norm_count))
norm_viral_count.df %>% filter(!Grouping == "SC2dFCS_Periphery_5dpi") %>% group_by(Tissue, Virus, Day) %>% shapiro_test(log_norm_count) # Sample size of SC2dFCS_Periphery_5dpi group is less than 3, cannot tested by Shapiro test
norm_viral_count.df %>% group_by(Tissue, Day) %>% levene_test(log_norm_count~Virus)

### Use Kruskal-Wallis since assumption of normality was violated
res.kruskal <- norm_viral_count.df %>% group_by(Tissue, Day) %>% kruskal_test(log_norm_count~Virus)
res.kruskal

norm_viral_count.df_stats <- norm_viral_count.df %>% group_by(Tissue, Day) %>% dunn_test(log_norm_count~Virus, p.adjust.method="fdr") # FDR (alias of BH)
norm_viral_count.df_stats <- norm_viral_count.df_stats %>% add_xy_position(x="Day", fun="max", dodge=0.8)
norm_viral_count.df_stats

# pdf("viral_count_hamster_lung.pdf", width=9, height=6)
norm_viral_count_plot.lung <- ggplot(norm_viral_count.df, aes(x=Day, y=log_norm_count, fill=Virus), alpha=0.5)
norm_viral_count_plot.lung <- norm_viral_count_plot.lung + geom_jitter(alpha=0.9, shape=16, size=3, position=position_jitterdodge(0.5)) 
norm_viral_count_plot.lung <- norm_viral_count_plot.lung + geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.75, position = position_dodge(0.8))
norm_viral_count_plot.lung <- norm_viral_count_plot.lung + labs(x="Day", y="log normalized counts") 
norm_viral_count_plot.lung <- norm_viral_count_plot.lung + theme_classic() 
norm_viral_count_plot.lung <- norm_viral_count_plot.lung + theme(text=element_text(size=24), title=element_text(size=24), legend.title=element_text(face="bold"), strip.text=element_text(face="bold")) 
norm_viral_count_plot.lung <- norm_viral_count_plot.lung + scale_y_continuous(limits=c(-2,7), breaks=c(-2,-1,0,1,2,3,4,5,6,7)) 
norm_viral_count_plot.lung <- norm_viral_count_plot.lung + facet_grid(.~Tissue) 
norm_viral_count_plot.lung <- norm_viral_count_plot.lung + stat_pvalue_manual(norm_viral_count.df_stats, label="{scales::pvalue(p.adj, accuracy=0.0001)}", inherit.aes=FALSE, size=8)
norm_viral_count_plot.lung
# dev.off()



########## Generate a table reporting transcriptome read alignment rates for airway epithelial cells and colon organoids (Supplementary Table 8)
samples.p <- samples
samples.p$star_stats_file <- paste("230323_STAR_alignment_stats_cell/", samples.p$Sample, "_S", rownames(samples.p), "_Log.final.out", sep="")
star_alignment_stats <- read.delim(file.path(dir, samples.p$star_stats_file[1]), header=FALSE)
colnames(star_alignment_stats) <- c("Parameter", samples.p$Sample[1])

for (a in 1:length(samples.p$star_stats_file)) {
  if (a != 1) {
    alignment_stats <- read.delim(file.path(dir, samples.p$star_stats_file[a]), header=FALSE)
    star_alignment_stats[samples.p$Sample[a]] <- alignment_stats[2]
  }
}

rownames(star_alignment_stats) <- star_alignment_stats[,1]
star_alignment_stats[,1] <- NULL
star_alignment_stats <- as.data.frame(t(star_alignment_stats))
star_alignment_stats$pct_mapped <- (as.numeric(star_alignment_stats[,8]) + as.numeric(star_alignment_stats[,23])) / as.numeric(star_alignment_stats[,5]) * 100
min(star_alignment_stats$pct_mapped)

star_alignment_stats.df <- star_alignment_stats[,c(5, 8, 23, 34)]
# write.table(star_alignment_stats.df, "table_s8_sample_quality_mapped_reads_hamster_lung.tsv", col.names=T, row.names=T, sep="\t", quote=F)



######################################## Analyze reads counts with separated viral gene reads ########################################
########## Read count table and record grouping information
count.df_deseq <- count.df
count.df_deseq <- count.df_deseq[,c(colnames(count.df_deseq)[1:6],passed_samples.v)]
count.df_no_viral_seq_deseq <- count.df_deseq[1:(nrow(count.df_deseq)-2),]



########## Run zFPKM analyses
##### Check distribution of FPKM in all samples
##### FPKM = Fragments Per Kilobase Million
fpkm.df <- count.df_no_viral_seq_deseq
fpkm.df[2:5] <- NULL
total_read_count_per_each_exp <- colSums(fpkm.df[,3:ncol(fpkm.df)])
fpkm.df[,3:ncol(fpkm.df)] <- (fpkm.df[,3:ncol(fpkm.df)]*1000000)/(fpkm.df[,2]/1000)
fpkm.df[,3:ncol(fpkm.df)] <- sweep(fpkm.df[,3:ncol(fpkm.df)],2,total_read_count_per_each_exp,FUN="/")
rownames(fpkm.df) <- fpkm.df$Geneid
fpkm.df[1:2] <- NULL
zfpkm.df <- zFPKM(fpkm.df)
zfpkm_plot <- zFPKMPlot(fpkm.df, FacetTitles=TRUE)
zfpkm_plot



######################################## Perform pairwise DEG analyses using Wald's test ########################################
########## Create FC and log2FC matrices
log2FC.mx <- data.frame()
padj.mx <- data.frame()
for (a in c("Hilum","Periphery")) {
  ##### Use "apeglm" method for LFC shrinkage
  for (b in c("2dpi","5dpi")) {
    filter_samples <- passed_samples %>% filter(Tissue==a, Day==b)
    if (b=="2dpi") {
      temp_saline_samples <- passed_samples %>% filter(Tissue==a, Virus=="Saline")
      temp_saline_samples$Day <- "2dpi"
      filter_samples <- rbind(filter_samples, temp_saline_samples)
    }
    dds_FC.mx <- DESeqDataSetFromMatrix(countData=count.df_no_viral_seq_deseq[,filter_samples$Sample], colData=filter_samples, design=~Virus)
    dds_FC.mx$Virus <- relevel(dds_FC.mx$Virus, "Saline")
    dds_FC.mx <- DESeq(dds_FC.mx)
    for (c in 2:length(resultsNames(dds_FC.mx))) {
      resLFC_FC.mx <- lfcShrink(dds_FC.mx, coef=resultsNames(dds_FC.mx)[c], type="apeglm")
      if (nrow(log2FC.mx)==0) {
        log2FC.mx <- data.frame(resLFC_FC.mx$log2FoldChange)
        rownames(log2FC.mx) <- rownames(resLFC_FC.mx)
        colnames(log2FC.mx) <- paste(a, str_split(str_split(resultsNames(dds_FC.mx)[c], "Virus_", simplify=T)[2], "_vs_", simplify=T)[1], b, sep=" ")
      } else {
        log2FC.mx[paste(a, str_split(str_split(resultsNames(dds_FC.mx)[c], "Virus_", simplify=T)[2], "_vs_", simplify=T)[1], b, sep=" ")] <- resLFC_FC.mx$log2FoldChange
      }
      if (nrow(padj.mx)==0) {
        padj.mx <- data.frame(resLFC_FC.mx$padj)
        rownames(padj.mx) <- rownames(resLFC_FC.mx)
        colnames(padj.mx) <- paste(a, str_split(str_split(resultsNames(dds_FC.mx)[c], "Virus_", simplify=T)[2], "_vs_", simplify=T)[1], b, sep=" ")
      } else {
        padj.mx[paste(a, str_split(str_split(resultsNames(dds_FC.mx)[c], "Virus_", simplify=T)[2], "_vs_", simplify=T)[1], b, sep=" ")] <- resLFC_FC.mx$padj
      }
    }
  }
}



########## Create PCA plots showing the expression change profiles of hamster lungs (Figure 4E)
log2FC.mx_no_any_na <- log2FC.mx[!apply(log2FC.mx, 1, function(x) any (is.na(x))),]
samples_FC_pca <- data.frame(colnames(log2FC.mx_no_any_na), row.names=colnames(log2FC.mx_no_any_na))
colnames(samples_FC_pca) <- c("Grouping")
samples_FC_pca$Tissue <- str_split(samples_FC_pca$Grouping, " ", simplify=T)[,1]
samples_FC_pca$Virus <- str_split(samples_FC_pca$Grouping, " ", simplify=T)[,2]
samples_FC_pca$Day <- str_split(samples_FC_pca$Grouping, " ", simplify=T)[,3]
samples_FC_pca$Virus <- factor(samples_FC_pca$Virus, levels=c("SC2","SC2dFCS","B236"))

df_pca <- pca(log2FC.mx_no_any_na, metadata=samples_FC_pca, center=TRUE, scale=TRUE) # Standardization = Center then Scale

# pdf("pca_fc_lungs.pdf", width=10, height=10)
log2_fc_pca_plot <- biplot(df_pca, colby="Day", shape="Virus", xlab=paste0("PC1, ", round(df_pca$variance["PC1"], digits=2), "% variation"), ylab=paste0("PC2, ", round(df_pca$variance["PC2"], digits=2), "% variation"), lab=NULL, labSize=3.5, pointSize=3, widthConnectors=0.5, legendPosition="right") + theme_classic() + theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.title=element_text(size=14), legend.text=element_text(size=14), strip.text=element_text(size=14), plot.margin = margin(0, 0.5, 0, 0.5, "cm")) + scale_shape_manual(values=15:17) + scale_colour_manual(values=c("lightblue","violet")) + geom_path(aes(group=df_pca$metadata$Virus), arrow=arrow(angle=22.5, ends="last", type="closed", length = unit(0.1, "inches")), linewidth=0.1) + coord_fixed() + facet_wrap(.~df_pca$metadata$Tissue)
log2_fc_pca_plot
#dev.off()



### Assign inactive/actively expressed genes
saline_hil_2d_inactive_genes <- rownames(zfpkm.df %>% select((passed_samples %>% filter(Grouping=="Saline_Hilum_5dpi"))$Sample) %>% filter(!rowMeans(.) > -3)) ##### Assume no significant changes in the Saline group from 2 to 5 d.p.i.
SC2_hil_2d_inactive_genes <- rownames(zfpkm.df %>% select((passed_samples %>% filter(Grouping=="SC2_Hilum_2dpi"))$Sample) %>% filter(!rowMeans(.) > -3))
SC2dFCS_hil_2d_inactive_genes <- rownames(zfpkm.df %>% select((passed_samples %>% filter(Grouping=="SC2dFCS_Hilum_2dpi"))$Sample) %>% filter(!rowMeans(.) > -3))
B236_hil_2d_inactive_genes <- rownames(zfpkm.df %>% select((passed_samples %>% filter(Grouping=="B236_Hilum_2dpi"))$Sample) %>% filter(!rowMeans(.) > -3))

saline_B236_hil_2d_inactive_genes <- intersect(saline_hil_2d_inactive_genes, B236_hil_2d_inactive_genes)
saline_SC2_hil_2d_inactive_genes <- intersect(saline_hil_2d_inactive_genes, SC2_hil_2d_inactive_genes)
saline_SC2dFCS_hil_2d_inactive_genes <- intersect(saline_hil_2d_inactive_genes, SC2dFCS_hil_2d_inactive_genes)

saline_ppr_2d_inactive_genes <- rownames(zfpkm.df %>% select((passed_samples %>% filter(Grouping=="Saline_Periphery_5dpi"))$Sample) %>% filter(!rowMeans(.) > -3)) ##### Assume no significant changes in the Saline group from 2 to 5 d.p.i.
SC2_ppr_2d_inactive_genes <- rownames(zfpkm.df %>% select((passed_samples %>% filter(Grouping=="SC2_Periphery_2dpi"))$Sample) %>% filter(!rowMeans(.) > -3))
SC2dFCS_ppr_2d_inactive_genes <- rownames(zfpkm.df %>% select((passed_samples %>% filter(Grouping=="SC2dFCS_Periphery_2dpi"))$Sample) %>% filter(!rowMeans(.) > -3))
B236_ppr_2d_inactive_genes <- rownames(zfpkm.df %>% select((passed_samples %>% filter(Grouping=="B236_Periphery_2dpi"))$Sample) %>% filter(!rowMeans(.) > -3))

saline_SC2_ppr_2d_inactive_genes <- intersect(saline_ppr_2d_inactive_genes, SC2_ppr_2d_inactive_genes)
saline_SC2dFCS_ppr_2d_inactive_genes <- intersect(saline_ppr_2d_inactive_genes, SC2dFCS_ppr_2d_inactive_genes)
saline_B236_ppr_2d_inactive_genes <- intersect(saline_ppr_2d_inactive_genes, B236_ppr_2d_inactive_genes)

saline_hil_5d_inactive_genes <- rownames(zfpkm.df %>% select((passed_samples %>% filter(Grouping=="Saline_Hilum_5dpi"))$Sample) %>% filter(!rowMeans(.) > -3))
SC2_hil_5d_inactive_genes <- rownames(zfpkm.df %>% select((passed_samples %>% filter(Grouping=="SC2_Hilum_5dpi"))$Sample) %>% filter(!rowMeans(.) > -3))
SC2dFCS_hil_5d_inactive_genes <- rownames(zfpkm.df %>% select((passed_samples %>% filter(Grouping=="SC2dFCS_Hilum_5dpi"))$Sample) %>% filter(!rowMeans(.) > -3))
B236_hil_5d_inactive_genes <- rownames(zfpkm.df %>% select((passed_samples %>% filter(Grouping=="B236_Hilum_5dpi"))$Sample) %>% filter(!rowMeans(.) > -3))

saline_SC2_hil_5d_inactive_genes <- intersect(saline_hil_5d_inactive_genes, SC2_hil_5d_inactive_genes)
saline_SC2dFCS_hil_5d_inactive_genes <- intersect(saline_hil_5d_inactive_genes, SC2dFCS_hil_5d_inactive_genes)
saline_B236_hil_5d_inactive_genes <- intersect(saline_hil_5d_inactive_genes, B236_hil_5d_inactive_genes)

saline_ppr_5d_inactive_genes <- rownames(zfpkm.df %>% select((passed_samples %>% filter(Grouping=="Saline_Periphery_5dpi"))$Sample) %>% filter(!rowMeans(.) > -3))
SC2_ppr_5d_inactive_genes <- rownames(zfpkm.df %>% select((passed_samples %>% filter(Grouping=="SC2_Periphery_5dpi"))$Sample) %>% filter(!rowMeans(.) > -3))
SC2dFCS_ppr_5d_inactive_genes <- rownames(zfpkm.df %>% select((passed_samples %>% filter(Grouping=="SC2dFCS_Periphery_5dpi"))$Sample) %>% filter(!rowMeans(.) > -3))
B236_ppr_5d_inactive_genes <- rownames(zfpkm.df %>% select((passed_samples %>% filter(Grouping=="B236_Periphery_5dpi"))$Sample) %>% filter(!rowMeans(.) > -3))

saline_B236_ppr_5d_inactive_genes <- intersect(saline_ppr_5d_inactive_genes, B236_ppr_5d_inactive_genes)
saline_SC2_ppr_5d_inactive_genes <- intersect(saline_ppr_5d_inactive_genes, SC2_ppr_5d_inactive_genes)
saline_SC2dFCS_ppr_5d_inactive_genes <- intersect(saline_ppr_5d_inactive_genes, SC2dFCS_ppr_5d_inactive_genes)



##### Create a heatmap for log2 expression fold change of selected host genes across all samples
log2FC.mx_no_all_na <- log2FC.mx[!apply(log2FC.mx, 1, function(x) all (is.na(x))),]
padj.mx_no_all_na <- padj.mx[!apply(log2FC.mx, 1, function(x) all (is.na(x))),]

inactive_gene_list <- list(saline_B236_hil_2d_inactive_genes, saline_SC2_hil_2d_inactive_genes, saline_SC2dFCS_hil_2d_inactive_genes, 
                           saline_B236_hil_5d_inactive_genes, saline_SC2_hil_5d_inactive_genes, saline_SC2dFCS_hil_5d_inactive_genes, 
                           saline_B236_ppr_2d_inactive_genes, saline_SC2_ppr_2d_inactive_genes, saline_SC2dFCS_ppr_2d_inactive_genes, 
                           saline_B236_ppr_5d_inactive_genes, saline_SC2_ppr_5d_inactive_genes, saline_SC2dFCS_ppr_5d_inactive_genes)
inactive.mx_no_all_na <- padj.mx_no_all_na

for (a in 1:ncol(inactive.mx_no_all_na)) {inactive.mx_no_all_na[,a] <- as.integer(rownames(inactive.mx_no_all_na) %in% inactive_gene_list[[a]])}



########## Create a heatmap for log2FC of selected ISGs in the infected hamster lungs (Figure 3H)
human_hamster_ortholog <- read.table(file.path(dir, "symbol_human_hamster_ortholog.txt"), header=TRUE) # Use the ortholog information to convert gene names of human and hamster

selected_isg_gene_list <- sort(c("IFNA1","IFNB1","IFNAR1","IFNAR2","IFITM3","ISG15","JAK1","STAT2","IRF3","IRF7","NFKB1","MAVS","DDX58","TBK1","OAS1","PLSCR1","LY6E","IFIH1")) # RIGI is DDX58 in this GFF3 file
selected_isg_gene_list <- human_hamster_ortholog[match(selected_isg_gene_list, human_hamster_ortholog$symbol.human),]$symbol.hamster
selected_isg_gene_list[3] <- "LOC101840863" ### LOC101840863 (IFITM3)
selected_isg_gene_list[15] <- "LOC106021992" ### LOC106021992 (OAS1)
selected_isg_gene_list[16] <- "LOC101834040" ### LOC101834040 (PLSCR1)

selected_isg_gene_list <- selected_isg_gene_list[which(!is.na(selected_isg_gene_list))]

# pdf("selected_isg_deg_heatmap_hamster_lung.pdf", width=7.501, height=5.455)
selected_isg_log2FC_heatmap_ha <- rowAnnotation(Virus=rep(c("SC2","SC2dFCS","B236"), 4), Day=rep(c(rep("2dpi", 3), rep("5dpi", 3)), 2), 
                                                col=list(Virus=c("SC2"="yellow", "SC2dFCS"="navyblue", "B236"="hotpink"), Day=c("2dpi"="lightblue", "5dpi"="violet")), 
                                                annotation_name_gp=gpar(fontsize=16, fontface="bold"), annotation_legend_param=list(title_gp=gpar(fontsize=16, fontface="bold"), labels_gp=gpar(fontsize=16)))
selected_isg_log2FC_heatmap <- Heatmap(t(log2FC.mx[selected_isg_gene_list,c(2,3,1,5,6,4,8,9,7,11,12,10)]), left_annotation=selected_isg_log2FC_heatmap_ha, 
                                       show_row_names=FALSE, row_title_rot=90, row_split=factor(c(rep("Hilum", 6), rep("Periphery", 6)), levels=c("Hilum","Periphery")), 
                                       row_names_gp=gpar(fontsize=24), row_gap=unit(2.5, "mm"), cluster_rows=FALSE, row_title_gp=gpar(fontface="bold", fontsize=18), clustering_method_rows="ward.D2", 
                                       column_names_gp=gpar(fontsize=14, fontface="italic"), column_gap=unit(2.5, "mm"), cluster_columns=FALSE, column_title_gp=gpar(fontface="bold", fontsize=18), clustering_method_columns="ward.D2", 
                                       col=colorRamp2(c(-2, 0, 2), c("blue","white","red")), na_col="lightgrey", heatmap_legend_param=list(title="log2FC", title_gp=gpar(fontface="bold", fontsize=16), labels_gp=gpar(fontsize=16)), 
                                       cell_fun=function(j, i, x, y, width, height, fill) {if (!is.na(padj.mx_no_all_na[selected_isg_gene_list,c(2,3,1,5,6,4,8,9,7,11,12,10)][j,i]) 
                                                                                               && padj.mx_no_all_na[selected_isg_gene_list,c(2,3,1,5,6,4,8,9,7,11,12,10)][j,i]<0.05 
                                                                                               && inactive.mx_no_all_na[selected_isg_gene_list,c(2,3,1,5,6,4,8,9,7,11,12,10)][j,i]==0) grid.text("*", x, y, gp=gpar(fontsize=16), vjust=+0.75)})
selected_isg_log2FC_heatmap
# dev.off()



########## Cellular response to type I interferon (GO:0071357) # identical to Type I interferon signaling pathway (GO:0060337) in 2021
GO_0071357_gene_list <- c("MX1","RSAD2","MX2","OAS1","TYK2","BST2","IRF1","OAS2","IRF4","IRF5","OAS3","IFI27","IRF2","IRF3","IRF8","IRF9","IRF6","IRF7","IFITM2","IFITM3","SP100","IFITM1","IFI6","ADAR","IRAK1","IP6K2","GBP2","EGR1","IFNA14","STAT1","IFNA16","STAT2","IFNA17","HLA-C","HLA-A","ISG15","HLA-B","HLA-G","HLA-E","HLA-F","ISG20","MYD88","XAF1","IFNAR1","IFNA10","IFNAR2","IFNA6","IFNA5","IFNA8","RNASEL","IFNA7","IFNA2","IFNA1","IFNA4","IFI35","IFIT5","IFIT2","IFIT1","OASL","SAMHD1","IFIT3","PSMB8","IFNB1","JAK1","IFNA21")
GO_0071357_gene_list <- sort(GO_0071357_gene_list)
GO_0071357_gene_list <- human_hamster_ortholog[match(GO_0071357_gene_list, human_hamster_ortholog$symbol.human),]$symbol.hamster

GO_0071357_gene_list[4] <- "LOC101824821" ### GBP2
GO_0071357_gene_list[14] <- "LOC106020791" ### IFIT1
GO_0071357_gene_list[18] <- "LOC106021464" ### IFITM1
GO_0071357_gene_list[19] <- "LOC121132778" ### IFITM2
GO_0071357_gene_list[20] <- "LOC101840863" ### IFITM3
GO_0071357_gene_list[53] <- "LOC106021992" ### OAS1
GO_0071357_gene_list[61] <- "LOC101840585" ### SP100

GO_0071357_gene_list <- GO_0071357_gene_list[which(!is.na(GO_0071357_gene_list))]

log2FC.mx_no_all_na.lung.isg <- as.data.frame(as.table(as.matrix(log2FC.mx_no_all_na[GO_0071357_gene_list,])))
colnames(log2FC.mx_no_all_na.lung.isg) <- c("Gene","Grouping","log2FC")

log2FC.mx_no_all_na.lung.isg$Virus <- str_split(log2FC.mx_no_all_na.lung.isg$Grouping, " ", simplify=T)[,2]
log2FC.mx_no_all_na.lung.isg$Tissue <- str_split(log2FC.mx_no_all_na.lung.isg$Grouping, " ", simplify=T)[,1]
log2FC.mx_no_all_na.lung.isg$Day <- str_split(log2FC.mx_no_all_na.lung.isg$Grouping, " ", simplify=T)[,3]

padj.mx_no_all_na.lung.isg <- as.data.frame(as.table(as.matrix(padj.mx_no_all_na[GO_0071357_gene_list,])))
colnames(padj.mx_no_all_na.lung.isg) <- c("Gene","Grouping","padj")
padj.mx_no_all_na.lung.isg$Virus <- rep(c(rep("B236",length(GO_0071357_gene_list)), rep("SC2",length(GO_0071357_gene_list)), rep("SC2dFCS",length(GO_0071357_gene_list))), 2)
padj.mx_no_all_na.lung.isg$Tissue <- c(rep("Hilum",3*length(GO_0071357_gene_list)), rep("Periphery",3*length(GO_0071357_gene_list)))

inactive.mx_no_all_na.lung.isg <- as.data.frame(as.table(as.matrix(inactive.mx_no_all_na[GO_0071357_gene_list,])))
colnames(inactive.mx_no_all_na.lung.isg) <- c("Gene","Grouping","inactive")
inactive.mx_no_all_na.lung.isg$Virus <- c(rep("B236",length(GO_0071357_gene_list)), rep("SC2",length(GO_0071357_gene_list)), rep("SC2dFCS",length(GO_0071357_gene_list)))

log2FC.mx_no_all_na.lung.isg$padj <- padj.mx_no_all_na.lung.isg$padj
log2FC.mx_no_all_na.lung.isg$inactive <- inactive.mx_no_all_na.lung.isg$inactive
log2FC.mx_no_all_na.lung.isg$DEG <- NA
for (a in 1:nrow(log2FC.mx_no_all_na.lung.isg)) {if (isTRUE(log2FC.mx_no_all_na.lung.isg[a,]$padj < 0.05 && log2FC.mx_no_all_na.lung.isg[a,]$inactive=="0")==TRUE) {log2FC.mx_no_all_na.lung.isg[a,]$DEG <- "DEG"} else {log2FC.mx_no_all_na.lung.isg[a,]$DEG <- "non-DEG"}}

log2FC.mx_no_all_na.lung.isg$Virus <- factor(log2FC.mx_no_all_na.lung.isg$Virus, levels=c("SC2","SC2dFCS","B236"))

log2FC.mx_no_all_na.lung.isg %>% filter(log2FC!=-Inf)  %>% group_by(Virus, Tissue, Day) %>% shapiro_test(log2FC)
log2FC.mx_no_all_na.lung.isg %>% filter(log2FC!=-Inf)  %>% group_by(Day, Tissue) %>% levene_test(log2FC~Virus)

### Use Kruskal-Wallist test since all assumptions are violated
res.kruskal <- log2FC.mx_no_all_na.lung.isg %>% filter(log2FC!=-Inf)  %>% group_by(Tissue, Day) %>% kruskal_test(log2FC~Virus)
res.kruskal

log2FC.mx_no_all_na.lung.isg_stats <- log2FC.mx_no_all_na.lung.isg %>% filter(log2FC!=-Inf)  %>% group_by(Tissue, Day) %>% dunn_test(log2FC~Virus, p.adjust.method="fdr") # FDR (alias of BH)
log2FC.mx_no_all_na.lung.isg_stats <- log2FC.mx_no_all_na.lung.isg_stats %>% add_xy_position(x="Day", fun="max", dodge=0.8)
log2FC.mx_no_all_na.lung.isg_stats

# pdf("violin_ifn_hamster_lung.pdf", width=9, height=6)
ifn_violin_plot <- ggplot(log2FC.mx_no_all_na.lung.isg %>% filter(log2FC!=-Inf) , aes(x=Day, y=log2FC, fill=Virus), alpha=0.5)
ifn_violin_plot <- ifn_violin_plot + geom_jitter(data=filter(log2FC.mx_no_all_na.lung.isg %>% filter(log2FC!=-Inf) , DEG=="DEG"), aes(col=DEG), alpha=0.5, shape=16, size=2, position=position_jitterdodge(jitter.width=0.5, dodge.width=0.8)) 
ifn_violin_plot <- ifn_violin_plot + geom_jitter(data=filter(log2FC.mx_no_all_na.lung.isg %>% filter(log2FC!=-Inf) , DEG=="non-DEG"), aes(col=DEG), alpha=0.5, shape=16, size=2, position=position_jitterdodge(jitter.width=0.5, dodge.width=0.8)) 
ifn_violin_plot <- ifn_violin_plot + geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.75, position = position_dodge(0.8)) 
ifn_violin_plot <- ifn_violin_plot + labs(title="GO:0071357", subtitle="cellular response to type I interferon", x="Day (d.p.i.)", y=bquote(log[2]~FC)) 
ifn_violin_plot <- ifn_violin_plot + guides() 
ifn_violin_plot <- ifn_violin_plot + theme_classic() 
ifn_violin_plot <- ifn_violin_plot + theme(text=element_text(size=16), plot.title=element_text(size=16), plot.subtitle=element_text(size=14), axis.title=element_text(size=16)) 
ifn_violin_plot <- ifn_violin_plot + scale_y_continuous(limits=c(-2,17), breaks=c(-2,0,2,4,6,8,10,12,14,16)) 
ifn_violin_plot <- ifn_violin_plot + facet_grid(.~Tissue) 
ifn_violin_plot <- ifn_violin_plot + stat_pvalue_manual(log2FC.mx_no_all_na.lung.isg_stats, label="{scales::pvalue(p.adj, accuracy=0.0001)}", inherit.aes=FALSE, size=4)
ifn_violin_plot
# dev.off()


######################################## Perform DEG analyses using Wald's test and likelihood ratio test ########################################
######################################## Hilum (2 d.p.i.) ########################################
########## 1. Analyze DEGs in hamster lung hila at 2 d.p.i.
hil_filter_samples <- rbind(passed_samples %>% filter(Tissue=="Hilum"), passed_samples %>% filter(Tissue=="Hilum", Virus=="Saline") %>% mutate(Day="2dpi", Grouping="Saline_Hilum_2dpi"))
hil_2d_filter_samples <- hil_filter_samples %>% filter(Day=="2dpi")
dds_hil_2d <- DESeqDataSetFromMatrix(countData=count.df_no_viral_seq_deseq[,hil_2d_filter_samples$Sample], colData=hil_2d_filter_samples, design=~Virus)
dds_hil_2d$Virus <- relevel(dds_hil_2d$Virus, "Saline")
dds_hil_2d_wald <- DESeq(dds_hil_2d) # Use for pairwise comparison
dds_hil_2d_lrt <- DESeq(dds_hil_2d, test="LRT", reduced=~1) # Use for multiple group comparison

####### 1.1 SC2 vs mock-infected
resLFC_B236_hil_2d_wald <- lfcShrink(dds_hil_2d_wald, coef="Virus_B236_vs_Saline", type="apeglm")
sorted_resLFC_B236_hil_2d_wald <- resLFC_B236_hil_2d_wald[order(resLFC_B236_hil_2d_wald$padj),]
sorted_sig_resLFC_B236_hil_2d_wald <- sorted_resLFC_B236_hil_2d_wald[which(sorted_resLFC_B236_hil_2d_wald$padj < 0.05),]
filtered_sorted_sig_resLFC_B236_hil_2d_wald <- sorted_sig_resLFC_B236_hil_2d_wald[which(!rownames(sorted_sig_resLFC_B236_hil_2d_wald) %in% saline_B236_hil_2d_inactive_genes),]
B236_hil_2d_all_down_gene <- rownames(filtered_sorted_sig_resLFC_B236_hil_2d_wald[which(filtered_sorted_sig_resLFC_B236_hil_2d_wald$log2FoldChange < 0),])
B236_hil_2d_all_up_gene <- rownames(filtered_sorted_sig_resLFC_B236_hil_2d_wald[which(filtered_sorted_sig_resLFC_B236_hil_2d_wald$log2FoldChange > 0),])
summary(filtered_sorted_sig_resLFC_B236_hil_2d_wald)

####### 1.2 SC2dFCS vs mock-infected
resLFC_SC2_hil_2d_wald <- lfcShrink(dds_hil_2d_wald, coef="Virus_SC2_vs_Saline", type="apeglm")
sorted_resLFC_SC2_hil_2d_wald <- resLFC_SC2_hil_2d_wald[order(resLFC_SC2_hil_2d_wald$padj),]
sorted_sig_resLFC_SC2_hil_2d_wald <- sorted_resLFC_SC2_hil_2d_wald[which(sorted_resLFC_SC2_hil_2d_wald$padj < 0.05),]
filtered_sorted_sig_resLFC_SC2_hil_2d_wald <- sorted_sig_resLFC_SC2_hil_2d_wald[which(!rownames(sorted_sig_resLFC_SC2_hil_2d_wald) %in% saline_SC2_hil_2d_inactive_genes),]
SC2_hil_2d_all_down_gene <- rownames(filtered_sorted_sig_resLFC_SC2_hil_2d_wald[which(filtered_sorted_sig_resLFC_SC2_hil_2d_wald$log2FoldChange < 0),])
SC2_hil_2d_all_up_gene <- rownames(filtered_sorted_sig_resLFC_SC2_hil_2d_wald[which(filtered_sorted_sig_resLFC_SC2_hil_2d_wald$log2FoldChange > 0),])
summary(filtered_sorted_sig_resLFC_SC2_hil_2d_wald)

####### 1.3 B236 vs mock-infected
resLFC_SC2dFCS_hil_2d_wald <- lfcShrink(dds_hil_2d_wald, coef="Virus_SC2dFCS_vs_Saline", type="apeglm")
sorted_resLFC_SC2dFCS_hil_2d_wald <- resLFC_SC2dFCS_hil_2d_wald[order(resLFC_SC2dFCS_hil_2d_wald$padj),]
sorted_sig_resLFC_SC2dFCS_hil_2d_wald <- sorted_resLFC_SC2dFCS_hil_2d_wald[which(sorted_resLFC_SC2dFCS_hil_2d_wald$padj < 0.05),]
filtered_sorted_sig_resLFC_SC2dFCS_hil_2d_wald <- sorted_sig_resLFC_SC2dFCS_hil_2d_wald[which(!rownames(sorted_sig_resLFC_SC2dFCS_hil_2d_wald) %in% saline_SC2dFCS_hil_2d_inactive_genes),]
SC2dFCS_hil_2d_all_down_gene <- rownames(filtered_sorted_sig_resLFC_SC2dFCS_hil_2d_wald[which(filtered_sorted_sig_resLFC_SC2dFCS_hil_2d_wald$log2FoldChange < 0),])
SC2dFCS_hil_2d_all_up_gene <- rownames(filtered_sorted_sig_resLFC_SC2dFCS_hil_2d_wald[which(filtered_sorted_sig_resLFC_SC2dFCS_hil_2d_wald$log2FoldChange > 0),])
summary(filtered_sorted_sig_resLFC_SC2dFCS_hil_2d_wald)

####### 1.4 Summarize the number of DEG
hil_2d_deg_summary.df <- data.frame(rep(c("B236","SC2","SC2dFCS"), 2), c(length(B236_hil_2d_all_up_gene), length(SC2_hil_2d_all_up_gene), length(SC2dFCS_hil_2d_all_up_gene), 
                                                                         length(B236_hil_2d_all_down_gene), length(SC2_hil_2d_all_down_gene), length(SC2dFCS_hil_2d_all_down_gene)), c(rep(c("Upregulated"), 3), rep(c("Downregulated"), 3)))
colnames(hil_2d_deg_summary.df) <- c("Virus","Genes","Direction")
hil_2d_deg_summary.df$Grouping <- "Hilum_2dpi"

####### 1.5 Analyze the difference between B236 infection and SC2 or SC2dFCS infection in lung hila at 2 d.p.i.
hil_2d_B236_SC2_filter_samples <- rbind(passed_samples %>% filter(Tissue=="Hilum", Day=="2dpi", !Virus=="SC2dFCS"), 
                                        passed_samples %>% filter(Tissue=="Hilum", Virus=="Saline") %>% mutate(Day="2dpi", Grouping="Saline_Hilum_2dpi"))
dds_hil_2d_B236_SC2 <- DESeqDataSetFromMatrix(countData=count.df_no_viral_seq_deseq[,hil_2d_B236_SC2_filter_samples$Sample], 
                                              colData=hil_2d_B236_SC2_filter_samples, design=~Virus)
dds_hil_2d_B236_SC2$Virus <- relevel(dds_hil_2d_B236_SC2$Virus, "Saline")
dds_hil_2d_B236_SC2_wald <- DESeq(dds_hil_2d_B236_SC2) # Use for pairwise comparison
dds_hil_2d_B236_SC2_lrt <- DESeq(dds_hil_2d_B236_SC2, test="LRT", reduced=~1) # Use for multiple group comparison

resLFC_B236_SC2_hil_2d_lrt <- lfcShrink(dds_hil_2d_B236_SC2_lrt, coef="Virus_B236_vs_Saline", type="apeglm")
sorted_resLFC_B236_SC2_hil_2d_lrt <- resLFC_B236_SC2_hil_2d_lrt[order(resLFC_B236_SC2_hil_2d_lrt$padj),]
sorted_sig_resLFC_B236_SC2_hil_2d_lrt <- sorted_resLFC_B236_SC2_hil_2d_lrt[which(sorted_resLFC_B236_SC2_hil_2d_lrt$padj < 0.05),]
summary(sorted_sig_resLFC_B236_SC2_hil_2d_lrt)
hil_2d_deg_list_lrt <- rownames(sorted_sig_resLFC_B236_SC2_hil_2d_lrt)

##### 1.5.1 Calculate disco score
hil_2d_disco <- data.frame(pvalue_B236_hil=resLFC_B236_hil_2d_wald$pvalue, padj_B236_hil=resLFC_B236_hil_2d_wald$padj, log2FC_B236_hil=resLFC_B236_hil_2d_wald$log2FoldChange, 
                           pvalue_SC2_hil=resLFC_SC2_hil_2d_wald$pvalue, padj_SC2_hil=resLFC_SC2_hil_2d_wald$padj, log2FC_SC2_hil=resLFC_SC2_hil_2d_wald$log2FoldChange, row.names=rownames(resLFC_B236_hil_2d_wald))
hil_2d_disco[,c("category","B236_plot_category","SC2_plot_category","pool_plot_category","disco_score")] <- NA

##### 1.5.2 Assign gene expression category based on directions of fold change
for (i in 1:nrow(hil_2d_disco)) {
  if (is.na(hil_2d_disco[i,]$pvalue_B236_hil) || is.na(hil_2d_disco[i,]$padj_B236_hil) || rownames(hil_2d_disco[i,]) %in% saline_B236_hil_2d_inactive_genes || hil_2d_disco[i,]$padj_B236_hil>=0.05) {
    B236_hil_2d_direction <- "X"} else if (hil_2d_disco[i,]$log2FC_B236_hil<0) {B236_hil_2d_direction="D"} else {B236_hil_2d_direction="U"}
  if (is.na(hil_2d_disco[i,]$pvalue_SC2_hil) || is.na(hil_2d_disco[i,]$padj_SC2_hil) || rownames(hil_2d_disco[i,]) %in% saline_SC2_hil_2d_inactive_genes || hil_2d_disco[i,]$padj_SC2_hil>=0.05) {
    SC2_hil_2d_direction <- "X"} else if (hil_2d_disco[i,]$log2FC_SC2_hil<0) {SC2_hil_2d_direction="D"} else {SC2_hil_2d_direction="U"}
  hil_2d_disco[i,]$category <- paste(B236_hil_2d_direction, SC2_hil_2d_direction, sep="")
  if (!rownames(hil_2d_disco[i,]) %in% hil_2d_deg_list_lrt) {hil_2d_disco[i,]$category <- "XX"}
  if (hil_2d_disco[i,]$category %in% c("DD")) {hil_2d_disco[i,]$B236_plot_category <- "B236-Down/SC2-Down"}
  if (hil_2d_disco[i,]$category %in% c("UU")) {hil_2d_disco[i,]$B236_plot_category <- "B236-Up/SC2-Up"}
  if (hil_2d_disco[i,]$category %in% c("DU","DX")) {hil_2d_disco[i,]$B236_plot_category <- "B236-Down/SC2-notDown"}
  if (hil_2d_disco[i,]$category %in% c("UD","UX")) {hil_2d_disco[i,]$B236_plot_category <- "B236-Up/SC2-notUp"}
  if (is.na(hil_2d_disco[i,]$B236_plot_category)) {hil_2d_disco[i,]$B236_plot_category <- "OOI"}
  if (hil_2d_disco[i,]$category %in% c("DD")) {hil_2d_disco[i,]$SC2_plot_category <- "B236-Down/SC2-Down"}
  if (hil_2d_disco[i,]$category %in% c("UU")) {hil_2d_disco[i,]$SC2_plot_category <- "B236-Up/SC2-Up"}
  if (hil_2d_disco[i,]$category %in% c("UD","XD")) {hil_2d_disco[i,]$SC2_plot_category <- "B236-notDown/SC2-Down"}
  if (hil_2d_disco[i,]$category %in% c("DU","XU")) {hil_2d_disco[i,]$SC2_plot_category <- "B236-notUp/SC2-Up"}
  if (is.na(hil_2d_disco[i,]$SC2_plot_category)) {hil_2d_disco[i,]$SC2_plot_category <- "OOI"}
}

hil_2d_disco$disco_score <- abs(hil_2d_disco$log2FC_B236_hil * hil_2d_disco$log2FC_SC2_hil * (log10(hil_2d_disco$pvalue_B236_hil) + log10(hil_2d_disco$pvalue_SC2_hil)))
sorted_hil_2d_disco <- hil_2d_disco[order(hil_2d_disco$disco_score, decreasing=TRUE),]

sorted_hil_2d_disco$human_gene_name <- human_hamster_ortholog[match(rownames(sorted_hil_2d_disco), human_hamster_ortholog$symbol.hamster),]$symbol.human

sorted_top_B236_hil_2d_disco <- sorted_hil_2d_disco %>% filter(!is.na(human_gene_name))
sorted_top_B236_hil_2d_disco <- Reduce(rbind, by(sorted_top_B236_hil_2d_disco, sorted_top_B236_hil_2d_disco["B236_plot_category"], head, n=100))
sorted_top_SC2_hil_2d_disco <- sorted_hil_2d_disco %>% filter(!is.na(human_gene_name))
sorted_top_SC2_hil_2d_disco <- Reduce(rbind, by(sorted_top_SC2_hil_2d_disco, sorted_top_SC2_hil_2d_disco["SC2_plot_category"], head, n=100))

######################################## Hilum (5 d.p.i.) ########################################
########## 2. Analyze DEGs in hamster lung hila at 5 d.p.i.
hil_5d_filter_samples <- hil_filter_samples %>% filter(Day=="5dpi")
dds_hil_5d <- DESeqDataSetFromMatrix(countData=count.df_no_viral_seq_deseq[,hil_5d_filter_samples$Sample], colData=hil_5d_filter_samples, design=~Virus)
dds_hil_5d$Virus <- relevel(dds_hil_5d$Virus, "Saline")
dds_hil_5d_wald <- DESeq(dds_hil_5d) # Use for pairwise comparison
dds_hil_5d_lrt <- DESeq(dds_hil_5d, test="LRT", reduced=~1) # Use for multiple group comparison

####### 2.1 SC2 vs mock-infected
resLFC_B236_hil_5d_wald <- lfcShrink(dds_hil_5d_wald, coef="Virus_B236_vs_Saline", type="apeglm")
sorted_resLFC_B236_hil_5d_wald <- resLFC_B236_hil_5d_wald[order(resLFC_B236_hil_5d_wald$padj),]
sorted_sig_resLFC_B236_hil_5d_wald <- sorted_resLFC_B236_hil_5d_wald[which(sorted_resLFC_B236_hil_5d_wald$padj < 0.05),]
filtered_sorted_sig_resLFC_B236_hil_5d_wald <- sorted_sig_resLFC_B236_hil_5d_wald[which(!rownames(sorted_sig_resLFC_B236_hil_5d_wald) %in% saline_B236_hil_5d_inactive_genes),]
B236_hil_5d_all_down_gene <- rownames(filtered_sorted_sig_resLFC_B236_hil_5d_wald[which(filtered_sorted_sig_resLFC_B236_hil_5d_wald$log2FoldChange < 0),])
B236_hil_5d_all_up_gene <- rownames(filtered_sorted_sig_resLFC_B236_hil_5d_wald[which(filtered_sorted_sig_resLFC_B236_hil_5d_wald$log2FoldChange > 0),])
summary(filtered_sorted_sig_resLFC_B236_hil_5d_wald)

####### 2.2 SC2dFCS vs mock-infected
resLFC_SC2_hil_5d_wald <- lfcShrink(dds_hil_5d_wald, coef="Virus_SC2_vs_Saline", type="apeglm")
sorted_resLFC_SC2_hil_5d_wald <- resLFC_SC2_hil_5d_wald[order(resLFC_SC2_hil_5d_wald$padj),]
sorted_sig_resLFC_SC2_hil_5d_wald <- sorted_resLFC_SC2_hil_5d_wald[which(sorted_resLFC_SC2_hil_5d_wald$padj < 0.05),]
filtered_sorted_sig_resLFC_SC2_hil_5d_wald <- sorted_sig_resLFC_SC2_hil_5d_wald[which(!rownames(sorted_sig_resLFC_SC2_hil_5d_wald) %in% saline_SC2_hil_5d_inactive_genes),]
SC2_hil_5d_all_down_gene <- rownames(filtered_sorted_sig_resLFC_SC2_hil_5d_wald[which(filtered_sorted_sig_resLFC_SC2_hil_5d_wald$log2FoldChange < 0),])
SC2_hil_5d_all_up_gene <- rownames(filtered_sorted_sig_resLFC_SC2_hil_5d_wald[which(filtered_sorted_sig_resLFC_SC2_hil_5d_wald$log2FoldChange > 0),])
summary(filtered_sorted_sig_resLFC_SC2_hil_5d_wald)

####### 2.3 B236 vs mock-infected
resLFC_SC2dFCS_hil_5d_wald <- lfcShrink(dds_hil_5d_wald, coef="Virus_SC2dFCS_vs_Saline", type="apeglm")
sorted_resLFC_SC2dFCS_hil_5d_wald <- resLFC_SC2dFCS_hil_5d_wald[order(resLFC_SC2dFCS_hil_5d_wald$padj),]
sorted_sig_resLFC_SC2dFCS_hil_5d_wald <- sorted_resLFC_SC2dFCS_hil_5d_wald[which(sorted_resLFC_SC2dFCS_hil_5d_wald$padj < 0.05),]
filtered_sorted_sig_resLFC_SC2dFCS_hil_5d_wald <- sorted_sig_resLFC_SC2dFCS_hil_5d_wald[which(!rownames(sorted_sig_resLFC_SC2dFCS_hil_5d_wald) %in% saline_SC2dFCS_hil_5d_inactive_genes),]
SC2dFCS_hil_5d_all_down_gene <- rownames(filtered_sorted_sig_resLFC_SC2dFCS_hil_5d_wald[which(filtered_sorted_sig_resLFC_SC2dFCS_hil_5d_wald$log2FoldChange < 0),])
SC2dFCS_hil_5d_all_up_gene <- rownames(filtered_sorted_sig_resLFC_SC2dFCS_hil_5d_wald[which(filtered_sorted_sig_resLFC_SC2dFCS_hil_5d_wald$log2FoldChange > 0),])
summary(filtered_sorted_sig_resLFC_SC2dFCS_hil_5d_wald)

####### 2.4 Summarize the number of DEG
hil_5d_deg_summary.df <- data.frame(rep(c("B236","SC2","SC2dFCS"), 2), c(length(B236_hil_5d_all_up_gene), length(SC2_hil_5d_all_up_gene), length(SC2dFCS_hil_5d_all_up_gene),
                                                                         length(B236_hil_5d_all_down_gene), length(SC2_hil_5d_all_down_gene), length(SC2dFCS_hil_5d_all_down_gene)), c(rep(c("Upregulated"), 3), rep(c("Downregulated"), 3)))
colnames(hil_5d_deg_summary.df) <- c("Virus","Genes","Direction")
hil_5d_deg_summary.df$Grouping <- "Hilum_5dpi"

####### 2.5 Analyze the difference between B236 infection and SC2 or SC2dFCS infection in lung hila at 5 d.p.i.
hil_5d_B236_SC2_filter_samples <- rbind(passed_samples %>% filter(Tissue=="Hilum", Day=="5dpi", !Virus=="SC2dFCS"), 
                                        passed_samples %>% filter(Tissue=="Hilum", Virus=="Saline") %>% mutate(Day="5dpi", Grouping="Saline_Hilum_5dpi"))
dds_hil_5d_B236_SC2 <- DESeqDataSetFromMatrix(countData=count.df_no_viral_seq_deseq[,hil_5d_B236_SC2_filter_samples$Sample], colData=hil_5d_B236_SC2_filter_samples, design=~Virus)

dds_hil_5d_B236_SC2$Virus <- relevel(dds_hil_5d_B236_SC2$Virus, "Saline")
dds_hil_5d_B236_SC2_wald <- DESeq(dds_hil_5d_B236_SC2) # Use for pairwise comparison
dds_hil_5d_B236_SC2_lrt <- DESeq(dds_hil_5d_B236_SC2, test="LRT", reduced=~1) # Use for multiple group comparison

resLFC_B236_SC2_hil_5d_lrt <- lfcShrink(dds_hil_5d_B236_SC2_lrt, coef="Virus_B236_vs_Saline", type="apeglm")
sorted_resLFC_B236_SC2_hil_5d_lrt <- resLFC_B236_SC2_hil_5d_lrt[order(resLFC_B236_SC2_hil_5d_lrt$padj),]
sorted_sig_resLFC_B236_SC2_hil_5d_lrt <- sorted_resLFC_B236_SC2_hil_5d_lrt[which(sorted_resLFC_B236_SC2_hil_5d_lrt$padj < 0.05),]
summary(sorted_sig_resLFC_B236_SC2_hil_5d_lrt)
hil_5d_deg_list_lrt <- rownames(sorted_sig_resLFC_B236_SC2_hil_5d_lrt)

######### 2.5.1 Calculate disco score
hil_5d_disco <- data.frame(pvalue_B236_hil=resLFC_B236_hil_5d_wald$pvalue, padj_B236_hil=resLFC_B236_hil_5d_wald$padj, log2FC_B236_hil=resLFC_B236_hil_5d_wald$log2FoldChange, 
                           pvalue_SC2_hil=resLFC_SC2_hil_5d_wald$pvalue, padj_SC2_hil=resLFC_SC2_hil_5d_wald$padj, log2FC_SC2_hil=resLFC_SC2_hil_5d_wald$log2FoldChange, row.names=rownames(resLFC_B236_hil_5d_wald))
hil_5d_disco[,c("category","B236_plot_category","SC2_plot_category","pool_plot_category","disco_score")] <- NA

######### 2.5.2 Assign gene expression category based on directions of fold change
for (i in 1:nrow(hil_5d_disco)) {
  if (is.na(hil_5d_disco[i,]$pvalue_B236_hil) || is.na(hil_5d_disco[i,]$padj_B236_hil) || rownames(hil_5d_disco[i,]) %in% saline_B236_hil_5d_inactive_genes || hil_5d_disco[i,]$padj_B236_hil>=0.05) {
    B236_hil_5d_direction <- "X"} else if (hil_5d_disco[i,]$log2FC_B236_hil<0) {B236_hil_5d_direction="D"} else {B236_hil_5d_direction="U"}
  if (is.na(hil_5d_disco[i,]$pvalue_SC2_hil) || is.na(hil_5d_disco[i,]$padj_SC2_hil) || rownames(hil_5d_disco[i,]) %in% saline_SC2_hil_5d_inactive_genes || hil_5d_disco[i,]$padj_SC2_hil>=0.05) {
    SC2_hil_5d_direction <- "X"} else if (hil_5d_disco[i,]$log2FC_SC2_hil<0) {SC2_hil_5d_direction="D"} else {SC2_hil_5d_direction="U"}
  hil_5d_disco[i,]$category <- paste(B236_hil_5d_direction, SC2_hil_5d_direction, sep="")
  if (!rownames(hil_5d_disco[i,]) %in% hil_5d_deg_list_lrt) {hil_5d_disco[i,]$category <- "XX"}
  if (hil_5d_disco[i,]$category %in% c("DD")) {hil_5d_disco[i,]$B236_plot_category <- "B236-Down/SC2-Down"}
  if (hil_5d_disco[i,]$category %in% c("UU")) {hil_5d_disco[i,]$B236_plot_category <- "B236-Up/SC2-Up"}
  if (hil_5d_disco[i,]$category %in% c("DU","DX")) {hil_5d_disco[i,]$B236_plot_category <- "B236-Down/SC2-notDown"}
  if (hil_5d_disco[i,]$category %in% c("UD","UX")) {hil_5d_disco[i,]$B236_plot_category <- "B236-Up/SC2-notUp"}
  if (is.na(hil_5d_disco[i,]$B236_plot_category)) {hil_5d_disco[i,]$B236_plot_category <- "OOI"}
  if (hil_5d_disco[i,]$category %in% c("DD")) {hil_5d_disco[i,]$SC2_plot_category <- "B236-Down/SC2-Down"}
  if (hil_5d_disco[i,]$category %in% c("UU")) {hil_5d_disco[i,]$SC2_plot_category <- "B236-Up/SC2-Up"}
  if (hil_5d_disco[i,]$category %in% c("UD","XD")) {hil_5d_disco[i,]$SC2_plot_category <- "B236-notDown/SC2-Down"}
  if (hil_5d_disco[i,]$category %in% c("DU","XU")) {hil_5d_disco[i,]$SC2_plot_category <- "B236-notUp/SC2-Up"}
  if (is.na(hil_5d_disco[i,]$SC2_plot_category)) {hil_5d_disco[i,]$SC2_plot_category <- "OOI"}
}

hil_5d_disco$disco_score <- abs(hil_5d_disco$log2FC_B236_hil * hil_5d_disco$log2FC_SC2_hil * (log10(hil_5d_disco$pvalue_B236_hil) + log10(hil_5d_disco$pvalue_SC2_hil)))
sorted_hil_5d_disco <- hil_5d_disco[order(hil_5d_disco$disco_score, decreasing=TRUE),]

sorted_hil_5d_disco$human_gene_name <- human_hamster_ortholog[match(rownames(sorted_hil_5d_disco), human_hamster_ortholog$symbol.hamster),]$symbol.human

sorted_top_B236_hil_5d_disco <- sorted_hil_5d_disco %>% filter(!is.na(human_gene_name))
sorted_top_B236_hil_5d_disco <- Reduce(rbind, by(sorted_top_B236_hil_5d_disco, sorted_top_B236_hil_5d_disco["B236_plot_category"], head, n=100))
sorted_top_SC2_hil_5d_disco <- sorted_hil_5d_disco %>% filter(!is.na(human_gene_name))
sorted_top_SC2_hil_5d_disco <- Reduce(rbind, by(sorted_top_SC2_hil_5d_disco, sorted_top_SC2_hil_5d_disco["SC2_plot_category"], head, n=100))

######################################## Periphery (2 d.p.i.) ########################################
########## 3. Analyze DEGs in hamster lung peripheries at 2 d.p.i.
ppr_filter_samples <- rbind(passed_samples %>% filter(Tissue=="Periphery"), passed_samples %>% filter(Tissue=="Periphery", Virus=="Saline") %>% mutate(Day="2dpi", Grouping="Saline_Periphery_2dpi"))
ppr_2d_filter_samples <- ppr_filter_samples %>% filter(Day=="2dpi")
dds_ppr_2d <- DESeqDataSetFromMatrix(countData=count.df_no_viral_seq_deseq[,ppr_2d_filter_samples$Sample], colData=ppr_2d_filter_samples, design=~Virus)
dds_ppr_2d$Virus <- relevel(dds_ppr_2d$Virus, "Saline")
dds_ppr_2d_wald <- DESeq(dds_ppr_2d) # Use for pairwise comparison
dds_ppr_2d_lrt <- DESeq(dds_ppr_2d, test="LRT", reduced=~1) # Use for multiple group comparison

####### 3.1 SC2 vs mock-infected
resLFC_B236_ppr_2d_wald <- lfcShrink(dds_ppr_2d_wald, coef="Virus_B236_vs_Saline", type="apeglm")
sorted_resLFC_B236_ppr_2d_wald <- resLFC_B236_ppr_2d_wald[order(resLFC_B236_ppr_2d_wald$padj),]
sorted_sig_resLFC_B236_ppr_2d_wald <- sorted_resLFC_B236_ppr_2d_wald[which(sorted_resLFC_B236_ppr_2d_wald$padj < 0.05),]
filtered_sorted_sig_resLFC_B236_ppr_2d_wald <- sorted_sig_resLFC_B236_ppr_2d_wald[which(!rownames(sorted_sig_resLFC_B236_ppr_2d_wald) %in% saline_B236_ppr_2d_inactive_genes),]
B236_ppr_2d_all_down_gene <- rownames(filtered_sorted_sig_resLFC_B236_ppr_2d_wald[which(filtered_sorted_sig_resLFC_B236_ppr_2d_wald$log2FoldChange < 0),])
B236_ppr_2d_all_up_gene <- rownames(filtered_sorted_sig_resLFC_B236_ppr_2d_wald[which(filtered_sorted_sig_resLFC_B236_ppr_2d_wald$log2FoldChange > 0),])
summary(filtered_sorted_sig_resLFC_B236_ppr_2d_wald)

####### 3.2 SC2dFCS vs mock-infected
resLFC_SC2_ppr_2d_wald <- lfcShrink(dds_ppr_2d_wald, coef="Virus_SC2_vs_Saline", type="apeglm")
sorted_resLFC_SC2_ppr_2d_wald <- resLFC_SC2_ppr_2d_wald[order(resLFC_SC2_ppr_2d_wald$padj),]
sorted_sig_resLFC_SC2_ppr_2d_wald <- sorted_resLFC_SC2_ppr_2d_wald[which(sorted_resLFC_SC2_ppr_2d_wald$padj < 0.05),]
filtered_sorted_sig_resLFC_SC2_ppr_2d_wald <- sorted_sig_resLFC_SC2_ppr_2d_wald[which(!rownames(sorted_sig_resLFC_SC2_ppr_2d_wald) %in% saline_SC2_ppr_2d_inactive_genes),]
SC2_ppr_2d_all_down_gene <- rownames(filtered_sorted_sig_resLFC_SC2_ppr_2d_wald[which(filtered_sorted_sig_resLFC_SC2_ppr_2d_wald$log2FoldChange < 0),])
SC2_ppr_2d_all_up_gene <- rownames(filtered_sorted_sig_resLFC_SC2_ppr_2d_wald[which(filtered_sorted_sig_resLFC_SC2_ppr_2d_wald$log2FoldChange > 0),])
summary(filtered_sorted_sig_resLFC_SC2_ppr_2d_wald)

####### 3.3 B236 vs mock-infected
resLFC_SC2dFCS_ppr_2d_wald <- lfcShrink(dds_ppr_2d_wald, coef="Virus_SC2dFCS_vs_Saline", type="apeglm")
sorted_resLFC_SC2dFCS_ppr_2d_wald <- resLFC_SC2dFCS_ppr_2d_wald[order(resLFC_SC2dFCS_ppr_2d_wald$padj),]
sorted_sig_resLFC_SC2dFCS_ppr_2d_wald <- sorted_resLFC_SC2dFCS_ppr_2d_wald[which(sorted_resLFC_SC2dFCS_ppr_2d_wald$padj < 0.05),]
filtered_sorted_sig_resLFC_SC2dFCS_ppr_2d_wald <- sorted_sig_resLFC_SC2dFCS_ppr_2d_wald[which(!rownames(sorted_sig_resLFC_SC2dFCS_ppr_2d_wald) %in% saline_SC2dFCS_ppr_2d_inactive_genes),]
SC2dFCS_ppr_2d_all_down_gene <- rownames(filtered_sorted_sig_resLFC_SC2dFCS_ppr_2d_wald[which(filtered_sorted_sig_resLFC_SC2dFCS_ppr_2d_wald$log2FoldChange < 0),])
SC2dFCS_ppr_2d_all_up_gene <- rownames(filtered_sorted_sig_resLFC_SC2dFCS_ppr_2d_wald[which(filtered_sorted_sig_resLFC_SC2dFCS_ppr_2d_wald$log2FoldChange > 0),])
summary(filtered_sorted_sig_resLFC_SC2dFCS_ppr_2d_wald)

####### 3.4 Summarize the number of DEG
ppr_2d_deg_summary.df <- data.frame(rep(c("B236","SC2","SC2dFCS"), 2), c(length(B236_ppr_2d_all_up_gene), length(SC2_ppr_2d_all_up_gene), length(SC2dFCS_ppr_2d_all_up_gene),
                                                                         length(B236_ppr_2d_all_down_gene), length(SC2_ppr_2d_all_down_gene), length(SC2dFCS_ppr_2d_all_down_gene)), c(rep(c("Upregulated"), 3), rep(c("Downregulated"), 3)))
colnames(ppr_2d_deg_summary.df) <- c("Virus","Genes","Direction")
ppr_2d_deg_summary.df$Grouping <- "Periphery_2dpi"

####### 3.5 Analyze the difference between B236 infection and SC2 or SC2dFCS infection in lung peripheries at 2 d.p.i.
ppr_2d_B236_SC2_filter_samples <- rbind(passed_samples %>% filter(Tissue=="Periphery", Day=="2dpi", !Virus=="SC2dFCS"), passed_samples %>% filter(Tissue=="Periphery", Virus=="Saline") %>% mutate(Day="2dpi", Grouping="Saline_periphery_2dpi"))
dds_ppr_2d_B236_SC2 <- DESeqDataSetFromMatrix(countData=count.df_no_viral_seq_deseq[,ppr_2d_B236_SC2_filter_samples$Sample], colData=ppr_2d_B236_SC2_filter_samples, design=~Virus)

dds_ppr_2d_B236_SC2$Virus <- relevel(dds_ppr_2d_B236_SC2$Virus, "Saline")
dds_ppr_2d_B236_SC2_wald <- DESeq(dds_ppr_2d_B236_SC2) # Use for pairwise comparison
dds_ppr_2d_B236_SC2_lrt <- DESeq(dds_ppr_2d_B236_SC2, test="LRT", reduced=~1) # Use for multiple group comparison

resLFC_B236_SC2_ppr_2d_lrt <- lfcShrink(dds_ppr_2d_B236_SC2_lrt, coef="Virus_B236_vs_Saline", type="apeglm")
sorted_resLFC_B236_SC2_ppr_2d_lrt <- resLFC_B236_SC2_ppr_2d_lrt[order(resLFC_B236_SC2_ppr_2d_lrt$padj),]
sorted_sig_resLFC_B236_SC2_ppr_2d_lrt <- sorted_resLFC_B236_SC2_ppr_2d_lrt[which(sorted_resLFC_B236_SC2_ppr_2d_lrt$padj < 0.05),]
summary(sorted_sig_resLFC_B236_SC2_ppr_2d_lrt)
ppr_2d_deg_list_lrt <- rownames(sorted_sig_resLFC_B236_SC2_ppr_2d_lrt)

######### 3.5.1 Calculate disco score
ppr_2d_disco <- data.frame(pvalue_B236_ppr=resLFC_B236_ppr_2d_wald$pvalue, padj_B236_ppr=resLFC_B236_ppr_2d_wald$padj, log2FC_B236_ppr=resLFC_B236_ppr_2d_wald$log2FoldChange, 
                           pvalue_SC2_ppr=resLFC_SC2_ppr_2d_wald$pvalue, padj_SC2_ppr=resLFC_SC2_ppr_2d_wald$padj, log2FC_SC2_ppr=resLFC_SC2_ppr_2d_wald$log2FoldChange, row.names=rownames(resLFC_B236_ppr_2d_wald))
ppr_2d_disco[,c("category","B236_plot_category","SC2_plot_category","pool_plot_category","disco_score")] <- NA

######### 3.5.2 Assign gene expression category based on directions of fold change
for (i in 1:nrow(ppr_2d_disco)) {
  if (is.na(ppr_2d_disco[i,]$pvalue_B236_ppr) || is.na(ppr_2d_disco[i,]$padj_B236_ppr) || rownames(ppr_2d_disco[i,]) %in% saline_B236_ppr_2d_inactive_genes || ppr_2d_disco[i,]$padj_B236_ppr>=0.05) {
    B236_ppr_2d_direction <- "X"} else if (ppr_2d_disco[i,]$log2FC_B236_ppr<0) {B236_ppr_2d_direction="D"} else {B236_ppr_2d_direction="U"}
  if (is.na(ppr_2d_disco[i,]$pvalue_SC2_ppr) || is.na(ppr_2d_disco[i,]$padj_SC2_ppr) || rownames(ppr_2d_disco[i,]) %in% saline_SC2_ppr_2d_inactive_genes || ppr_2d_disco[i,]$padj_SC2_ppr>=0.05) {
    SC2_ppr_2d_direction <- "X"} else if (ppr_2d_disco[i,]$log2FC_SC2_ppr<0) {SC2_ppr_2d_direction="D"} else {SC2_ppr_2d_direction="U"}
  ppr_2d_disco[i,]$category <- paste(B236_ppr_2d_direction, SC2_ppr_2d_direction, sep="")
  if (!rownames(ppr_2d_disco[i,]) %in% ppr_2d_deg_list_lrt) {ppr_2d_disco[i,]$category <- "XX"}
  if (ppr_2d_disco[i,]$category %in% c("DD")) {ppr_2d_disco[i,]$B236_plot_category <- "B236-Down/SC2-Down"}
  if (ppr_2d_disco[i,]$category %in% c("UU")) {ppr_2d_disco[i,]$B236_plot_category <- "B236-Up/SC2-Up"}
  if (ppr_2d_disco[i,]$category %in% c("DU","DX")) {ppr_2d_disco[i,]$B236_plot_category <- "B236-Down/SC2-notDown"}
  if (ppr_2d_disco[i,]$category %in% c("UD","UX")) {ppr_2d_disco[i,]$B236_plot_category <- "B236-Up/SC2-notUp"}
  if (is.na(ppr_2d_disco[i,]$B236_plot_category)) {ppr_2d_disco[i,]$B236_plot_category <- "OOI"}
  if (ppr_2d_disco[i,]$category %in% c("DD")) {ppr_2d_disco[i,]$SC2_plot_category <- "B236-Down/SC2-Down"}
  if (ppr_2d_disco[i,]$category %in% c("UU")) {ppr_2d_disco[i,]$SC2_plot_category <- "B236-Up/SC2-Up"}
  if (ppr_2d_disco[i,]$category %in% c("UD","XD")) {ppr_2d_disco[i,]$SC2_plot_category <- "B236-notDown/SC2-Down"}
  if (ppr_2d_disco[i,]$category %in% c("DU","XU")) {ppr_2d_disco[i,]$SC2_plot_category <- "B236-notUp/SC2-Up"}
  if (is.na(ppr_2d_disco[i,]$SC2_plot_category)) {ppr_2d_disco[i,]$SC2_plot_category <- "OOI"}
}

ppr_2d_disco$disco_score <- abs(ppr_2d_disco$log2FC_B236_ppr * ppr_2d_disco$log2FC_SC2_ppr * (log10(ppr_2d_disco$pvalue_B236_ppr) + log10(ppr_2d_disco$pvalue_SC2_ppr)))
sorted_ppr_2d_disco <- ppr_2d_disco[order(ppr_2d_disco$disco_score, decreasing=TRUE),]

sorted_ppr_2d_disco$human_gene_name <- human_hamster_ortholog[match(rownames(sorted_ppr_2d_disco), human_hamster_ortholog$symbol.hamster),]$symbol.human

sorted_top_B236_ppr_2d_disco <- sorted_ppr_2d_disco %>% filter(!is.na(human_gene_name))
sorted_top_B236_ppr_2d_disco <- Reduce(rbind, by(sorted_top_B236_ppr_2d_disco, sorted_top_B236_ppr_2d_disco["B236_plot_category"], head, n=100))
sorted_top_SC2_ppr_2d_disco <- sorted_ppr_2d_disco %>% filter(!is.na(human_gene_name))
sorted_top_SC2_ppr_2d_disco <- Reduce(rbind, by(sorted_top_SC2_ppr_2d_disco, sorted_top_SC2_ppr_2d_disco["SC2_plot_category"], head, n=100))

######################################## Periphery (5 d.p.i.) ########################################
########## 4. Analyze DEGs in hamster lung peripheries at 5 d.p.i.
ppr_5d_filter_samples <- ppr_filter_samples %>% filter(Day=="5dpi")
dds_ppr_5d <- DESeqDataSetFromMatrix(countData=count.df_no_viral_seq_deseq[,ppr_5d_filter_samples$Sample], colData=ppr_5d_filter_samples, design=~Virus)
dds_ppr_5d$Virus <- relevel(dds_ppr_5d$Virus, "Saline")
dds_ppr_5d_wald <- DESeq(dds_ppr_5d) # Use for pairwise comparison
dds_ppr_5d_lrt <- DESeq(dds_ppr_5d, test="LRT", reduced=~1) # Use for multiple group comparison

####### 4.1 SC2 vs mock-infected
resLFC_B236_ppr_5d_wald <- lfcShrink(dds_ppr_5d_wald, coef="Virus_B236_vs_Saline", type="apeglm")
sorted_resLFC_B236_ppr_5d_wald <- resLFC_B236_ppr_5d_wald[order(resLFC_B236_ppr_5d_wald$padj),]
sorted_sig_resLFC_B236_ppr_5d_wald <- sorted_resLFC_B236_ppr_5d_wald[which(sorted_resLFC_B236_ppr_5d_wald$padj < 0.05),]
filtered_sorted_sig_resLFC_B236_ppr_5d_wald <- sorted_sig_resLFC_B236_ppr_5d_wald[which(!rownames(sorted_sig_resLFC_B236_ppr_5d_wald) %in% saline_B236_ppr_5d_inactive_genes),]
B236_ppr_5d_all_down_gene <- rownames(filtered_sorted_sig_resLFC_B236_ppr_5d_wald[which(filtered_sorted_sig_resLFC_B236_ppr_5d_wald$log2FoldChange < 0),])
B236_ppr_5d_all_up_gene <- rownames(filtered_sorted_sig_resLFC_B236_ppr_5d_wald[which(filtered_sorted_sig_resLFC_B236_ppr_5d_wald$log2FoldChange > 0),])
summary(filtered_sorted_sig_resLFC_B236_ppr_5d_wald)

####### 4.2 SC2dFCS vs mock-infected
resLFC_SC2_ppr_5d_wald <- lfcShrink(dds_ppr_5d_wald, coef="Virus_SC2_vs_Saline", type="apeglm")
sorted_resLFC_SC2_ppr_5d_wald <- resLFC_SC2_ppr_5d_wald[order(resLFC_SC2_ppr_5d_wald$padj),]
sorted_sig_resLFC_SC2_ppr_5d_wald <- sorted_resLFC_SC2_ppr_5d_wald[which(sorted_resLFC_SC2_ppr_5d_wald$padj < 0.05),]
filtered_sorted_sig_resLFC_SC2_ppr_5d_wald <- sorted_sig_resLFC_SC2_ppr_5d_wald[which(!rownames(sorted_sig_resLFC_SC2_ppr_5d_wald) %in% saline_SC2_ppr_5d_inactive_genes),]
SC2_ppr_5d_all_down_gene <- rownames(filtered_sorted_sig_resLFC_SC2_ppr_5d_wald[which(filtered_sorted_sig_resLFC_SC2_ppr_5d_wald$log2FoldChange < 0),])
SC2_ppr_5d_all_up_gene <- rownames(filtered_sorted_sig_resLFC_SC2_ppr_5d_wald[which(filtered_sorted_sig_resLFC_SC2_ppr_5d_wald$log2FoldChange > 0),])
summary(filtered_sorted_sig_resLFC_SC2_ppr_5d_wald)

####### 4.3 B236 vs mock-infected
resLFC_SC2dFCS_ppr_5d_wald <- lfcShrink(dds_ppr_5d_wald, coef="Virus_SC2dFCS_vs_Saline", type="apeglm")
sorted_resLFC_SC2dFCS_ppr_5d_wald <- resLFC_SC2dFCS_ppr_5d_wald[order(resLFC_SC2dFCS_ppr_5d_wald$padj),]
sorted_sig_resLFC_SC2dFCS_ppr_5d_wald <- sorted_resLFC_SC2dFCS_ppr_5d_wald[which(sorted_resLFC_SC2dFCS_ppr_5d_wald$padj < 0.05),]
filtered_sorted_sig_resLFC_SC2dFCS_ppr_5d_wald <- sorted_sig_resLFC_SC2dFCS_ppr_5d_wald[which(!rownames(sorted_sig_resLFC_SC2dFCS_ppr_5d_wald) %in% saline_SC2dFCS_ppr_5d_inactive_genes),]
SC2dFCS_ppr_5d_all_down_gene <- rownames(filtered_sorted_sig_resLFC_SC2dFCS_ppr_5d_wald[which(filtered_sorted_sig_resLFC_SC2dFCS_ppr_5d_wald$log2FoldChange < 0),])
SC2dFCS_ppr_5d_all_up_gene <- rownames(filtered_sorted_sig_resLFC_SC2dFCS_ppr_5d_wald[which(filtered_sorted_sig_resLFC_SC2dFCS_ppr_5d_wald$log2FoldChange > 0),])
summary(filtered_sorted_sig_resLFC_SC2dFCS_ppr_5d_wald)

####### 4.4 Summarize the number of DEG
ppr_5d_deg_summary.df <- data.frame(rep(c("B236","SC2","SC2dFCS"), 2), c(length(B236_ppr_5d_all_up_gene), length(SC2_ppr_5d_all_up_gene), length(SC2dFCS_ppr_5d_all_up_gene), 
                                                                         length(B236_ppr_5d_all_down_gene), length(SC2_ppr_5d_all_down_gene), length(SC2dFCS_ppr_5d_all_down_gene)), c(rep(c("Upregulated"), 3), rep(c("Downregulated"), 3)))
colnames(ppr_5d_deg_summary.df) <- c("Virus","Genes","Direction")
ppr_5d_deg_summary.df$Grouping <- "Periphery_5dpi"

####### 4.5 Analyze the difference between B236 infection and SC2 or SC2dFCS infection in lung peripheries at 5 d.p.i.
ppr_5d_B236_SC2_filter_samples <- rbind(passed_samples %>% filter(Tissue=="Periphery", Day=="5dpi", !Virus=="SC2dFCS"), passed_samples %>% filter(Tissue=="Periphery", Virus=="Saline") %>% mutate(Day="5dpi", Grouping="Saline_periphery_5dpi"))
dds_ppr_5d_B236_SC2 <- DESeqDataSetFromMatrix(countData=count.df_no_viral_seq_deseq[,ppr_5d_B236_SC2_filter_samples$Sample], colData=ppr_5d_B236_SC2_filter_samples, design=~Virus)

dds_ppr_5d_B236_SC2$Virus <- relevel(dds_ppr_5d_B236_SC2$Virus, "Saline")
dds_ppr_5d_B236_SC2_wald <- DESeq(dds_ppr_5d_B236_SC2) # Use for pairwise comparison
dds_ppr_5d_B236_SC2_lrt <- DESeq(dds_ppr_5d_B236_SC2, test="LRT", reduced=~1) # Use for multiple group comparison

resLFC_B236_SC2_ppr_5d_lrt <- lfcShrink(dds_ppr_5d_B236_SC2_lrt, coef="Virus_B236_vs_Saline", type="apeglm")
sorted_resLFC_B236_SC2_ppr_5d_lrt <- resLFC_B236_SC2_ppr_5d_lrt[order(resLFC_B236_SC2_ppr_5d_lrt$padj),]
sorted_sig_resLFC_B236_SC2_ppr_5d_lrt <- sorted_resLFC_B236_SC2_ppr_5d_lrt[which(sorted_resLFC_B236_SC2_ppr_5d_lrt$padj < 0.05),]
summary(sorted_sig_resLFC_B236_SC2_ppr_5d_lrt)
ppr_5d_deg_list_lrt <- rownames(sorted_sig_resLFC_B236_SC2_ppr_5d_lrt)

######### 4.5.1 Calculate disco score
ppr_5d_disco <- data.frame(pvalue_B236_ppr=resLFC_B236_ppr_5d_wald$pvalue, padj_B236_ppr=resLFC_B236_ppr_5d_wald$padj, log2FC_B236_ppr=resLFC_B236_ppr_5d_wald$log2FoldChange, pvalue_SC2_ppr=resLFC_SC2_ppr_5d_wald$pvalue, padj_SC2_ppr=resLFC_SC2_ppr_5d_wald$padj, log2FC_SC2_ppr=resLFC_SC2_ppr_5d_wald$log2FoldChange, row.names=rownames(resLFC_B236_ppr_5d_wald))
ppr_5d_disco[,c("category","B236_plot_category","SC2_plot_category","pool_plot_category","disco_score")] <- NA

######### 4.5.2 Assign gene expression category based on directions of fold change
for (i in 1:nrow(ppr_5d_disco)) {
  if (is.na(ppr_5d_disco[i,]$pvalue_B236_ppr) || is.na(ppr_5d_disco[i,]$padj_B236_ppr) || rownames(ppr_5d_disco[i,]) %in% saline_B236_ppr_5d_inactive_genes || ppr_5d_disco[i,]$padj_B236_ppr>=0.05) {
    B236_ppr_5d_direction <- "X"} else if (ppr_5d_disco[i,]$log2FC_B236_ppr<0) {B236_ppr_5d_direction="D"} else {B236_ppr_5d_direction="U"}
  if (is.na(ppr_5d_disco[i,]$pvalue_SC2_ppr) || is.na(ppr_5d_disco[i,]$padj_SC2_ppr) || rownames(ppr_5d_disco[i,]) %in% saline_SC2_ppr_5d_inactive_genes || ppr_5d_disco[i,]$padj_SC2_ppr>=0.05) {
    SC2_ppr_5d_direction <- "X"} else if (ppr_5d_disco[i,]$log2FC_SC2_ppr<0) {SC2_ppr_5d_direction="D"} else {SC2_ppr_5d_direction="U"}
  ppr_5d_disco[i,]$category <- paste(B236_ppr_5d_direction, SC2_ppr_5d_direction, sep="")
  if (!rownames(ppr_5d_disco[i,]) %in% ppr_5d_deg_list_lrt) {ppr_5d_disco[i,]$category <- "XX"}
  if (ppr_5d_disco[i,]$category %in% c("DD")) {ppr_5d_disco[i,]$B236_plot_category <- "B236-Down/SC2-Down"}
  if (ppr_5d_disco[i,]$category %in% c("UU")) {ppr_5d_disco[i,]$B236_plot_category <- "B236-Up/SC2-Up"}
  if (ppr_5d_disco[i,]$category %in% c("DU","DX")) {ppr_5d_disco[i,]$B236_plot_category <- "B236-Down/SC2-notDown"}
  if (ppr_5d_disco[i,]$category %in% c("UD","UX")) {ppr_5d_disco[i,]$B236_plot_category <- "B236-Up/SC2-notUp"}
  if (is.na(ppr_5d_disco[i,]$B236_plot_category)) {ppr_5d_disco[i,]$B236_plot_category <- "OOI"}
  if (ppr_5d_disco[i,]$category %in% c("DD")) {ppr_5d_disco[i,]$SC2_plot_category <- "B236-Down/SC2-Down"}
  if (ppr_5d_disco[i,]$category %in% c("UU")) {ppr_5d_disco[i,]$SC2_plot_category <- "B236-Up/SC2-Up"}
  if (ppr_5d_disco[i,]$category %in% c("UD","XD")) {ppr_5d_disco[i,]$SC2_plot_category <- "B236-notDown/SC2-Down"}
  if (ppr_5d_disco[i,]$category %in% c("DU","XU")) {ppr_5d_disco[i,]$SC2_plot_category <- "B236-notUp/SC2-Up"}
  if (is.na(ppr_5d_disco[i,]$SC2_plot_category)) {ppr_5d_disco[i,]$SC2_plot_category <- "OOI"}
}

ppr_5d_disco$disco_score <- abs(ppr_5d_disco$log2FC_B236_ppr * ppr_5d_disco$log2FC_SC2_ppr * (log10(ppr_5d_disco$pvalue_B236_ppr) + log10(ppr_5d_disco$pvalue_SC2_ppr)))
sorted_ppr_5d_disco <- ppr_5d_disco[order(ppr_5d_disco$disco_score, decreasing=TRUE),]

sorted_ppr_5d_disco$human_gene_name <- human_hamster_ortholog[match(rownames(sorted_ppr_5d_disco), human_hamster_ortholog$symbol.hamster),]$symbol.human

sorted_top_B236_ppr_5d_disco <- sorted_ppr_5d_disco %>% filter(!is.na(human_gene_name))
sorted_top_B236_ppr_5d_disco <- Reduce(rbind, by(sorted_top_B236_ppr_5d_disco, sorted_top_B236_ppr_5d_disco["B236_plot_category"], head, n=100))
sorted_top_SC2_ppr_5d_disco <- sorted_ppr_5d_disco %>% filter(!is.na(human_gene_name))
sorted_top_SC2_ppr_5d_disco <- Reduce(rbind, by(sorted_top_SC2_ppr_5d_disco, sorted_top_SC2_ppr_5d_disco["SC2_plot_category"], head, n=100))

######### Plot numbers of DEGs in hamster lung hila and peripheries (Figure 3G)
pool_deg_summary.df <- rbind(hil_2d_deg_summary.df, hil_5d_deg_summary.df, ppr_2d_deg_summary.df, ppr_5d_deg_summary.df)
pool_deg_summary.df <- pool_deg_summary.df %>% mutate(Tissue = str_split(Grouping, "_", simplify=T)[,1], Day = str_split(Grouping, "_", simplify=T)[,2])
pool_deg_summary.df$Virus <- factor(pool_deg_summary.df$Virus, levels=c("SC2","SC2dFCS","B236"))
pool_deg_summary.df$Direction <- factor(pool_deg_summary.df$Direction, levels=c("Upregulated","Downregulated"))

# pdf("pool_deg_num_line_hamster_lung.pdf", width=12, height=6)
pool_deg_line_plot <- ggplot(data=pool_deg_summary.df, aes(x=Day, y=Genes, color=Virus, group=Virus)) 
pool_deg_line_plot <- pool_deg_line_plot + geom_line() 
pool_deg_line_plot <- pool_deg_line_plot + geom_point() 
pool_deg_line_plot <- pool_deg_line_plot + theme_classic() 
pool_deg_line_plot <- pool_deg_line_plot + theme(text=element_text(size=24), strip.text=element_text(size=20), strip.placement="outside") 
pool_deg_line_plot <- pool_deg_line_plot + scale_y_continuous(limits=c(0,6000), breaks=c(0,1000,2000,3000,4000,5000,6000)) 
pool_deg_line_plot <- pool_deg_line_plot + facet_wrap(Tissue~Direction, nrow=1)
pool_deg_line_plot
# dev.off()

####### Create the table of DEG (Supplementary Table 3)
hamster_deg_list <- list(rownames(filtered_sorted_sig_resLFC_SC2_hil_2d_wald), rownames(filtered_sorted_sig_resLFC_SC2dFCS_hil_2d_wald), rownames(filtered_sorted_sig_resLFC_B236_hil_2d_wald), 
                         rownames(filtered_sorted_sig_resLFC_SC2_hil_5d_wald), rownames(filtered_sorted_sig_resLFC_SC2dFCS_hil_5d_wald), rownames(filtered_sorted_sig_resLFC_B236_hil_5d_wald), 
                         rownames(filtered_sorted_sig_resLFC_SC2_ppr_2d_wald), rownames(filtered_sorted_sig_resLFC_SC2dFCS_ppr_2d_wald), rownames(filtered_sorted_sig_resLFC_B236_ppr_2d_wald), 
                         rownames(filtered_sorted_sig_resLFC_SC2_ppr_5d_wald), rownames(filtered_sorted_sig_resLFC_SC2dFCS_ppr_5d_wald), rownames(filtered_sorted_sig_resLFC_B236_ppr_5d_wald))
mx_hamster_deg_list <- max(lengths(hamster_deg_list))
hamster_deg.df <- as.data.frame(as.data.table(do.call(cbind, lapply(hamster_deg_list, `length<-`, mx_hamster_deg_list))))
colnames(hamster_deg.df) <- c("Hilum_2dpi_SC2", "Hilum_2dpi_SC2dFCS", "Hilum_2dpi_B236", "Hilum_5dpi_SC2", "Hilum_5dpi_SC2dFCS", "Hilum_5dpi_B236", 
                              "Periphery_2dpi_SC2", "Periphery_2dpi_SC2dFCS", "Periphery_2dpi_B236", "Periphery_5dpi_SC2", "Periphery_5dpi_SC2dFCS", "Periphery_5dpi_B236")
# write.table(hamster_deg.df, "table_s3_hamster_lung_deg.tsv", col.names=T, row.names=F, sep="\t", quote=F)



########## 5.Perform enrichment analyses
####### 5.1 Create lists of the top 100 genes in the infected hamster lungs with highest disco scores
highest_hil_2d_disco_dd_gene <- (sorted_top_B236_hil_2d_disco %>% filter(B236_plot_category=="B236-Down/SC2-Down"))$human_gene_name
highest_hil_2d_disco_uu_gene <- (sorted_top_B236_hil_2d_disco %>% filter(B236_plot_category=="B236-Up/SC2-Up"))$human_gene_name
highest_hil_2d_disco_du_dx_gene <- (sorted_top_B236_hil_2d_disco %>% filter(B236_plot_category=="B236-Down/SC2-notDown"))$human_gene_name
highest_hil_2d_disco_ud_ux_gene <- (sorted_top_B236_hil_2d_disco %>% filter(B236_plot_category=="B236-Up/SC2-notUp"))$human_gene_name
highest_hil_2d_disco_ud_xd_gene <- (sorted_top_SC2_hil_2d_disco %>% filter(SC2_plot_category=="B236-notDown/SC2-Down"))$human_gene_name
highest_hil_2d_disco_du_xu_gene <- (sorted_top_SC2_hil_2d_disco %>% filter(SC2_plot_category=="B236-notUp/SC2-Up"))$human_gene_name
highest_hil_5d_disco_dd_gene <- (sorted_top_B236_hil_5d_disco %>% filter(B236_plot_category=="B236-Down/SC2-Down"))$human_gene_name
highest_hil_5d_disco_uu_gene <- (sorted_top_B236_hil_5d_disco %>% filter(B236_plot_category=="B236-Up/SC2-Up"))$human_gene_name
highest_hil_5d_disco_du_dx_gene <- (sorted_top_B236_hil_5d_disco %>% filter(B236_plot_category=="B236-Down/SC2-notDown"))$human_gene_name
highest_hil_5d_disco_ud_ux_gene <- (sorted_top_B236_hil_5d_disco %>% filter(B236_plot_category=="B236-Up/SC2-notUp"))$human_gene_name
highest_hil_5d_disco_ud_xd_gene <- (sorted_top_SC2_hil_5d_disco %>% filter(SC2_plot_category=="B236-notDown/SC2-Down"))$human_gene_name
highest_hil_5d_disco_du_xu_gene <- (sorted_top_SC2_hil_5d_disco %>% filter(SC2_plot_category=="B236-notUp/SC2-Up"))$human_gene_name
highest_ppr_2d_disco_dd_gene <- (sorted_top_B236_ppr_2d_disco %>% filter(B236_plot_category=="B236-Down/SC2-Down"))$human_gene_name
highest_ppr_2d_disco_uu_gene <- (sorted_top_B236_ppr_2d_disco %>% filter(B236_plot_category=="B236-Up/SC2-Up"))$human_gene_name
highest_ppr_2d_disco_du_dx_gene <- (sorted_top_B236_ppr_2d_disco %>% filter(B236_plot_category=="B236-Down/SC2-notDown"))$human_gene_name
highest_ppr_2d_disco_ud_ux_gene <- (sorted_top_B236_ppr_2d_disco %>% filter(B236_plot_category=="B236-Up/SC2-notUp"))$human_gene_name
highest_ppr_2d_disco_ud_xd_gene <- (sorted_top_SC2_ppr_2d_disco %>% filter(SC2_plot_category=="B236-notDown/SC2-Down"))$human_gene_name
highest_ppr_2d_disco_du_xu_gene <- (sorted_top_SC2_ppr_2d_disco %>% filter(SC2_plot_category=="B236-notUp/SC2-Up"))$human_gene_name
highest_ppr_5d_disco_dd_gene <- (sorted_top_B236_ppr_5d_disco %>% filter(B236_plot_category=="B236-Down/SC2-Down"))$human_gene_name
highest_ppr_5d_disco_uu_gene <- (sorted_top_B236_ppr_5d_disco %>% filter(B236_plot_category=="B236-Up/SC2-Up"))$human_gene_name
highest_ppr_5d_disco_du_dx_gene <- (sorted_top_B236_ppr_5d_disco %>% filter(B236_plot_category=="B236-Down/SC2-notDown"))$human_gene_name
highest_ppr_5d_disco_ud_ux_gene <- (sorted_top_B236_ppr_5d_disco %>% filter(B236_plot_category=="B236-Up/SC2-notUp"))$human_gene_name
highest_ppr_5d_disco_ud_xd_gene <- (sorted_top_SC2_ppr_5d_disco %>% filter(SC2_plot_category=="B236-notDown/SC2-Down"))$human_gene_name
highest_ppr_5d_disco_du_xu_gene <- (sorted_top_SC2_ppr_5d_disco %>% filter(SC2_plot_category=="B236-notUp/SC2-Up"))$human_gene_name

mDIR <- "gmt_files/"
gmt_files <- paste0(mDIR, list.files(mDIR, pattern = ".gmt"))

EnrichmentAnalysis_Pool <- function(gmt_file, genelist_file, DIR) {
  # gmt_file: an item in gmt_files
  # genelist_file: gene list of interest
  # DIR: directory to export csv and pdf files
  
  # import gmt file
  enrichr_gmt <- read.gmt(gmt_file)
  
  # get gmt file name
  gmt_file_strip <- str_split(gmt_file, "/", simplify = T)
  gmt_file_name <- str_split(gmt_file_strip[length(gmt_file_strip)], "\\.", simplify = T)[1]
  title_name <- paste(str_split(gmt_file_name, "_", simplify = T), collapse=' ')
  print(paste0("calc enrichment with ", gmt_file_name, " ..."))
  
  # calc. enrichment score by clusterProfiler
  # enricher() is a function for hypergeometric test
  # see http://yulab-smu.top/clusterProfiler-book/chapter3.html
  CPfile <- compareCluster(geneCluster=genelist_file, fun='enricher', minGSSize=1, maxGSSize=3000, TERM2GENE=enrichr_gmt, pvalueCutoff=0.05)
  if (grepl("GO_", gmt_file)==TRUE) {
    CPfile@compareClusterResult$Description <- str_split(CPfile@compareClusterResult$Description, " \\(GO:", simplify=T)[,1]
  }
  if (grepl("MGI_Mammalian_Phenotype_Level_4_2019", gmt_file)==TRUE) {
    CPfile@compareClusterResult$Description <- unlist(lapply(strsplit(CPfile@compareClusterResult$Description, " "), FUN=function(x){paste(x[2:length(x)], collapse=" ")}))
  }
  if (grepl("MGI_Mammalian_Phenotype_Level_4_2021", gmt_file)==TRUE) {
    CPfile@compareClusterResult$Description <- str_split(CPfile@compareClusterResult$Description, " MP:", simplify=T)[,1]
  }
  CPfile@compareClusterResult$group <- NA
  CPfile@compareClusterResult$Category <- NA
  for (i in 1:nrow(CPfile@compareClusterResult)) {
    if (substr(CPfile@compareClusterResult$Cluster[i],1,2)=="H2") {CPfile@compareClusterResult$group[i] <- "Hilum 2 d.p.i."}
    if (substr(CPfile@compareClusterResult$Cluster[i],1,2)=="H5") {CPfile@compareClusterResult$group[i] <- "Hilum 5 d.p.i."}
    if (substr(CPfile@compareClusterResult$Cluster[i],1,2)=="P2") {CPfile@compareClusterResult$group[i] <- "Periphery 2 d.p.i."}
    if (substr(CPfile@compareClusterResult$Cluster[i],1,2)=="P5") {CPfile@compareClusterResult$group[i] <- "Periphery 5 d.p.i."}
  }
  CPfile@compareClusterResult$dummy_category <- str_split(CPfile@compareClusterResult$Cluster, ":", simplify=T)[,2]
  CPfile@compareClusterResult$Category <- paste(str_split(CPfile@compareClusterResult$Cluster, ":", simplify=T)[,2], "\n", "(", as.character(str_split(CPfile@compareClusterResult$GeneRatio, "/", simplify=T)[,2]), ")", sep="")
  CPfile@compareClusterResult$group <- factor(CPfile@compareClusterResult$group, levels=c("Hilum","Hilum 2 d.p.i.","Hilum 5 d.p.i.","Periphery","Periphery 2 d.p.i.","Periphery 5 d.p.i."))
  desired_level_factor <- intersect(c("B236-Up/SC2-Up","B236-Down/SC2-Down","B236-Up/SC2-notUp","B236-Down/SC2-notDown","B236-notUp/SC2-Up","B236-notDown/SC2-Down"), unique(CPfile@compareClusterResult$dummy_category))
  plot_level_factor <- unique(CPfile@compareClusterResult$Category)
  CPfile@compareClusterResult$Category <- factor(CPfile@compareClusterResult$Category, levels=plot_level_factor)
  
  # export graphs showCategory = 10, 20, 30 and 50
  # Set label_format for the maximum number characters to wrap text
  for (j in c(10)) {
    w <- myWidthAlgorithm(CPfile, j)
    g <- dotplot(CPfile, x="Category", showCategory=j, font.size=14, label_format=100) 
    g <- g + ggtitle(title_name) 
    g <- g + theme(text=element_text(size=14), plot.title=element_text(face="bold", size=14), legend.title=element_text(size=14), legend.text=element_text(size=14), legend.justification="top", 
                   axis.text.x=element_text(angle=90, vjust=0.5, hjust=1), strip.text.x=element_text(size=14, face="bold"), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) 
    g <- g + guides(size=guide_legend(title="Gene ratio", order=1), colour=guide_colourbar(title=bquote(Adjusted~italic(p)), order=2, reverse=TRUE)) 
    g <- g + facet_grid(~group, space="free", scales="free")
    ggsave(plot=g, filename=paste0(DIR, "/", gmt_file_name, "_show_", j, ".pdf"), width=w, h=6)
  }
  
  # export table
  write_csv(as.data.frame(CPfile), paste0(DIR, "/", gmt_file_name,".csv"))
  return(g)
}

myWidthAlgorithm <- function(CPfile, n) {
  df <- CPfile@compareClusterResult
  df_show <- df %>% group_by(Cluster) %>% top_n(-qvalue, n = n)
  maxnchar <- max(nchar(df_show$ID))
  nclu <- df$Cluster %>% unique() %>% length()
  return(nclu*1.5 + maxnchar/10) 
}

####### 5.2 Perform functional enrichment analysis with top 100 genes in the infected hamster lungs with highest disco scores (Supplementary Figure 3)
compare_samples <- list(
  "H2:B236-Up/SC2-Up" = highest_hil_2d_disco_uu_gene,
  "H2:B236-Down/SC2-Down" = highest_hil_2d_disco_dd_gene,
  "H2:B236-Up/SC2-notUp" = highest_hil_2d_disco_ud_ux_gene,
  "H2:B236-Down/SC2-notDown" = highest_hil_2d_disco_du_dx_gene,
  "H2:B236-notUp/SC2-Up" = highest_hil_2d_disco_du_xu_gene,
  "H2:B236-notDown/SC2-Down" = highest_hil_2d_disco_ud_xd_gene,
  "H5:B236-Up/SC2-Up" = highest_hil_5d_disco_uu_gene,
  "H5:B236-Down/SC2-Down" = highest_hil_5d_disco_dd_gene,
  "H5:B236-Up/SC2-notUp" = highest_hil_5d_disco_ud_ux_gene,
  "H5:B236-Down/SC2-notDown" = highest_hil_5d_disco_du_dx_gene,
  "H5:B236-notUp/SC2-Up" = highest_hil_5d_disco_du_xu_gene,
  "H5:B236-notDown/SC2-Down" = highest_hil_5d_disco_ud_xd_gene,
  "P2:B236-Up/SC2-Up" = highest_ppr_2d_disco_uu_gene,
  "P2:B236-Down/SC2-Down" = highest_ppr_2d_disco_dd_gene,
  "P2:B236-Up/SC2-notUp" = highest_ppr_2d_disco_ud_ux_gene,
  "P2:B236-Down/SC2-notDown" = highest_ppr_2d_disco_du_dx_gene,
  "P2:B236-notUp/SC2-Up" = highest_ppr_2d_disco_du_xu_gene,
  "P2:B236-notDown/SC2-Down" = highest_ppr_2d_disco_ud_xd_gene,
  "P5:B236-Up/SC2-Up" = highest_ppr_5d_disco_uu_gene,
  "P5:B236-Down/SC2-Down" = highest_ppr_5d_disco_dd_gene,
  "P5:B236-Up/SC2-notUp" = highest_ppr_5d_disco_ud_ux_gene,
  "P5:B236-Down/SC2-notDown" = highest_ppr_5d_disco_du_dx_gene,
  "P5:B236-notUp/SC2-Up" = highest_ppr_5d_disco_du_xu_gene,
  "P5:B236-notDown/SC2-Down" = highest_ppr_5d_disco_ud_xd_gene
)

DIR <- "230424_compareCluster_hilum_periphery_time_sep"
dir.create(DIR)

myplots <- list()
start.time <- proc.time() # start time
for (i in 1:length(gmt_files)){
  skip_to_next <- FALSE
  gmt_file <- gmt_files[i]
  tryCatch(myplots[[i]] <- EnrichmentAnalysis_Pool(gmt_file, compare_samples, DIR), error = function(e) {skip_to_next <<- TRUE})
  if(skip_to_next) { next } 
}
end.time <- proc.time() # finish time
(end.time-start.time) # passed time

# pdf("enrichment_top_100_hamster_lung.pdf", width=20, height=50)
enrich_plot.hamster_lung <- plot_grid(myplots[[3]], myplots[[11]], labels=c("C","D"), label_size=32, ncol=1, align="v", axis="lr", rel_heights=c(0.50,0.50))
enrich_plot.hamster_lung
# dev.off()
