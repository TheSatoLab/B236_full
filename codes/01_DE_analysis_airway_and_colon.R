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
dir <- "/Users/chainorato/Desktop/BANAL_project/230105_Airway_and_Colon/"
setwd(dir)

######################################## Analyze read counts with total viral reads ########################################
########## Read count table and record grouping information
count.df <- read.table(file.path(dir, "B236_airway_lung_colon_counts.txt"), header=TRUE)
rownames(count.df) <- count.df$Geneid

samples <- read.delim(file.path(dir, "B236_airway_lung_colon_samples.txt"), header=TRUE)
samples <- samples %>% mutate(Virus = ifelse(Virus=="Uninfected", "Mock_infected",
                                      ifelse(Virus=="WK-521", "SC2",
                                      ifelse(Virus=="WK-521 delFCS", "SC2dFCS",
                                      ifelse(Virus=="BANAL-236", "B236", NA)))))
samples <- samples %>% mutate(cell.type = ifelse(cell.type=="HiLung_Airway", "Airway",
                                          ifelse(cell.type=="HiLung_Lung", "Lung",
                                          ifelse(cell.type=="Colon organoid", "Colon_organoid", NA))))
samples <- samples %>% mutate(grouping = paste(Virus, cell.type, sep="."))
samples <- samples %>% mutate(Name = ave(grouping, grouping, FUN = function(i) paste0(i, '_', seq_along(i))))
samples <- samples %>% mutate(fastq_f_dir = ifelse(cell.type=="Airway", paste("X221021_A01297_0339_AHN3TFDRX2.bam.", fastq_f, ".", fastq_f, "_Aligned.out.bam", sep=""),
                                            ifelse(cell.type=="Lung", paste("X221124_A00691_0609_BHLKHCDRX2.bam.", fastq_f, ".", fastq_f, "_Aligned.out.bam", sep=""),
                                            ifelse(cell.type=="Colon_organoid", paste("X221202_A01297_0354_AHNJ7YDRX2.bam.", str_replace_all(fastq_f, "-", "."), ".", str_replace_all(fastq_f, "-", "."), "_Aligned.out.bam", sep=""), NA))))
sample_dict <- unique(samples[,c("fastq_f_dir","Name")])

colnames(count.df) <- c(colnames(count.df)[1:6], sample_dict$Name[match(colnames(count.df[7:ncol(count.df)]), sample_dict$fastq_f_dir)])

count.df_deseq <- count.df
count.df_deseq[1:6] <- NULL
count.df_deseq <- count.df_deseq[samples$Name]

##### Since the genome length of B236 (29,844 bp) and that of SC2 (29,903 bp) are similar, it's unnecessary to normalize counts with virus genome length. #####

total_count_per_each_exp <- colSums(count.df_deseq)
total_viral_count.df <- count.df_deseq[(nrow(count.df_deseq)-1):nrow(count.df_deseq),]
norm_viral_count.mx <- sweep(total_viral_count.df*1000000, 2, total_count_per_each_exp, FUN="/")
norm_viral_count.df <- as.data.frame.table(as.matrix(norm_viral_count.mx))
colnames(norm_viral_count.df) <- c("Genome","grouping","norm_count")

norm_viral_count.df$Genome <- str_replace(norm_viral_count.df$Genome, "WK-521", "SC2")
norm_viral_count.df$Genome <- str_replace(norm_viral_count.df$Genome, "MZ937003.2", "B236")
norm_viral_count.df$grouping <- substr(as.vector(norm_viral_count.df$grouping), 1, nchar(as.vector(norm_viral_count.df$grouping))-2)
norm_viral_count.df$Virus <- str_split(norm_viral_count.df$grouping, "\\.", simplify=T)[,1]
norm_viral_count.df$cell_type <- str_split(norm_viral_count.df$grouping, "\\.", simplify=T)[,2]
norm_viral_count.df$log_norm_count <- log10(norm_viral_count.df$norm_count)

norm_viral_count.df <- norm_viral_count.df %>% filter(Virus!="Mock_infected") %>% filter((Genome=="SC2" & Virus=="SC2") | (Genome=="SC2" & Virus=="SC2dFCS") | (Genome=="B236" & Virus=="B236"))
norm_viral_count.df$Virus <- factor(norm_viral_count.df$Virus, levels=c("SC2","SC2dFCS","B236"))
norm_viral_count.df$Genome <- factor(norm_viral_count.df$Genome, levels=c("SC2","B236"))



########## Calculate the normalized counts of SC2 and B236 in airway epithelial cells (Figure 1F)
norm_viral_count.df.airway <- norm_viral_count.df %>% filter(cell_type=="Airway") 
norm_viral_count.df.airway %>% filter(log_norm_count!=-Inf) %>% group_by(Virus) %>% shapiro_test(log_norm_count)
norm_viral_count.df.airway %>% filter(log_norm_count!=-Inf) %>% levene_test(log_norm_count~Virus)

### Use Welch's ANOVA since variance are not equal in all groups
res.welch.anova <- norm_viral_count.df.airway %>% filter(log_norm_count!=-Inf) %>% welch_anova_test(log_norm_count~Virus)
res.welch.anova

norm_viral_count.df_stats.airway <- norm_viral_count.df.airway %>% filter(log_norm_count!=-Inf) %>% games_howell_test(log_norm_count~Virus)
norm_viral_count.df_stats.airway <- norm_viral_count.df_stats.airway %>% add_xy_position(x="Virus", fun="max")
norm_viral_count.df_stats.airway

# pdf("viral_count_airway.pdf", width=3, height=6)
norm_viral_count_plot.airway <- ggplot(norm_viral_count.df.airway, aes(x=Virus, y=log_norm_count, fill=Virus), alpha=0.5)
norm_viral_count_plot.airway <- norm_viral_count_plot.airway + geom_jitter(aes(col=Virus), alpha=0.9, shape=16, size=3, position=position_jitterdodge(0.5))
norm_viral_count_plot.airway <- norm_viral_count_plot.airway + geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.75)
norm_viral_count_plot.airway <- norm_viral_count_plot.airway + labs(x="Virus", y="log normalized counts")
norm_viral_count_plot.airway <- norm_viral_count_plot.airway + guides(fill="none", col='none')
norm_viral_count_plot.airway <- norm_viral_count_plot.airway + theme_classic()
norm_viral_count_plot.airway <- norm_viral_count_plot.airway + theme(text=element_text(size=24), title=element_text(size=24), legend.title=element_text(face="bold"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
norm_viral_count_plot.airway <- norm_viral_count_plot.airway + ylim(3,7)
norm_viral_count_plot.airway <- norm_viral_count_plot.airway + stat_pvalue_manual(norm_viral_count.df_stats.airway, label="{scales::pvalue(p.adj, accuracy=0.0001)}", inherit.aes=FALSE, size=8)
norm_viral_count_plot.airway
# dev.off()

########## Calculate the normalized counts of SC2 and B236 in colon organoids (Figure 4C)
norm_viral_count.df.colon <- norm_viral_count.df %>% filter(cell_type=="Colon_organoid")
norm_viral_count.df.colon %>% filter(cell_type=="Colon_organoid") %>% group_by(Virus) %>% shapiro_test(log_norm_count)
norm_viral_count.df.colon %>% filter(cell_type=="Colon_organoid") %>% levene_test(log_norm_count~Virus)

### Use one-way ANOVA since all assumptions are valid
res.anova <- norm_viral_count.df.colon  %>% anova_test(log_norm_count~Virus)
res.anova

norm_viral_count.df_stats.colon <- norm_viral_count.df.colon %>% filter(log_norm_count!=-Inf) %>% tukey_hsd(log_norm_count~Virus, p.adjust.method="fdr") # FDR (alias of BH)
norm_viral_count.df_stats.colon <- norm_viral_count.df_stats.colon %>% add_xy_position(x="Virus", fun="max")
norm_viral_count.df_stats.colon

# pdf("viral_count_colon_organoid.pdf", width=3, height=6)
norm_viral_count_plot.colon <- ggplot(norm_viral_count.df.colon, aes(x=Virus, y=log_norm_count, fill=Virus), alpha=0.5)
norm_viral_count_plot.colon <- norm_viral_count_plot.colon + geom_jitter(aes(col=Virus), alpha=0.9, shape=16, size=3, position=position_jitterdodge(0.5))
norm_viral_count_plot.colon <- norm_viral_count_plot.colon + geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.75)
norm_viral_count_plot.colon <- norm_viral_count_plot.colon + labs(x="Virus", y="log normalized counts")
norm_viral_count_plot.colon <- norm_viral_count_plot.colon + guides(fill="none", col='none')
norm_viral_count_plot.colon <- norm_viral_count_plot.colon + theme_classic()
norm_viral_count_plot.colon <- norm_viral_count_plot.colon + theme(text=element_text(size=24), title=element_text(size=24), legend.title=element_text(face="bold"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
norm_viral_count_plot.colon <- norm_viral_count_plot.colon + ylim(3,7)
norm_viral_count_plot.colon <- norm_viral_count_plot.colon + stat_pvalue_manual(norm_viral_count.df_stats.colon, label="{scales::pvalue(p.adj, accuracy=0.0001)}", inherit.aes=FALSE, size=8)
norm_viral_count_plot.colon
# dev.off()



########## Generate a table reporting transcriptome read alignment rates for airway epithelial cells and colon organoids (Supplementary Table 8)
samples.p <- samples %>% filter(cell.type %in% c("Airway","Colon_organoid"))
samples.p$star_stats_file <- paste("221202_STAR_alignment_stats_cell/", samples.p$fastq_f, "_Log.final.out", sep="")
star_alignment_stats <- read.delim(file.path(dir, samples.p$star_stats_file[1]), header=FALSE)
colnames(star_alignment_stats) <- c("Parameter", samples.p$fastq_f[1])

for (a in 1:length(samples.p$star_stats_file)) {
  if (a != 1) {
    alignment_stats <- read.delim(file.path(dir, samples.p$star_stats_file[a]), header=FALSE)
    star_alignment_stats[samples.p$fastq_f[a]] <- alignment_stats[2]
  }
}

rownames(star_alignment_stats) <- star_alignment_stats[,1]
star_alignment_stats[,1] <- NULL
star_alignment_stats <- as.data.frame(t(star_alignment_stats))
star_alignment_stats$pct_mapped <- (as.numeric(star_alignment_stats[,8]) + as.numeric(star_alignment_stats[,23])) / as.numeric(star_alignment_stats[,5]) * 100
min(star_alignment_stats$pct_mapped)

star_alignment_stats.df <- star_alignment_stats[,c(5, 8, 23, 34)]
# write.table(star_alignment_stats.df, "table_s8_sample_quality_mapped_reads_airway_colon.tsv", col.names=T, row.names=T, sep="\t", quote=F)



######################################## Analyze reads counts with separated viral gene reads ########################################
########## Read count table and record grouping information
viral_sep_count.df <- read.table(file.path(dir, "B236_airway_lung_colon_counts_viral_genes_sep.txt"), header=TRUE)
viral_sep_count.df_no_viral_seq <- viral_sep_count.df[1:(nrow(viral_sep_count.df)-22),]

samples_viral_sep <- samples %>% mutate(fastq_f_dir = ifelse(cell.type=="Airway", paste("X221021_A01297_0339_AHN3TFDRX2.bam_edited.", fastq_f, ".", fastq_f, "_Aligned.out.bam", sep=""),
                                                      ifelse(cell.type=="Lung", paste("X221124_A00691_0609_BHLKHCDRX2.bam_edited.", fastq_f, ".", fastq_f, "_Aligned.out.bam", sep=""),
                                                      ifelse(cell.type=="Colon_organoid", paste("X221202_A01297_0354_AHNJ7YDRX2.bam_edited.", str_replace_all(fastq_f, "-", "."), ".", str_replace_all(fastq_f, "-", "."), "_Aligned.out.bam", sep=""), NA))))
sample_dict_viral_sep <- unique(samples_viral_sep[,c("fastq_f_dir","Name")])

colnames(viral_sep_count.df) <- c(colnames(viral_sep_count.df)[1:6],sample_dict_viral_sep$Name[match(colnames(viral_sep_count.df[7:ncol(viral_sep_count.df)]), sample_dict_viral_sep$fastq_f_dir)])
colnames(viral_sep_count.df_no_viral_seq) <- c(colnames(viral_sep_count.df_no_viral_seq)[1:6],sample_dict_viral_sep$Name[match(colnames(viral_sep_count.df_no_viral_seq[7:ncol(viral_sep_count.df_no_viral_seq)]), sample_dict_viral_sep$fastq_f_dir)])

viral_sep_count.df_deseq <- viral_sep_count.df
rownames(viral_sep_count.df_deseq) <- viral_sep_count.df_deseq$Geneid
viral_sep_count.df_deseq[1:6] <- NULL
viral_sep_count.df_deseq <- viral_sep_count.df_deseq[samples_viral_sep$Name]
viral_sep_count.df_no_viral_seq_deseq <- viral_sep_count.df_deseq[1:(nrow(viral_sep_count.df_deseq)-22),]



########## Run zFPKM analysis
##### Check distribution of FPKM in all samples
##### FPKM = Fragments Per Kilobase Million
fpkm.df <- viral_sep_count.df_no_viral_seq
fpkm.df[2:5] <- NULL
total_read_count_per_each_exp <- colSums(fpkm.df[,3:ncol(fpkm.df)])
fpkm.df[,3:ncol(fpkm.df)] <- (fpkm.df[,3:ncol(fpkm.df)]*1000000)/(fpkm.df[,2]/1000)
fpkm.df[,3:ncol(fpkm.df)] <- sweep(fpkm.df[,3:ncol(fpkm.df)],2,total_read_count_per_each_exp,FUN="/")
rownames(fpkm.df) <- fpkm.df$Geneid
fpkm.df[1:2] <- NULL
zfpkm.df <- zFPKM(fpkm.df)
zfpkm_plot <- zFPKMPlot(fpkm.df, FacetTitles=TRUE)
zfpkm_plot

########## Create TPM table
tpm.df <- fpkm.df
total_fpkm_per_each_exp <- colSums(tpm.df)
tpm.df <- sweep((tpm.df*1000000),2,total_fpkm_per_each_exp,FUN="/")

########## Heatmap showing TPM of SC2 entry factor genes in the mock-infected airway epithelial cells and colon organoids (Supplementary Figure 4B)
average_tpm_mock.df <- data.frame(Airway = rowMeans(tpm.df[,c(13:16)]), Colon = rowMeans(tpm.df[,c(33:35)]))

entry_factor_gene_list <- c("ACE2","CTSB","CTSL","NRP1","TMEM106B","TMPRSS2")
entry_factor_gene_category.v <- c(rep("SC2 entry factor", length(entry_factor_gene_list)))
entry_factor_gene_category_plot_order <- c("SC2 entry factor")

# pdf("entry_factor_tpm_heatmap.pdf", width=3.307, height=2.2168)
entry_factor_tpm_heatmap <- Heatmap(t(log10(average_tpm_mock.df[entry_factor_gene_list,])), 
                                    show_row_names=FALSE, row_title_rot=90, row_names_gp=gpar(fontsize=24), 
                                    row_split=factor(c("Airway","Colon organoid"), levels=c("Airway","Colon organoid")), row_gap=unit(2.5, "mm"), cluster_rows=FALSE, row_title_gp=gpar(fontface="bold", fontsize=18), 
                                    column_split=factor(entry_factor_gene_category.v, levels=entry_factor_gene_category_plot_order), column_gap=unit(2.5, "mm"), cluster_columns=FALSE, clustering_method_columns="ward.D2", 
                                    column_title_gp=gpar(fontface="bold", fontsize=18), column_names_gp=gpar(fontsize=14, fontface="italic"), 
                                    col=colorRamp2(c(0, 3), c("white","red")), na_col="lightgrey", heatmap_legend_param=list(title="log10 TPM", title_gp=gpar(fontface="bold", fontsize=16), labels_gp=gpar(fontsize=16)))
entry_factor_tpm_heatmap
# dev.off()

########## Heatmap showing TPM of protease genes in the mock-infected airway epithelial cells and colon organoids (Supplementary Figure 4C)
tpm_gene_list_mock <- c("HPN","TMPRSS2","TMPRSS3","TMPRSS4","TMPRSS5","TMPRSS13","TMPRSS15","TMPRSS11A","TMPRSS11B","TMPRSS11D","TMPRSS11E","TMPRSS11F","TMPRSS6","TMPRSS7","TMPRSS9","ST14","PRSS8","PRSS21","TPSG1","CORIN",
                        "MMP14","MMP15","MMP16","MMP17","MMP24","MMP25",
                        "ADAM8","ADAM9","ADAM10","ADAM12","ADAM15","ADAM17","ADAM19","ADAM33")
tpm_gene_category.v_mock <- c(rep("Serine protease", 20), rep("MT-MMP", 6), rep("ADAM", 8))
tpm_gene_category_plot_order_mock <- c("Serine protease","MT-MMP","ADAM")

# pdf("protease_tpm_heatmap_airway_colon.pdf", width=16, height=2.25)
tpm_heatmap_mock <- Heatmap(t(log10(average_tpm_mock.df[tpm_gene_list_mock,])), 
                            show_row_names=FALSE, row_title_rot = 0, row_names_side="left", row_names_gp=gpar(fontsize=18), 
                            row_split=factor(c("Airway","Colon organoid"), levels=c("Airway","Colon organoid")), row_gap=unit(2.5, "mm"), cluster_rows=FALSE, clustering_method_rows="ward.D2", 
                            column_split=factor(tpm_gene_category.v_mock, levels=tpm_gene_category_plot_order_mock), column_gap=unit(2.5, "mm"), cluster_columns=FALSE, clustering_method_columns="ward.D2", 
                            column_title_gp=gpar(fontface="bold", fontsize=14), column_names_gp=gpar(fontsize=10, fontface="italic"), 
                            col=colorRamp2(c(0, 3), c("white","red")), na_col="lightgrey", heatmap_legend_param=list(title="log10 TPM", title_gp=gpar(fontface="bold"), labels_gp=gpar(fontsize=12)))
tpm_heatmap_mock
# dev.off()



######################################## Perform pairwise DEG analyses using Wald's test ########################################
########## Create FC and log2FC matrices
log2FC.mx <- data.frame()
padj.mx <- data.frame()
for (a in c("Airway","Colon_organoid")) {
	##### Use "apeglm" method for LFC shrinkage
	dds_FC.mx <- DESeqDataSetFromMatrix(countData=viral_sep_count.df_no_viral_seq_deseq[,(samples_viral_sep %>% filter(cell.type==a))$Name], colData=(samples_viral_sep %>% filter(cell.type==a)), design=~Virus)
	dds_FC.mx$Virus <- relevel(dds_FC.mx$Virus, "Mock_infected")
	dds_FC.mx <- DESeq(dds_FC.mx)
	for (b in 2:length(resultsNames(dds_FC.mx))) {
		resLFC_FC.mx <- lfcShrink(dds_FC.mx, coef=resultsNames(dds_FC.mx)[b], type="apeglm")
		if (nrow(log2FC.mx)==0) {
			log2FC.mx <- data.frame(resLFC_FC.mx$log2FoldChange)
			rownames(log2FC.mx) <- rownames(resLFC_FC.mx)
			colnames(log2FC.mx) <- paste(a, str_split(str_split(resultsNames(dds_FC.mx)[b], "Virus_", simplify=T)[2], "_vs_", simplify=T)[1], sep=" ")
		} else {
			log2FC.mx[paste(a, str_split(str_split(resultsNames(dds_FC.mx)[b], "Virus_", simplify=T)[2], "_vs_", simplify=T)[1], sep=" ")] <- resLFC_FC.mx$log2FoldChange
		}
		if (nrow(padj.mx)==0) {
			padj.mx <- data.frame(resLFC_FC.mx$padj)
			rownames(padj.mx) <- rownames(resLFC_FC.mx)
			colnames(padj.mx) <- paste(a, str_split(str_split(resultsNames(dds_FC.mx)[b], "Virus_", simplify=T)[2], "_vs_", simplify=T)[1], sep=" ")
		} else {
			padj.mx[paste(a, str_split(str_split(resultsNames(dds_FC.mx)[b], "Virus_", simplify=T)[2], "_vs_", simplify=T)[1], sep=" ")] <- resLFC_FC.mx$padj
		}
	}
}



########## Create PCA plots showing the expression change profiles of airway epithelial cells and colon organoids (Figure 4E)
log2FC.mx_no_any_na <- log2FC.mx[!apply(log2FC.mx, 1, function(x) any (is.na(x))),]
samples_viral_sep_FC_pca <- data.frame(colnames(log2FC.mx_no_any_na), row.names=colnames(log2FC.mx_no_any_na))
samples_viral_sep_FC_pca$Tissue <- c(rep("Airway", 3), rep("Colon organoid", 3))
samples_viral_sep_FC_pca$Virus <- c(rep(c("B236","SC2","SC2dFCS"), 2))
samples_viral_sep_FC_pca$Virus <- factor(samples_viral_sep_FC_pca$Virus, levels=c("SC2","SC2dFCS","B236"))

df_pca <- pca(log2FC.mx_no_any_na, metadata=samples_viral_sep_FC_pca, center=TRUE, scale=TRUE) # Standardization = Center then Scale

# pdf("pca_fc_airway_and_colon.pdf", width=10, height=10)
log2_fc_pca_plot <- biplot(df_pca, colby="Virus", shape="Tissue", xlab=paste0("PC1, ", round(df_pca$variance["PC1"], digits=2), "% variation"), ylab=paste0("PC2, ", round(df_pca$variance["PC2"], digits=2), "% variation"), 
                                 lab=NULL, legendTitleSize=24, legendLabSize=24, axisLabSize=24, pointSize=6, widthConnectors=1, legendPosition="right")
log2_fc_pca_plot <- log2_fc_pca_plot + theme(legend.title=element_text(face="bold"), plot.margin = margin(0, 0, 0, 1, "cm"), panel.grid.major=element_blank(), panel.grid.minor=element_blank())
log2_fc_pca_plot <- log2_fc_pca_plot + scale_shape_manual(values=15:17)
log2_fc_pca_plot <- log2_fc_pca_plot + coord_fixed()
log2_fc_pca_plot
# dev.off()



########## Assign inactive/actively expressed genes
mock_airway_inactive_genes <- rownames(zfpkm.df %>% select((samples_viral_sep %>% filter(grouping=="Mock_infected.Airway"))$Name) %>% filter(!rowMeans(.) > -3))
B236_airway_inactive_genes <- rownames(zfpkm.df %>% select((samples_viral_sep %>% filter(grouping=="B236.Airway"))$Name) %>% filter(!rowMeans(.) > -3))
SC2_airway_inactive_genes <- rownames(zfpkm.df %>% select((samples_viral_sep %>% filter(grouping=="SC2.Airway"))$Name) %>% filter(!rowMeans(.) > -3))
SC2dFCS_airway_inactive_genes <- rownames(zfpkm.df %>% select((samples_viral_sep %>% filter(grouping=="SC2dFCS.Airway"))$Name) %>% filter(!rowMeans(.) > -3))

mock_B236_airway_inactive_genes <- intersect(mock_airway_inactive_genes, B236_airway_inactive_genes)
mock_SC2_airway_inactive_genes <- intersect(mock_airway_inactive_genes, SC2_airway_inactive_genes)
mock_SC2dFCS_airway_inactive_genes <- intersect(mock_airway_inactive_genes, SC2dFCS_airway_inactive_genes)

mock_colon_inactive_genes <- rownames(zfpkm.df %>% select((samples_viral_sep %>% filter(grouping=="Mock_infected.Colon_organoid"))$Name) %>% filter(!rowMeans(.) > -3))
B236_colon_inactive_genes <- rownames(zfpkm.df %>% select((samples_viral_sep %>% filter(grouping=="B236.Colon"))$Name) %>% filter(!rowMeans(.) > -3))
SC2_colon_inactive_genes <- rownames(zfpkm.df %>% select((samples_viral_sep %>% filter(grouping=="SC2.Colon"))$Name) %>% filter(!rowMeans(.) > -3))
SC2dFCS_colon_inactive_genes <- rownames(zfpkm.df %>% select((samples_viral_sep %>% filter(grouping=="SC2dFCS.Colon"))$Name) %>% filter(!rowMeans(.) > -3))

mock_B236_colon_inactive_genes <- intersect(mock_colon_inactive_genes, B236_colon_inactive_genes)
mock_SC2_colon_inactive_genes <- intersect(mock_colon_inactive_genes, SC2_colon_inactive_genes)
mock_SC2dFCS_colon_inactive_genes <- intersect(mock_colon_inactive_genes, SC2dFCS_colon_inactive_genes)



########## Create a heatmap for log2FC of selected host genes across all samples
log2FC.mx_no_all_na <- log2FC.mx[!apply(log2FC.mx, 1, function(x) all (is.na(x))),]
padj.mx_no_all_na <- padj.mx[!apply(log2FC.mx, 1, function(x) all (is.na(x))),]

inactive_gene_list <- list(mock_B236_airway_inactive_genes, mock_SC2_airway_inactive_genes, mock_SC2dFCS_airway_inactive_genes, 
                           mock_B236_colon_inactive_genes, mock_SC2_colon_inactive_genes, mock_SC2dFCS_colon_inactive_genes)
inactive.mx_no_all_na <- padj.mx_no_all_na

for (a in 1:ncol(inactive.mx_no_all_na)) {inactive.mx_no_all_na[,a] <- as.integer(rownames(inactive.mx_no_all_na) %in% inactive_gene_list[[a]])}



########## Create a heatmap for log2FC of selected ISGs in the infected airway epithelial cells and colon organoids (Figures 1H, 4F)
selected_isg_gene_list <- sort(c("IFNA1","IFNB1","IFNAR1","IFNAR2","IFITM3","ISG15","JAK1","STAT2","IRF3","IRF7","NFKB1","MAVS","DDX58","TBK1","OAS1","PLSCR1","LY6E","IFIH1")) # RIGI is DDX58 in this GFF3 file

# pdf("selected_isg_deg_heatmap_airway_colon.pdf", width=5.35, height=9.785)
selected_isg_log2FC_heatmap_ha <- columnAnnotation(Virus=rep(c("SC2","SC2dFCS","B236"), 2), col=list(Virus=c("SC2"="yellow", "SC2dFCS"="navyblue", "B236"="hotpink")), 
                                                annotation_name_gp=gpar(fontsize=16, fontface="bold"), annotation_legend_param=list(title_gp=gpar(fontsize=16, fontface="bold"), at=c("SC2","SC2dFCS","B236"), labels_gp=gpar(fontsize=16)))
selected_isg_log2FC_heatmap <- Heatmap(log2FC.mx_no_all_na[,c(2,3,1,5,6,4)][selected_isg_gene_list,], top_annotation=selected_isg_log2FC_heatmap_ha, 
                                       row_names_side="left", row_names_gp=gpar(fontsize=14, fontface="italic"), cluster_rows=FALSE, row_title_gp=gpar(fontface="bold", fontsize=18), clustering_method_rows="ward.D2", 
                                       show_column_names=TRUE, column_title_rot=0, column_split=factor(c(rep("Airway", 3), rep("Colon organoid", 3)), levels=c("Airway","Colon organoid")), 
                                       column_names_side="top", column_names_gp=gpar(fontsize=24), cluster_columns=FALSE, column_title_gp=gpar(fontface="bold", fontsize=18), clustering_method_columns="ward.D2", 
                                       col=colorRamp2(c(-2, 0, 2), c("blue","white","red")), na_col="lightgrey", heatmap_legend_param=list(title="log2FC", title_gp=gpar(fontface="bold", fontsize=16), labels_gp=gpar(fontsize=16)), 
                                       cell_fun=function(j, i, x, y, width, height, fill) {if (!is.na(padj.mx_no_all_na[,c(2,3,1,5,6,4)][selected_isg_gene_list,][i,j]) 
                                                                                               && padj.mx_no_all_na[,c(2,3,1,5,6,4)][selected_isg_gene_list,][i,j]<0.05 
                                                                                               && inactive.mx_no_all_na[,c(2,3,1,5,6,4)][selected_isg_gene_list,][i,j]==0) grid.text("*", x, y, gp=gpar(fontsize=16), vjust=+0.75)})
selected_isg_log2FC_heatmap
# dev.off()



########## Cellular response to type I interferon (GO:0071357) # identical to Type I interferon signaling pathway (GO:0060337) in 2021
GO_0071357_gene_list <- c("MX1","RSAD2","MX2","OAS1","TYK2","BST2","IRF1","OAS2","IRF4","IRF5","OAS3","IFI27","IRF2","IRF3","IRF8","IRF9","IRF6","IRF7","IFITM2","IFITM3","SP100","IFITM1","IFI6","ADAR","IRAK1","IP6K2","GBP2","EGR1",
                          "IFNA14","STAT1","IFNA16","STAT2","IFNA17","HLA-C","HLA-A","ISG15","HLA-B","HLA-G","HLA-E","HLA-F","ISG20","MYD88","XAF1","IFNAR1","IFNA10","IFNAR2","IFNA6","IFNA5","IFNA8","RNASEL","IFNA7","IFNA2","IFNA1",
                          "IFNA4","IFI35","IFIT5","IFIT2","IFIT1","OASL","SAMHD1","IFIT3","PSMB8","IFNB1","JAK1","IFNA21")

########## Plot log2FC of ISGs in the infected airway epithelial cells (Figure 1G)
log2FC.mx_no_all_na.airway.isg <- as.data.frame(as.table(as.matrix(log2FC.mx_no_all_na[GO_0071357_gene_list,1:3])))
colnames(log2FC.mx_no_all_na.airway.isg) <- c("Gene","Virus","log2FC")
log2FC.mx_no_all_na.airway.isg$Virus <- c(rep("B236",length(GO_0071357_gene_list)), rep("SC2",length(GO_0071357_gene_list)), rep("SC2dFCS",length(GO_0071357_gene_list)))

padj.mx_no_all_na.airway.isg <- as.data.frame(as.table(as.matrix(padj.mx_no_all_na[GO_0071357_gene_list,1:3])))
colnames(padj.mx_no_all_na.airway.isg) <- c("Gene","Virus","padj")
padj.mx_no_all_na.airway.isg$Virus <- c(rep("B236",length(GO_0071357_gene_list)), rep("SC2",length(GO_0071357_gene_list)), rep("SC2dFCS",length(GO_0071357_gene_list)))

inactive.mx_no_all_na.airway.isg <- as.data.frame(as.table(as.matrix(inactive.mx_no_all_na[GO_0071357_gene_list,1:3])))
colnames(inactive.mx_no_all_na.airway.isg) <- c("Gene","Virus","inactive")
inactive.mx_no_all_na.airway.isg$Virus <- c(rep("B236",length(GO_0071357_gene_list)), rep("SC2",length(GO_0071357_gene_list)), rep("SC2dFCS",length(GO_0071357_gene_list)))

log2FC.mx_no_all_na.airway.isg$padj <- padj.mx_no_all_na.airway.isg$padj
log2FC.mx_no_all_na.airway.isg$inactive <- inactive.mx_no_all_na.airway.isg$inactive
log2FC.mx_no_all_na.airway.isg$DEG <- NA
for (a in 1:nrow(log2FC.mx_no_all_na.airway.isg)) {if (isTRUE(log2FC.mx_no_all_na.airway.isg[a,]$padj < 0.05 && log2FC.mx_no_all_na.airway.isg[a,]$inactive=="0")==TRUE) {log2FC.mx_no_all_na.airway.isg[a,]$DEG <- "DEG"} else {log2FC.mx_no_all_na.airway.isg[a,]$DEG <- "non-DEG"}}

log2FC.mx_no_all_na.airway.isg$Virus <- factor(log2FC.mx_no_all_na.airway.isg$Virus, levels=c("SC2","SC2dFCS","B236"))

log2FC.mx_no_all_na.airway.isg %>% filter(log2FC!=-Inf) %>% group_by(Virus) %>% shapiro_test(log2FC)
log2FC.mx_no_all_na.airway.isg %>% filter(log2FC!=-Inf) %>% levene_test(log2FC~Virus)

### Use Kruskal-Wallist test since all assumptions are violated
res.kruskal <- log2FC.mx_no_all_na.airway.isg %>% filter(log2FC!=-Inf) %>% kruskal_test(log2FC~Virus)
res.kruskal

log2FC.mx_no_all_na.airway.isg_stats <- log2FC.mx_no_all_na.airway.isg %>% filter(log2FC!=-Inf) %>% dunn_test(log2FC~Virus, p.adjust.method="fdr") # FDR (alias of BH)
log2FC.mx_no_all_na.airway.isg_stats <- log2FC.mx_no_all_na.airway.isg_stats %>% add_xy_position(x="Virus", fun="max")
log2FC.mx_no_all_na.airway.isg_stats

# pdf("violin_ifn_airway.pdf", width=7.5, height=10)
ifn_violin_plot.airway.isg <- ggplot((log2FC.mx_no_all_na.airway.isg %>% filter(log2FC!=-Inf)), aes(x=Virus, y=log2FC, fill=Virus), alpha=0.5) 
ifn_violin_plot.airway.isg <- ifn_violin_plot.airway.isg + geom_jitter(aes(col=DEG), alpha=0.5, shape=16, size=3, position=position_jitter(0.25))
ifn_violin_plot.airway.isg <- ifn_violin_plot.airway.isg + geom_boxplot(width=0.1, outlier.shape=NA, alpha=0.75)
ifn_violin_plot.airway.isg <- ifn_violin_plot.airway.isg + labs(title="GO:0071357", subtitle="cellular response to type I interferon", x="Virus", y=bquote(log[2]~FC))
ifn_violin_plot.airway.isg <- ifn_violin_plot.airway.isg + guides(fill="none")
ifn_violin_plot.airway.isg <- ifn_violin_plot.airway.isg + theme_classic()
ifn_violin_plot.airway.isg <- ifn_violin_plot.airway.isg + theme(text=element_text(size=24), plot.title=element_text(size=24), plot.subtitle=element_text(size=14), axis.title=element_text(size=24), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ifn_violin_plot.airway.isg <- ifn_violin_plot.airway.isg + ylim(-1,20)
ifn_violin_plot.airway.isg <- ifn_violin_plot.airway.isg + stat_pvalue_manual(log2FC.mx_no_all_na.airway.isg_stats, label="{scales::pvalue(p.adj, accuracy=0.0001)}", inherit.aes=FALSE, size=8)
ifn_violin_plot.airway.isg
# dev.off()

########## Plot log2FC of ISGs in the infected colon organoids (Supplementary Figure 5A)
log2FC.mx_no_all_na.colon.isg <- as.data.frame(as.table(as.matrix(log2FC.mx_no_all_na[GO_0071357_gene_list,4:6])))
colnames(log2FC.mx_no_all_na.colon.isg) <- c("Gene","Virus","log2FC")
log2FC.mx_no_all_na.colon.isg$Virus <- c(rep("B236",length(GO_0071357_gene_list)), rep("SC2",length(GO_0071357_gene_list)), rep("SC2dFCS",length(GO_0071357_gene_list)))
log2FC.mx_no_all_na.colon.isg$Virus <- factor(log2FC.mx_no_all_na.colon.isg$Virus, levels=c("SC2","SC2dFCS","B236"))

padj.mx_no_all_na.colon.isg <- as.data.frame(as.table(as.matrix(padj.mx_no_all_na[GO_0071357_gene_list,4:6])))
colnames(padj.mx_no_all_na.colon.isg) <- c("Gene","Virus","padj")
padj.mx_no_all_na.colon.isg$Virus <- c(rep("B236",length(GO_0071357_gene_list)), rep("SC2",length(GO_0071357_gene_list)), rep("SC2dFCS",length(GO_0071357_gene_list)))

inactive.mx_no_all_na.colon.isg <- as.data.frame(as.table(as.matrix(inactive.mx_no_all_na[GO_0071357_gene_list,4:6])))
colnames(inactive.mx_no_all_na.colon.isg) <- c("Gene","Virus","inactive")
inactive.mx_no_all_na.colon.isg$Virus <- c(rep("B236",length(GO_0071357_gene_list)), rep("SC2",length(GO_0071357_gene_list)), rep("SC2dFCS",length(GO_0071357_gene_list)))

log2FC.mx_no_all_na.colon.isg$padj <- padj.mx_no_all_na.colon.isg$padj
log2FC.mx_no_all_na.colon.isg$inactive <- inactive.mx_no_all_na.colon.isg$inactive
log2FC.mx_no_all_na.colon.isg$DEG <- NA
for (a in 1:nrow(log2FC.mx_no_all_na.colon.isg)) {if (isTRUE(log2FC.mx_no_all_na.colon.isg[a,]$padj < 0.05 && log2FC.mx_no_all_na.colon.isg[a,]$inactive=="0")==TRUE) {log2FC.mx_no_all_na.colon.isg[a,]$DEG <- "DEG"} else {log2FC.mx_no_all_na.colon.isg[a,]$DEG <- "non-DEG"}}

log2FC.mx_no_all_na.colon.isg %>% filter(log2FC!=-Inf) %>% group_by(Virus) %>% shapiro_test(log2FC)
log2FC.mx_no_all_na.colon.isg %>% filter(log2FC!=-Inf) %>% levene_test(log2FC~Virus)

### Use Kruskal-Wallist test since all assumptions are violated
res.kruskal <- log2FC.mx_no_all_na.colon.isg %>% filter(log2FC!=-Inf) %>% kruskal_test(log2FC~Virus)
res.kruskal

log2FC.mx_no_all_na.colon.isg_stats <- log2FC.mx_no_all_na.colon.isg %>% filter(log2FC!=-Inf) %>% dunn_test(log2FC~Virus, p.adjust.method="fdr") # FDR (alias of BH)
log2FC.mx_no_all_na.colon.isg_stats <- log2FC.mx_no_all_na.colon.isg_stats %>% add_xy_position(x="Virus", fun="max")
log2FC.mx_no_all_na.colon.isg_stats

# pdf("violin_ifn_colon.pdf", width=7.5, height=10)
ifn_violin_plot.colon.isg <- ggplot((log2FC.mx_no_all_na.colon.isg %>% filter(log2FC!=-Inf)), aes(x=Virus, y=log2FC, fill=Virus), alpha=0.5) 
ifn_violin_plot.colon.isg <- ifn_violin_plot.colon.isg + geom_jitter(aes(col=DEG), alpha=0.5, shape=16, size=3, position=position_jitter(0.25)) 
ifn_violin_plot.colon.isg <- ifn_violin_plot.colon.isg + geom_boxplot(width=0.1, outlier.shape=NA, alpha=0.75) 
ifn_violin_plot.colon.isg <- ifn_violin_plot.colon.isg + labs(title="GO:0071357", subtitle="cellular response to type I interferon", x="Virus", y=bquote(log[2]~FC)) 
ifn_violin_plot.colon.isg <- ifn_violin_plot.colon.isg + guides(fill="none") 
ifn_violin_plot.colon.isg <- ifn_violin_plot.colon.isg + theme_classic() 
ifn_violin_plot.colon.isg <- ifn_violin_plot.colon.isg + theme(text=element_text(size=24), plot.title=element_text(size=24), plot.subtitle=element_text(size=14), axis.title=element_text(size=24), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
ifn_violin_plot.colon.isg <- ifn_violin_plot.colon.isg + scale_y_continuous(limits=c(-4,12), breaks=c(-4,0,4,8,12)) 
ifn_violin_plot.colon.isg <- ifn_violin_plot.colon.isg + stat_pvalue_manual(log2FC.mx_no_all_na.colon.isg_stats, label="{scales::pvalue(p.adj, accuracy=0.0001)}", inherit.aes=FALSE, size=8) 
ifn_violin_plot.colon.isg
# dev.off()



########## Plot log2FC of intestinal marker genes (Figure 4G)
colon_marker_gene_list <- c("SI","SLC51A","SLC51B","VIL1","MUC2","TFF3","CHGA","REG4","GP2","LGR5")
colon_marker_gene_category.v <- c(rep("Colonocyte", 4), rep("Goblet cell", 2), rep("Enteroendocrine cell", 2), "Tuft cell", "Colonic stem cell")
colon_marker_gene_category_plot_order <- c("Goblet cell","Colonocyte","Enteroendocrine cell","Tuft cell","Colonic stem cell")

# pdf("heatmap_colon_organoid_marker_genes.pdf", width=6, height=6.753)
log2FC_heatmap_no_all_na.colon.marker <- Heatmap(log2FC.mx_no_all_na[,c(5,6,4)][colon_marker_gene_list,], 
                                                 row_names_side="left", row_split=factor(colon_marker_gene_category.v, levels=colon_marker_gene_category_plot_order), row_gap=unit(2.5, "mm"), row_title_gp=gpar(fontsize=18), row_title_rot=0, 
                                                 row_names_gp=gpar(fontsize=18, fontface="italic"), cluster_rows=FALSE, clustering_method_rows="ward.D2", 
                                                 column_names_side="top", show_column_names=TRUE, column_labels=c("SC2","SC2dFCS","B236"), column_names_gp=gpar(fontsize=18), cluster_columns=FALSE, clustering_method_columns="ward.D2", 
                                                 col=colorRamp2(c(-2, 0, 2), c("blue","white","red")), na_col="lightgrey", heatmap_legend_param=list(title="log2FC", title_gp=gpar(fontface="bold", fontsize=18), labels_gp=gpar(fontsize=18)), 
                                                 cell_fun=function(i, j, x, y, width, height, fill) {if (!is.na(padj.mx_no_all_na[,c(5,6,4)][colon_marker_gene_list,][j,i]) 
                                                                                                         && padj.mx_no_all_na[,c(5,6,4)][colon_marker_gene_list,][j,i]<0.05 
                                                                                                         && inactive.mx_no_all_na[,c(5,6,4)][colon_marker_gene_list,][j,i]==0) grid.text("*", x, y, gp=gpar(fontsize=24), vjust=+0.75)})
log2FC_heatmap_no_all_na.colon.marker
# dev.off()



######################################## Perform DEG analyses using Wald's test and likelihood ratio test ########################################\
########## 1. Analyze DEGs in airway experiment
dds_airway <- DESeqDataSetFromMatrix(countData=viral_sep_count.df_no_viral_seq_deseq %>% select((samples_viral_sep %>% filter(cell.type=="Airway"))$Name), colData=samples_viral_sep %>% filter(cell.type=="Airway"), design=~Virus)
dds_airway$Virus <- relevel(dds_airway$Virus, "Mock_infected")
dds_airway_wald <- DESeq(dds_airway) # Use for pairwise comparison
dds_airway_lrt <- DESeq(dds_airway, test="LRT", reduced=~1) # Use for multiple group comparison

####### 1.1 SC2 vs mock-infected
resLFC_SC2_airway_wald <- lfcShrink(dds_airway_wald, coef="Virus_SC2_vs_Mock_infected", type="apeglm")
sorted_resLFC_SC2_airway_wald <- resLFC_SC2_airway_wald[order(resLFC_SC2_airway_wald$padj),]
sorted_sig_resLFC_SC2_airway_wald <- sorted_resLFC_SC2_airway_wald[which(sorted_resLFC_SC2_airway_wald$padj < 0.05),]
filtered_sorted_sig_resLFC_SC2_airway_wald <- sorted_sig_resLFC_SC2_airway_wald[which(!rownames(sorted_sig_resLFC_SC2_airway_wald) %in% mock_SC2_airway_inactive_genes),]
SC2_airway_all_down_gene <- rownames(filtered_sorted_sig_resLFC_SC2_airway_wald[which(filtered_sorted_sig_resLFC_SC2_airway_wald$log2FoldChange < 0),])
SC2_airway_all_up_gene <- rownames(filtered_sorted_sig_resLFC_SC2_airway_wald[which(filtered_sorted_sig_resLFC_SC2_airway_wald$log2FoldChange > 0),])
summary(filtered_sorted_sig_resLFC_SC2_airway_wald)

####### 1.2 SC2dFCS vs mock-infected
resLFC_SC2dFCS_airway_wald <- lfcShrink(dds_airway_wald, coef="Virus_SC2dFCS_vs_Mock_infected", type="apeglm")
sorted_resLFC_SC2dFCS_airway_wald <- resLFC_SC2dFCS_airway_wald[order(resLFC_SC2dFCS_airway_wald$padj),]
sorted_sig_resLFC_SC2dFCS_airway_wald <- sorted_resLFC_SC2dFCS_airway_wald[which(sorted_resLFC_SC2dFCS_airway_wald$padj < 0.05),]
filtered_sorted_sig_resLFC_SC2dFCS_airway_wald <- sorted_sig_resLFC_SC2dFCS_airway_wald[which(!rownames(sorted_sig_resLFC_SC2dFCS_airway_wald) %in% mock_SC2dFCS_airway_inactive_genes),]
SC2dFCS_airway_all_down_gene <- rownames(filtered_sorted_sig_resLFC_SC2dFCS_airway_wald[which(filtered_sorted_sig_resLFC_SC2dFCS_airway_wald$log2FoldChange < 0),])
SC2dFCS_airway_all_up_gene <- rownames(filtered_sorted_sig_resLFC_SC2dFCS_airway_wald[which(filtered_sorted_sig_resLFC_SC2dFCS_airway_wald$log2FoldChange > 0),])
summary(filtered_sorted_sig_resLFC_SC2dFCS_airway_wald)

####### 1.3 B236 vs mock-infected
resLFC_B236_airway_wald <- lfcShrink(dds_airway_wald, coef="Virus_B236_vs_Mock_infected", type="apeglm")
sorted_resLFC_B236_airway_wald <- resLFC_B236_airway_wald[order(resLFC_B236_airway_wald$padj),]
sorted_sig_resLFC_B236_airway_wald <- sorted_resLFC_B236_airway_wald[which(sorted_resLFC_B236_airway_wald$padj < 0.05),]
filtered_sorted_sig_resLFC_B236_airway_wald <- sorted_sig_resLFC_B236_airway_wald[which(!rownames(sorted_sig_resLFC_B236_airway_wald) %in% mock_B236_airway_inactive_genes),]
B236_airway_all_down_gene <- rownames(filtered_sorted_sig_resLFC_B236_airway_wald[which(filtered_sorted_sig_resLFC_B236_airway_wald$log2FoldChange < 0),])
B236_airway_all_up_gene <- rownames(filtered_sorted_sig_resLFC_B236_airway_wald[which(filtered_sorted_sig_resLFC_B236_airway_wald$log2FoldChange > 0),])
summary(filtered_sorted_sig_resLFC_B236_airway_wald)

####### 1.4 Summarize the number of DEG to a bar plot (Supplementary Figure 2A)
airway_deg_summary.df <- data.frame(rep(c("B236","SC2","SC2dFCS"), 2), c(length(B236_airway_all_up_gene), length(SC2_airway_all_up_gene), length(SC2dFCS_airway_all_up_gene), 
                                                                         length(B236_airway_all_down_gene), length(SC2_airway_all_down_gene), length(SC2dFCS_airway_all_down_gene)), c(rep(c("Up"), 3), rep(c("Down"), 3)))
colnames(airway_deg_summary.df) <- c("Virus","Genes","Direction")
airway_deg_summary.df$Virus <- factor(airway_deg_summary.df$Virus, levels=c("SC2","SC2dFCS","B236"))
airway_deg_summary.df$Direction <- factor(airway_deg_summary.df$Direction, levels=c("Up","Down"))

# pdf("deg_num_airway.pdf", width=6, height=8)
airway_deg_summary_plot <- ggplot(data=airway_deg_summary.df, aes(x=Virus, y=Genes, fill=Direction)) 
airway_deg_summary_plot <- airway_deg_summary_plot + geom_bar(stat="identity", position="dodge", width=0.75) 
airway_deg_summary_plot <- airway_deg_summary_plot + theme_classic() 
airway_deg_summary_plot <- airway_deg_summary_plot + theme(text=element_text(size=24), strip.text=element_text(size=20, angle=270), strip.placement="outside", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
airway_deg_summary_plot <- airway_deg_summary_plot + scale_y_continuous(limits=c(0,7000), breaks=c(0,1000,2000,3000,4000,5000,6000,7000))
airway_deg_summary_plot
# dev.off()

####### 1.5 Create the table of DEG (Supplementary Table 1)
airway_deg_list <- list(rownames(filtered_sorted_sig_resLFC_SC2_airway_wald), rownames(filtered_sorted_sig_resLFC_SC2dFCS_airway_wald), rownames(filtered_sorted_sig_resLFC_B236_airway_wald))
mx_airway_deg_list <- max(lengths(airway_deg_list))
airway_deg.df <- as.data.frame(as.data.table(do.call(cbind, lapply(airway_deg_list, `length<-`, mx_airway_deg_list))))
colnames(airway_deg.df) <- c("SC2", "SC2dFCS", "B236")
# write.table(airway_deg.df, "table_s1_airway_deg.tsv", col.names=T, row.names=F, sep="\t", quote=F)

####### 1.6 Analyze the difference between B236 infection and SC2 or SC2dFCS infection in airway epithelial cells
dds_airway_B236_SC2 <- DESeqDataSetFromMatrix(countData=viral_sep_count.df_no_viral_seq_deseq[,((samples_viral_sep %>% filter(cell.type=="Airway", Virus!="SC2dFCS"))$Name)], 
                                              colData=samples_viral_sep %>% filter(cell.type=="Airway", Virus!="SC2dFCS"), design=~Virus)
dds_airway_B236_SC2$Virus <- relevel(dds_airway_B236_SC2$Virus, "Mock_infected")
dds_airway_B236_SC2_wald <- DESeq(dds_airway_B236_SC2) # Use for pairwise comparison
dds_airway_B236_SC2_lrt <- DESeq(dds_airway_B236_SC2, test="LRT", reduced=~1) # Use for multiple group comparison

resLFC_B236_SC2_airway_lrt <- lfcShrink(dds_airway_B236_SC2_lrt, coef="Virus_B236_vs_Mock_infected", type="apeglm")
sorted_resLFC_B236_SC2_airway_lrt <- resLFC_B236_SC2_airway_lrt[order(resLFC_B236_SC2_airway_lrt$padj),]
sorted_sig_resLFC_B236_SC2_airway_lrt <- sorted_resLFC_B236_SC2_airway_lrt[which(sorted_resLFC_B236_SC2_airway_lrt$padj < 0.05),]
summary(sorted_sig_resLFC_B236_SC2_airway_lrt)
airway_deg_list_lrt <- rownames(sorted_sig_resLFC_B236_SC2_airway_lrt)

##### 1.6.1 Calculate disco score
airway_disco <- data.frame(pvalue_B236_airway=resLFC_B236_airway_wald$pvalue, padj_B236_airway=resLFC_B236_airway_wald$padj, log2FC_B236_airway=resLFC_B236_airway_wald$log2FoldChange, 
                           pvalue_SC2_airway=resLFC_SC2_airway_wald$pvalue, padj_SC2_airway=resLFC_SC2_airway_wald$padj, log2FC_SC2_airway=resLFC_SC2_airway_wald$log2FoldChange, row.names=rownames(resLFC_B236_airway_wald))
airway_disco[,c("category","B236_plot_category","SC2_plot_category","pool_plot_category","disco_score")] <- NA

##### 1.6.2 Assign gene expression category based on directions of fold change
for (i in 1:nrow(airway_disco)) {
	if (is.na(airway_disco[i,]$pvalue_B236_airway) || is.na(airway_disco[i,]$padj_B236_airway) || rownames(airway_disco[i,]) %in% mock_B236_airway_inactive_genes || airway_disco[i,]$padj_B236_airway>=0.05) {
		B236_airway_direction <- "X"} else if (airway_disco[i,]$log2FC_B236_airway<0) {B236_airway_direction="D"} else {B236_airway_direction="U"}
	if (is.na(airway_disco[i,]$pvalue_SC2_airway) || is.na(airway_disco[i,]$padj_SC2_airway) || rownames(airway_disco[i,]) %in% mock_SC2_airway_inactive_genes || airway_disco[i,]$padj_SC2_airway>=0.05) {
		SC2_airway_direction <- "X"} else if (airway_disco[i,]$log2FC_SC2_airway<0) {SC2_airway_direction="D"} else {SC2_airway_direction="U"}
	airway_disco[i,]$category <- paste(B236_airway_direction, SC2_airway_direction, sep="")
	if (!rownames(airway_disco[i,]) %in% airway_deg_list_lrt) {airway_disco[i,]$category <- "XX"}
	if (airway_disco[i,]$category %in% c("DD")) {airway_disco[i,]$B236_plot_category <- "B236-Down/SC2-Down"}
	if (airway_disco[i,]$category %in% c("UU")) {airway_disco[i,]$B236_plot_category <- "B236-Up/SC2-Up"}
	if (airway_disco[i,]$category %in% c("DU","DX")) {airway_disco[i,]$B236_plot_category <- "B236-Down/SC2-notDown"}
	if (airway_disco[i,]$category %in% c("UD","UX")) {airway_disco[i,]$B236_plot_category <- "B236-Up/SC2-notUp"}
	if (is.na(airway_disco[i,]$B236_plot_category)) {airway_disco[i,]$B236_plot_category <- "OOI"}
	if (airway_disco[i,]$category %in% c("DD")) {airway_disco[i,]$SC2_plot_category <- "B236-Down/SC2-Down"}
	if (airway_disco[i,]$category %in% c("UU")) {airway_disco[i,]$SC2_plot_category <- "B236-Up/SC2-Up"}
	if (airway_disco[i,]$category %in% c("UD","XD")) {airway_disco[i,]$SC2_plot_category <- "B236-notDown/SC2-Down"}
	if (airway_disco[i,]$category %in% c("DU","XU")) {airway_disco[i,]$SC2_plot_category <- "B236-notUp/SC2-Up"}
	if (is.na(airway_disco[i,]$SC2_plot_category)) {airway_disco[i,]$SC2_plot_category <- "OOI"}
	airway_disco[i,]$pool_plot_category <- airway_disco[i,]$B236_plot_category
	if (airway_disco[i,]$SC2_plot_category=="B236-notDown/SC2-Down") {
		if (airway_disco[i,]$B236_plot_category=="B236-Up/SC2-notUp") {airway_disco[i,]$pool_plot_category <- "B236-Up/SC2-notUp or B236-notDown/SC2-Down"}
		if (airway_disco[i,]$B236_plot_category=="OOI") {airway_disco[i,]$pool_plot_category <- airway_disco[i,]$SC2_plot_category}
	}
	if (airway_disco[i,]$SC2_plot_category=="B236-notUp/SC2-Up") {
		if (airway_disco[i,]$B236_plot_category=="B236-Down/SC2-notDown") {airway_disco[i,]$pool_plot_category <- "B236-Down/SC2-notDown or B236-notUp/SC2-Up"}
		if (airway_disco[i,]$B236_plot_category=="OOI") {airway_disco[i,]$pool_plot_category <- airway_disco[i,]$SC2_plot_category}
	}
}

airway_disco$disco_score <- abs(airway_disco$log2FC_B236_airway * airway_disco$log2FC_SC2_airway * (log10(airway_disco$pvalue_B236_airway) + log10(airway_disco$pvalue_SC2_airway)))
sorted_airway_disco <- airway_disco[order(airway_disco$disco_score, decreasing=TRUE),]

sorted_top_B236_airway_disco <- Reduce(rbind, by(sorted_airway_disco, sorted_airway_disco["B236_plot_category"], head, n=100))
sorted_top_SC2_airway_disco <- Reduce(rbind, by(sorted_airway_disco, sorted_airway_disco["SC2_plot_category"], head, n=100))

sorted_airway_disco$pool_plot_color_category <- sorted_airway_disco$pool_plot_category
for (i in 1:nrow(sorted_airway_disco)) {if (!rownames(sorted_airway_disco[i,]) %in% c(rownames(sorted_top_B236_airway_disco), rownames(sorted_top_SC2_airway_disco))) {sorted_airway_disco[i,]$pool_plot_color_category <- "Out of interest"}}

sorted_airway_disco$pool_plot_color_category <- factor(sorted_airway_disco$pool_plot_color_category, levels=c("B236-Up/SC2-Up","B236-Down/SC2-Down","B236-Up/SC2-notUp","B236-Down/SC2-notDown","B236-notUp/SC2-Up","B236-notDown/SC2-Down",
                                                                                                              "B236-Up/SC2-notUp and B236-notDown/SC2-Down","B236-Down/SC2-notDown and B236-notUp/SC2-Up","Out of interest","OOI"))
sorted_airway_disco <- rbind(sorted_airway_disco[which(sorted_airway_disco$pool_plot_color_category=="Out of interest"),], sorted_airway_disco[which(sorted_airway_disco$pool_plot_color_category!="Out of interest"),])

##### 1.6.3 Create a representative plot showing gene expression category (Supplementary Figure 2B)
# pdf("deg_scatter_airway.pdf", width=10, height=10)
airway_fc_scatter_plot <- ggplot(sorted_airway_disco %>% filter(!pool_plot_color_category %in% c("Out of interest","OOI")), aes(x=log2FC_B236_airway, y=log2FC_SC2_airway, shape=pool_plot_color_category, color=pool_plot_color_category)) 
airway_fc_scatter_plot <- airway_fc_scatter_plot + geom_vline(xintercept=0, linewidth=0.5, alpha=0.1) 
airway_fc_scatter_plot <- airway_fc_scatter_plot + geom_hline(yintercept=0, linewidth=0.5, alpha=0.1) 
airway_fc_scatter_plot <- airway_fc_scatter_plot + geom_point(alpha=0.75, size=3) 
airway_fc_scatter_plot <- airway_fc_scatter_plot + xlab(bquote(~log[2] ~ "fold change B236")) 
airway_fc_scatter_plot <- airway_fc_scatter_plot + ylab(bquote(~log[2] ~ "fold change SC2")) 
airway_fc_scatter_plot <- airway_fc_scatter_plot + labs(shape="Category", color="Category") 
airway_fc_scatter_plot <- airway_fc_scatter_plot + scale_shape_manual(values=c(seq(15,18), c(21:23), 4)) 
airway_fc_scatter_plot <- airway_fc_scatter_plot + scale_color_manual(values=c(brewer.pal(6,"Set2")[c(1,2,3,5,4,6,7)], "lightgrey")) 
airway_fc_scatter_plot <- airway_fc_scatter_plot + theme_classic() 
airway_fc_scatter_plot <- airway_fc_scatter_plot + theme(axis.text=element_text(size=24), legend.position="right", text=element_text(size=24), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) 
airway_fc_scatter_plot <- airway_fc_scatter_plot + scale_x_continuous(limits=c(-4,12), breaks=c(-4,-2,0,2,4,6,8,10,12)) 
airway_fc_scatter_plot <- airway_fc_scatter_plot + scale_y_continuous(limits=c(-10,12), breaks=c(-10,-8,-6,-4,-2,0,2,4,6,8,10,12)) 
airway_fc_scatter_plot <- airway_fc_scatter_plot + coord_fixed()
airway_fc_scatter_plot
# dev.off()



########## 2. Analyze DEGs in the infected colon organoids
dds_colon <- DESeqDataSetFromMatrix(countData=viral_sep_count.df_no_viral_seq_deseq %>% select((samples_viral_sep %>% filter(cell.type=="Colon_organoid"))$Name), colData=samples_viral_sep %>% filter(cell.type=="Colon_organoid"), design=~Virus)
dds_colon$Virus <- relevel(dds_colon$Virus, "Mock_infected")
dds_colon_wald <- DESeq(dds_colon) # Use for pairwise comparison
dds_colon_lrt <- DESeq(dds_colon, test="LRT", reduced=~1) # Use for multiple group comparison

####### 2.1 SC2 vs mock-infected
resLFC_SC2_colon_wald <- lfcShrink(dds_colon_wald, coef="Virus_SC2_vs_Mock_infected", type="apeglm")
sorted_resLFC_SC2_colon_wald <- resLFC_SC2_colon_wald[order(resLFC_SC2_colon_wald$padj),]
sorted_sig_resLFC_SC2_colon_wald <- sorted_resLFC_SC2_colon_wald[which(sorted_resLFC_SC2_colon_wald$padj < 0.05),]
filtered_sorted_sig_resLFC_SC2_colon_wald <- sorted_sig_resLFC_SC2_colon_wald[which(!rownames(sorted_sig_resLFC_SC2_colon_wald) %in% mock_SC2_colon_inactive_genes),]
SC2_colon_all_down_gene <- rownames(filtered_sorted_sig_resLFC_SC2_colon_wald[which(filtered_sorted_sig_resLFC_SC2_colon_wald$log2FoldChange < 0),])
SC2_colon_all_up_gene <- rownames(filtered_sorted_sig_resLFC_SC2_colon_wald[which(filtered_sorted_sig_resLFC_SC2_colon_wald$log2FoldChange > 0),])
summary(filtered_sorted_sig_resLFC_SC2_colon_wald)

####### 2.2 SC2dFCS vs mock-infected
resLFC_SC2dFCS_colon_wald <- lfcShrink(dds_colon_wald, coef="Virus_SC2dFCS_vs_Mock_infected", type="apeglm")
sorted_resLFC_SC2dFCS_colon_wald <- resLFC_SC2dFCS_colon_wald[order(resLFC_SC2dFCS_colon_wald$padj),]
sorted_sig_resLFC_SC2dFCS_colon_wald <- sorted_resLFC_SC2dFCS_colon_wald[which(sorted_resLFC_SC2dFCS_colon_wald$padj < 0.05),]
filtered_sorted_sig_resLFC_SC2dFCS_colon_wald <- sorted_sig_resLFC_SC2dFCS_colon_wald[which(!rownames(sorted_sig_resLFC_SC2dFCS_colon_wald) %in% mock_SC2dFCS_colon_inactive_genes),]
SC2dFCS_colon_all_down_gene <- rownames(filtered_sorted_sig_resLFC_SC2dFCS_colon_wald[which(filtered_sorted_sig_resLFC_SC2dFCS_colon_wald$log2FoldChange < 0),])
SC2dFCS_colon_all_up_gene <- rownames(filtered_sorted_sig_resLFC_SC2dFCS_colon_wald[which(filtered_sorted_sig_resLFC_SC2dFCS_colon_wald$log2FoldChange > 0),])
summary(filtered_sorted_sig_resLFC_SC2dFCS_colon_wald)

####### 2.3 B236 vs mock-infected
resLFC_B236_colon_wald <- lfcShrink(dds_colon_wald, coef="Virus_B236_vs_Mock_infected", type="apeglm")
sorted_resLFC_B236_colon_wald <- resLFC_B236_colon_wald[order(resLFC_B236_colon_wald$padj),]
sorted_sig_resLFC_B236_colon_wald <- sorted_resLFC_B236_colon_wald[which(sorted_resLFC_B236_colon_wald$padj < 0.05),]
filtered_sorted_sig_resLFC_B236_colon_wald <- sorted_sig_resLFC_B236_colon_wald[which(!rownames(sorted_sig_resLFC_B236_colon_wald) %in% mock_B236_colon_inactive_genes),]
B236_colon_all_down_gene <- rownames(filtered_sorted_sig_resLFC_B236_colon_wald[which(filtered_sorted_sig_resLFC_B236_colon_wald$log2FoldChange < 0),])
B236_colon_all_up_gene <- rownames(filtered_sorted_sig_resLFC_B236_colon_wald[which(filtered_sorted_sig_resLFC_B236_colon_wald$log2FoldChange > 0),])
summary(filtered_sorted_sig_resLFC_B236_colon_wald)

####### 2.4 Summarize the number of DEG to a bar plot (Figure 4D)
colon_deg_summary.df <- data.frame(rep(c("B236","SC2","SC2dFCS"), 2), c(length(B236_colon_all_up_gene), length(SC2_colon_all_up_gene), length(SC2dFCS_colon_all_up_gene),
                                                                        length(B236_colon_all_down_gene), length(SC2_colon_all_down_gene), length(SC2dFCS_colon_all_down_gene)), c(rep(c("Up"), 3), rep(c("Down"), 3)))
colnames(colon_deg_summary.df) <- c("Virus","Genes","Direction")
colon_deg_summary.df$Virus <- factor(colon_deg_summary.df$Virus, levels=c("SC2","SC2dFCS","B236"))
colon_deg_summary.df$Direction <- factor(colon_deg_summary.df$Direction, levels=c("Up","Down"))

# pdf("deg_num_colon.pdf", width=6, height=8)
colon_deg_summary_plot <- ggplot(data=colon_deg_summary.df, aes(x=Virus, y=Genes, fill=Direction))
colon_deg_summary_plot <- colon_deg_summary_plot + geom_bar(stat="identity", position="dodge", width=0.75) 
colon_deg_summary_plot <- colon_deg_summary_plot + theme_classic() 
colon_deg_summary_plot <- colon_deg_summary_plot + theme(text=element_text(size=24), strip.text=element_text(size=20, angle=270), strip.placement="outside", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
colon_deg_summary_plot <- colon_deg_summary_plot + scale_y_continuous(limits=c(0,9000), breaks=c(0,1000,2000,3000,4000,5000,6000,7000,8000,9000))
colon_deg_summary_plot
# dev.off()

####### 2.5 Create the table of DEG (Supplementary Table 5)
colon_deg_list <- list(rownames(filtered_sorted_sig_resLFC_SC2_colon_wald), rownames(filtered_sorted_sig_resLFC_SC2dFCS_colon_wald), rownames(filtered_sorted_sig_resLFC_B236_colon_wald))
mx_colon_deg_list <- max(lengths(colon_deg_list))
colon_deg.df <- as.data.frame(as.data.table(do.call(cbind, lapply(colon_deg_list, `length<-`, mx_colon_deg_list))))
colnames(colon_deg.df) <- c("SC2", "SC2dFCS", "B236")
# write.table(colon_deg.df, "table_s1_colon_deg.tsv", col.names=T, row.names=F, sep="\t", quote=F)

####### 2.6 Analyze the difference between B236 infection and SC2 or SC2dFCS infection in colon organoids
dds_colon_B236_SC2 <- DESeqDataSetFromMatrix(countData=viral_sep_count.df_no_viral_seq_deseq[,((samples_viral_sep %>% filter(cell.type=="Colon_organoid", Virus!="SC2dFCS"))$Name)], 
                                             colData=samples_viral_sep %>% filter(cell.type=="Colon_organoid", Virus!="SC2dFCS"), design=~Virus)
dds_colon_B236_SC2$Virus <- relevel(dds_colon_B236_SC2$Virus, "Mock_infected")
dds_colon_B236_SC2_wald <- DESeq(dds_colon_B236_SC2) # Use for pairwise comparison
dds_colon_B236_SC2_lrt <- DESeq(dds_colon_B236_SC2, test="LRT", reduced=~1) # Use for multiple group comparison

resLFC_B236_SC2_colon_lrt <- lfcShrink(dds_colon_B236_SC2_lrt, coef="Virus_B236_vs_Mock_infected", type="apeglm")
sorted_resLFC_B236_SC2_colon_lrt <- resLFC_B236_SC2_colon_lrt[order(resLFC_B236_SC2_colon_lrt$padj),]
sorted_sig_resLFC_B236_SC2_colon_lrt <- sorted_resLFC_B236_SC2_colon_lrt[which(sorted_resLFC_B236_SC2_colon_lrt$padj < 0.05),]
summary(sorted_sig_resLFC_B236_SC2_colon_lrt)
colon_deg_list_lrt <- rownames(sorted_sig_resLFC_B236_SC2_colon_lrt)

##### 2.6.1 Calculate disco score
colon_disco <- data.frame(pvalue_B236_colon=resLFC_B236_colon_wald$pvalue, padj_B236_colon=resLFC_B236_colon_wald$padj, log2FC_B236_colon=resLFC_B236_colon_wald$log2FoldChange, 
                          pvalue_SC2_colon=resLFC_SC2_colon_wald$pvalue, padj_SC2_colon=resLFC_SC2_colon_wald$padj, log2FC_SC2_colon=resLFC_SC2_colon_wald$log2FoldChange, row.names=rownames(resLFC_B236_colon_wald))
colon_disco[,c("category","B236_plot_category","SC2_plot_category","pool_plot_category","disco_score")] <- NA

##### 2.6.2 Assign gene expression category based on directions of fold change
for (i in 1:nrow(colon_disco)) {
  if (is.na(colon_disco[i,]$pvalue_B236_colon) || is.na(colon_disco[i,]$padj_B236_colon) || rownames(colon_disco[i,]) %in% mock_B236_colon_inactive_genes || colon_disco[i,]$padj_B236_colon>=0.05) {
    B236_colon_direction <- "X"} else if (colon_disco[i,]$log2FC_B236_colon<0) {B236_colon_direction="D"} else {B236_colon_direction="U"}
  if (is.na(colon_disco[i,]$pvalue_SC2_colon) || is.na(colon_disco[i,]$padj_SC2_colon) || rownames(colon_disco[i,]) %in% mock_SC2_colon_inactive_genes || colon_disco[i,]$padj_SC2_colon>=0.05) {
    SC2_colon_direction <- "X"} else if (colon_disco[i,]$log2FC_SC2_colon<0) {SC2_colon_direction="D"} else {SC2_colon_direction="U"}
  colon_disco[i,]$category <- paste(B236_colon_direction, SC2_colon_direction, sep="")
  if (!rownames(colon_disco[i,]) %in% colon_deg_list_lrt) {colon_disco[i,]$category <- "XX"}
  if (colon_disco[i,]$category %in% c("DD")) {colon_disco[i,]$B236_plot_category <- "B236-Down/SC2-Down"}
  if (colon_disco[i,]$category %in% c("UU")) {colon_disco[i,]$B236_plot_category <- "B236-Up/SC2-Up"}
  if (colon_disco[i,]$category %in% c("DU","DX")) {colon_disco[i,]$B236_plot_category <- "B236-Down/SC2-notDown"}
  if (colon_disco[i,]$category %in% c("UD","UX")) {colon_disco[i,]$B236_plot_category <- "B236-Up/SC2-notUp"}
  if (is.na(colon_disco[i,]$B236_plot_category)) {colon_disco[i,]$B236_plot_category <- "OOI"}
  if (colon_disco[i,]$category %in% c("DD")) {colon_disco[i,]$SC2_plot_category <- "B236-Down/SC2-Down"}
  if (colon_disco[i,]$category %in% c("UU")) {colon_disco[i,]$SC2_plot_category <- "B236-Up/SC2-Up"}
  if (colon_disco[i,]$category %in% c("UD","XD")) {colon_disco[i,]$SC2_plot_category <- "B236-notDown/SC2-Down"}
  if (colon_disco[i,]$category %in% c("DU","XU")) {colon_disco[i,]$SC2_plot_category <- "B236-notUp/SC2-Up"}
  if (is.na(colon_disco[i,]$SC2_plot_category)) {colon_disco[i,]$SC2_plot_category <- "OOI"}
  colon_disco[i,]$pool_plot_category <- colon_disco[i,]$B236_plot_category
  if (colon_disco[i,]$SC2_plot_category=="B236-notDown/SC2-Down") {
    if (colon_disco[i,]$B236_plot_category=="B236-Up/SC2-notUp") {colon_disco[i,]$pool_plot_category <- "B236-Up/SC2-notUp or B236-notDown/SC2-Down"}
    if (colon_disco[i,]$B236_plot_category=="OOI") {colon_disco[i,]$pool_plot_category <- colon_disco[i,]$SC2_plot_category}
  }
  if (colon_disco[i,]$SC2_plot_category=="B236-notUp/SC2-Up") {
    if (colon_disco[i,]$B236_plot_category=="B236-Down/SC2-notDown") {colon_disco[i,]$pool_plot_category <- "B236-Down/SC2-notDown or B236-notUp/SC2-Up"}
    if (colon_disco[i,]$B236_plot_category=="OOI") {colon_disco[i,]$pool_plot_category <- colon_disco[i,]$SC2_plot_category}
  }
}

colon_disco$disco_score <- abs(colon_disco$log2FC_B236_colon * colon_disco$log2FC_SC2_colon * (log10(colon_disco$pvalue_B236_colon) + log10(colon_disco$pvalue_SC2_colon)))
sorted_colon_disco <- colon_disco[order(colon_disco$disco_score, decreasing=TRUE),]

sorted_top_B236_colon_disco <- Reduce(rbind, by(sorted_colon_disco, sorted_colon_disco["B236_plot_category"], head, n=100))
sorted_top_SC2_colon_disco <- Reduce(rbind, by(sorted_colon_disco, sorted_colon_disco["SC2_plot_category"], head, n=100))

sorted_colon_disco$pool_plot_color_category <- sorted_colon_disco$pool_plot_category
for (i in 1:nrow(sorted_colon_disco)) {if (!rownames(sorted_colon_disco[i,]) %in% c(rownames(sorted_top_B236_colon_disco), rownames(sorted_top_SC2_colon_disco))) {sorted_colon_disco[i,]$pool_plot_color_category <- "Out of interest"}}

sorted_colon_disco$pool_plot_color_category <- factor(sorted_colon_disco$pool_plot_color_category, levels=c("B236-Up/SC2-Up","B236-Down/SC2-Down","B236-Up/SC2-notUp","B236-Down/SC2-notDown","B236-notUp/SC2-Up","B236-notDown/SC2-Down",
                                                                                                            "B236-Up/SC2-notUp and B236-notDown/SC2-Down","B236-Down/SC2-notDown and B236-notUp/SC2-Up","Out of interest","OOI"))
sorted_colon_disco <- rbind(sorted_colon_disco[which(sorted_colon_disco$pool_plot_color_category=="Out of interest"),], sorted_colon_disco[which(sorted_colon_disco$pool_plot_color_category!="Out of interest"),])

# pdf("deg_scatter_colon.pdf", width=10, height=10)
colon_fc_scatter_plot <- ggplot(sorted_colon_disco %>% filter(!pool_plot_color_category %in% c("Out of interest","OOI")), aes(x=log2FC_B236_colon, y=log2FC_SC2_colon, shape=pool_plot_color_category, color=pool_plot_color_category)) 
colon_fc_scatter_plot <- colon_fc_scatter_plot + geom_vline(xintercept=0, linewidth=0.5, alpha=0.1) 
colon_fc_scatter_plot <- colon_fc_scatter_plot + geom_hline(yintercept=0, linewidth=0.5, alpha=0.1) 
colon_fc_scatter_plot <- colon_fc_scatter_plot + geom_point(alpha=0.75, size=3) 
colon_fc_scatter_plot <- colon_fc_scatter_plot + xlab(bquote(~log[2] ~ "fold change B236")) 
colon_fc_scatter_plot <- colon_fc_scatter_plot + ylab(bquote(~log[2] ~ "fold change SC2")) 
colon_fc_scatter_plot <- colon_fc_scatter_plot + labs(shape="Category", color="Category") 
colon_fc_scatter_plot <- colon_fc_scatter_plot + scale_shape_manual(values=c(15, 17, 18, 4)) 
colon_fc_scatter_plot <- colon_fc_scatter_plot + scale_color_manual(values=c(brewer.pal(6,"Set2")[c(1,3,5)], "lightgrey")) 
colon_fc_scatter_plot <- colon_fc_scatter_plot + theme_classic() 
colon_fc_scatter_plot <- colon_fc_scatter_plot + theme(axis.text=element_text(size=24), legend.position="right", text=element_text(size=24), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) 
colon_fc_scatter_plot <- colon_fc_scatter_plot + scale_x_continuous(limits=c(-8,12), breaks=c(-8,-6,-4,-2,0,2,4,6,8,10,12)) 
colon_fc_scatter_plot <- colon_fc_scatter_plot + scale_y_continuous(limits=c(-2,6), breaks=c(-2,0,2,4,6)) 
colon_fc_scatter_plot <- colon_fc_scatter_plot + coord_fixed()
colon_fc_scatter_plot
# dev.off()



########## 3. Perform enrichment analyses
####### 3.1 Create lists of the top 100 genes in the infected airway epithelial cells with highest disco scores
highest_airway_disco_dd_gene <- row.names(sorted_top_B236_airway_disco %>% filter(B236_plot_category=="B236-Down/SC2-Down"))
highest_airway_disco_uu_gene <- row.names(sorted_top_B236_airway_disco %>% filter(B236_plot_category=="B236-Up/SC2-Up"))
highest_airway_disco_du_dx_gene <- row.names(sorted_top_B236_airway_disco %>% filter(B236_plot_category=="B236-Down/SC2-notDown"))
highest_airway_disco_ud_ux_gene <- row.names(sorted_top_B236_airway_disco %>% filter(B236_plot_category=="B236-Up/SC2-notUp"))
highest_airway_disco_ud_xd_gene <- row.names(sorted_top_SC2_airway_disco %>% filter(SC2_plot_category=="B236-notDown/SC2-Down"))
highest_airway_disco_du_xu_gene <- row.names(sorted_top_SC2_airway_disco %>% filter(SC2_plot_category=="B236-notUp/SC2-Up"))
highest_colon_disco_dd_gene <- row.names(sorted_top_B236_colon_disco %>% filter(B236_plot_category=="B236-Down/SC2-Down"))
highest_colon_disco_uu_gene <- row.names(sorted_top_B236_colon_disco %>% filter(B236_plot_category=="B236-Up/SC2-Up"))
highest_colon_disco_du_dx_gene <- row.names(sorted_top_B236_colon_disco %>% filter(B236_plot_category=="B236-Down/SC2-notDown"))
highest_colon_disco_ud_ux_gene <- row.names(sorted_top_B236_colon_disco %>% filter(B236_plot_category=="B236-Up/SC2-notUp"))
highest_colon_disco_ud_xd_gene <- row.names(sorted_top_SC2_colon_disco %>% filter(SC2_plot_category=="B236-notDown/SC2-Down"))
highest_colon_disco_du_xu_gene <- row.names(sorted_top_SC2_colon_disco %>% filter(SC2_plot_category=="B236-notUp/SC2-Up"))

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
  	if (substr(CPfile@compareClusterResult$Cluster[i],1,1)=="C") {CPfile@compareClusterResult$group[i] <- "Colon organoid"}
  	if (substr(CPfile@compareClusterResult$Cluster[i],1,1)=="A") {CPfile@compareClusterResult$group[i] <- "Airway"}
  }
  CPfile@compareClusterResult$dummy_category <- str_split(CPfile@compareClusterResult$Cluster, ":", simplify=T)[,2]
  CPfile@compareClusterResult$Category <- paste(str_split(CPfile@compareClusterResult$Cluster, ":", simplify=T)[,2], "\n", "(", as.character(str_split(CPfile@compareClusterResult$GeneRatio, "/", simplify=T)[,2]), ")", sep="")
  CPfile@compareClusterResult$group <- factor(CPfile@compareClusterResult$group, levels=c("Colon organoid","Airway"))
  desired_level_factor <- intersect(c("B236-Up/SC2-Up","B236-Down/SC2-Down","B236-Up/SC2-notUp","B236-Down/SC2-notDown","B236-notUp/SC2-Up","B236-notDown/SC2-Down"), unique(CPfile@compareClusterResult$dummy_category))
  plot_level_factor <- unique(CPfile@compareClusterResult$Category)
  plot_level_factor <- plot_level_factor[sapply(desired_level_factor, function(x) { grep(x, plot_level_factor)})]
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

####### 3.2 Perform functional enrichment analysis with top 100 genes in the infected airway epithelial cells with highest disco scores (Supplementary Figure 2C)
compare_samples <- list(
  "A:B236-Up/SC2-Up" = highest_airway_disco_uu_gene,
  "A:B236-Down/SC2-Down" = highest_airway_disco_dd_gene,
  "A:B236-Up/SC2-notUp" = highest_airway_disco_ud_ux_gene,
  "A:B236-Down/SC2-notDown" = highest_airway_disco_du_dx_gene,
  "A:B236-notUp/SC2-Up" = highest_airway_disco_du_xu_gene,
  "A:B236-notDown/SC2-Down" = highest_airway_disco_ud_xd_gene
)

DIR <- "230309_compareCluster_airway_BvsW_top_100"
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

# pdf("enrichment_top_100_airway.pdf", width=20, height=15)
enrich_plot.airway_1 <- plot_grid(myplots[[3]], labels=c("B"), label_size=32, ncol=1)
enrich_plot.airway_2 <- plot_grid(myplots[[5]], myplots[[10]], labels=c("C","D"), label_size=32, ncol=1, align="v", axis="lr", rel_heights=c(0.32,0.68))
enrich_plot.airway <- plot_grid(enrich_plot.airway_1, enrich_plot.airway_2, nrow=1, rel_widths=c(0.54,0.46), align="hv")
enrich_plot.airway
# dev.off()

####### 3.3 Perform functional enrichment analysis with top 100 genes in the infected colon organoids with highest disco scores (Supplementary Figure 5B-D)
compare_samples <- list(
  "C:B236-Up/SC2-Up" = highest_colon_disco_uu_gene,
  "C:B236-Down/SC2-Down" = highest_colon_disco_dd_gene,
  "C:B236-Up/SC2-notUp" = highest_colon_disco_ud_ux_gene,
  "C:B236-Down/SC2-notDown" = highest_colon_disco_du_dx_gene,
  "C:B236-notUp/SC2-Up" = highest_colon_disco_du_xu_gene,
  "C:B236-notDown/SC2-Down" = highest_colon_disco_ud_xd_gene
)

DIR <- "230309_compareCluster_colon_BvsW_top_100"
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

# pdf("enrichment_top_100_colon.pdf", width=20, height=15)
enrich_plot.colon_1 <- plot_grid(myplots[[3]], labels=c("B"), label_size=32, ncol=1)
enrich_plot.colon_2 <- plot_grid(myplots[[5]], myplots[[10]], labels=c("C","D"), label_size=32, ncol=1, align="v", axis="lr", rel_heights=c(0.3875,0.6125))
enrich_plot.colon <- plot_grid(enrich_plot.colon_1, enrich_plot.colon_2, nrow=1, rel_widths=c(0.54,0.46), align="hv")
enrich_plot.colon
# dev.off()
