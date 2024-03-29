library(Biostrings)
library(dplyr)
library(stringr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

### Change working directory
dir <- "~/Desktop/BANAL_project/240307_MSA_identity/"
setwd(dir)

gene_name.v <- c("ORF1ab","ORF1a","S","ORF3a","E","M","ORF6","ORF7a","ORF7b","ORF8","N","ORF10")
cds_pid.v <- c()
prot_pid.v <- c()

################################################## Calculate PID between SC2 and B236 CDSs of each gene ##################################################
SC2_cds_fasta <- readDNAStringSet("240228_SC2_CDS.fasta", "fasta")
B236_cds_fasta <- readDNAStringSet("240228_B236_CDS.fasta", "fasta")

# names(SC2_cds_fasta) <- str_split(str_split(names(SC2_cds_fasta), "protein=", simplify = TRUE)[,2], "] ", simplify = TRUE)[,1]
# names(B236_cds_fasta) <- str_split(str_split(names(B236_cds_fasta), "protein=", simplify = TRUE)[,2], "] ", simplify = TRUE)[,1]
# names(SC2_cds_fasta) == names(B236_cds_fasta) ### Check order of proteins --> SAME!

names(SC2_cds_fasta) <- paste("SC2_", gene_name.v, sep="")
names(B236_cds_fasta) <- paste("B236_", gene_name.v, sep="")

for (i in 1:length(SC2_cds_fasta)) {
  cds_name <- gsub("SC2_", "", names(SC2_cds_fasta[i]), fixed=TRUE)
  cds_pair_fasta <- c(SC2_cds_fasta[names(SC2_cds_fasta[i])], B236_cds_fasta[names(B236_cds_fasta[i])])
  data_fasta <- data.frame(label=names(cds_pair_fasta), seq=paste(cds_pair_fasta))
  data_fasta$label <- paste('>', data_fasta$label, sep="")
  
  output_fasta <- paste("pair_", cds_name, "_cds.fasta", sep="")
  # write.table(data_fasta, output_fasta, col.names=F, row.names=F, sep="\n", quote=F)
  
  output_mafft <- paste("pair_", cds_name, "_cds.mafft", sep="")
  mafft_cmd <- paste("mafft --maxiterate 1000 --globalpair ", output_fasta, " > ", output_mafft, sep="")
  # system(mafft_cmd)
  
  mafft_fasta <- readDNAStringSet(output_mafft, "fasta")
  data_mafft <- data.frame(label=names(mafft_fasta), seq=paste(mafft_fasta))
  
  cds_pid.v <- append(cds_pid.v, mapply(function(x,y) sum(x==y), strsplit(data_mafft[1,2], ""), strsplit(data_mafft[2,2], ""))*100/(width(mafft_fasta)[1])) 
}

################################################## Calculate PID between SC2 and B236 proteins of each gene ##################################################
SC2_prot_fasta <- readAAStringSet("240228_SC2_prot.fasta", "fasta")
B236_prot_fasta <- readAAStringSet("240228_B236_prot.fasta", "fasta")

# names(SC2_prot_fasta) <- str_split(str_split(names(SC2_prot_fasta), "protein=", simplify = TRUE)[,2], "] ", simplify = TRUE)[,1]
# names(B236_prot_fasta) <- str_split(str_split(names(B236_prot_fasta), "protein=", simplify = TRUE)[,2], "] ", simplify = TRUE)[,1]
# names(SC2_prot_fasta) == names(B236_prot_fasta) ### Check order of proteins --> SAME!

names(SC2_prot_fasta) <- paste("SC2_", gene_name.v, sep="")
names(B236_prot_fasta) <- paste("B236_", gene_name.v, sep="")

for (i in 1:length(SC2_prot_fasta)) {
  prot_name <- gsub("SC2_", "", names(SC2_prot_fasta[i]), fixed=TRUE)
  prot_pair_fasta <- c(SC2_prot_fasta[names(SC2_prot_fasta[i])], B236_prot_fasta[names(B236_prot_fasta[i])])
  data_fasta <- data.frame(label=names(prot_pair_fasta), seq=paste(prot_pair_fasta))
  data_fasta$label <- paste('>', data_fasta$label, sep="")
  
  output_fasta <- paste("pair_", prot_name, "_prot.fasta", sep="")
  # write.table(data_fasta, output_fasta, col.names=F, row.names=F, sep="\n", quote=F)
  
  output_mafft <- paste("pair_", prot_name, "_prot.mafft", sep="")
  mafft_cmd <- paste("mafft --maxiterate 1000 --globalpair ", output_fasta, " > ", output_mafft, sep="")
  # system(mafft_cmd)
  
  mafft_fasta <- readAAStringSet(output_mafft, "fasta")
  data_mafft <- data.frame(label=names(mafft_fasta), seq=paste(mafft_fasta))
  
  prot_pid.v <- append(prot_pid.v, mapply(function(x,y) sum(x==y), strsplit(data_mafft[1,2], ""), strsplit(data_mafft[2,2], ""))*100/(width(mafft_fasta)[1])) 
}

pid.df <- data.frame(cds_pid=cds_pid.v, prot_pid=prot_pid.v, row.names=gene_name.v) %>% arrange(prot_pid)

# pdf("prot_identity_heatmap.pdf", width=3, height=8.025)
pid_heatmap <- Heatmap(pid.df[,2,drop=FALSE], row_names_gp=gpar(fontsize=20), row_gap=unit(2.5, "mm"), cluster_rows=FALSE, row_title_gp=gpar(fontface="bold", fontsize=18), column_gap=unit(2.5, "mm"), cluster_columns=FALSE, col=colorRamp2(c(90, 100), c("white","red")),
                       column_title_gp=gpar(fontface="bold", fontsize=18), column_names_gp=gpar(fontsize=14, fontface="italic"), clustering_method_rows="ward.D2", clustering_method_columns="ward.D2", na_col="lightgrey",
                       heatmap_legend_param=list(title="% identity", title_gp=gpar(fontface="bold", fontsize=16), labels_gp=gpar(fontsize=16)), cell_fun=function(j, i, x, y, width, height, fill) grid.text(round(pid.df[i,2,drop=FALSE], 0), x, y, gp=gpar(fontsize=16), vjust=+0.75))
pid_heatmap
# dev.off()
