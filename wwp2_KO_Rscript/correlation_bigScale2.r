# 
# Reference methods:
#  1. Single-cell transcriptomics unveils gene regulatory network plasticity
#  https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1713-4

options(stringsAsFactors = FALSE)
options(digits = 4)
#options(future.globals.maxSize = 16 *1000 * 1024^2) #memory

library(Seurat)
library(ggplot2)
library(pheatmap)
library(SingleCellExperiment)
#library(scater)
#library(scran)
library(dplyr)
#library(dismay)
#library(propr)
library(xlsx)
source("plot_pdf_png.r")



project_dir <- "project_path"

path_output <- paste0(project_dir,"results/")
#system(paste0("mkdir ",path_output))
path_output_fb <- paste0(path_output,"1.correlation_bigScale/") 
# system(paste0("mkdir ",path_output_fb))

path_deposit <- paste0(project_dir,"dataDeposit/")

path_input <- paste0(project_dir, "input_Robjects/")


fb12 <- readRDS(file = paste0(path_input, "fb12.RDS"))
#p <- DimPlot(fb, group.by = "final_celltype_refined", label = T)
#pdf_and_png(p, file_name = paste0(path_output_fb, "umap_fb"), h_to_w_ratio = 0.6)

#Idents(fb) <- "final_celltype_refined"


counts  <- GetAssayData(object = fb12 , slot = "counts")
#counts.wt  <- GetAssayData(object = fb12.WT , slot = "counts")
#counts.ko  <- GetAssayData(object = fb12.KO , slot = "counts")
gene.symbols <- data.matrix(rownames(counts))



# Run only once, from recursive results
################################################
library(bigSCale)

results.recursive = compute.network(expr.data = counts, gene.names = gene.symbols, clustering = "recursive")
results.direct = compute.network(expr.data = counts, gene.names = gene.symbols, clustering = "direct")

correlations =as.data.frame(as.double(results$correlations))


results.wt=compute.network(expr.data = counts.wt, gene.names = gene.symbols, clustering = "direct")
correlations.wt =as.data.frame(as.double(results.wt$correlations))


results.ko=compute.network(expr.data = counts.ko, gene.names = gene.symbols, clustering = "direct")
correlations.ko =as.data.frame(as.double(results.ko$correlations))


#maintain almost equal size of centrality nodes and edges between wt and ko 
output=homogenize.networks(list(results.wt,results.ko))
results.wt2=output[[1]]
results.ko2=output[[2]]

comparison=compare.centrality(list(results.wt2$centrality,results.ko2$centrality),c('WT','KO'))
for(i in names(comparison)){
  print(i)
  t <- as.data.frame(comparison[[i]])
  write.xlsx(t, file = paste0(path_output_fb, "ko_wt_comparison.xlsx"), sheetName = i, col.names = T, row.names = T, append = T)
}


save(results.wt, results.ko, results.wt2, results.ko2, comparison, file = paste0(path_output_fb, "bigScale_results_direct.RData"))
save(results.direct, results.recursive, paste0(path_deposit, "bigScale_fb12_all.RData"))

################################################
# note that correlation of gene pairs are identical between results.wt$correlation and results.wt2$correlation

library(float) # install of float package needs development tool (Xcode) in mac
#load(paste0(path_deposit,"bigScale_fb12_all.RData"))

#correlations.all = as.data.frame(as.double(results.recursive$correlations))
#saveRDS(correlations.all, paste0(path_deposit,"bigScale_fb12_correlation.RDS"))
correlations.all <- readRDS(paste0(path_deposit,"bigScale_fb12_correlation.RDS"))

ppargc1a.cor <- correlations.all[, "Ppargc1a"]
names(ppargc1a.cor) <- rownames(correlations.all)
ppargc1a.cor <- sort(ppargc1a.cor, decreasing = T)

ko.degenes <- readRDS(paste0(path_deposit, "KO_vs_WT.degenes.RDS"))

## test the difference between direct and recursive model
#correlations.direct.all <-  as.data.frame(as.double(results.direct$correlations))
#ppargc1a.direct.cor <- correlations.direct.all[, "Ppargc1a"]
#names(ppargc1a.direct.cor) <- rownames(correlations.direct.all)
#ppargc1a.direct.cor <- ppargc1a.direct.cor[names(ppargc1a.cor)]

#cor.test(ppargc1a.direct.cor,ppargc1a.cor)
## 0.9628, p-value <2e-16

ppargc1a.cor.df <- data.frame(cor  = ppargc1a.cor)
write.xlsx(ppargc1a.cor.df, file = paste0(path_output_fb, "ppargc1a.cor.xlsx"), row.names = T, col.names = T)
  
library(Hmisc)
zscore <- data.matrix(results.recursive$tot.scores)
colnames(zscore) <- colnames(correlations.all)
zscore.pearson <- rcorr(zscore,type = "pearson")
zscore.spearman <- rcorr(zscore, type = "spearman") 

  
################################################

ppargc1a.cor.df <- 


ppargc1a.cor.all <- data.frame(wt = ppargc1a.wt.cor[ppargc1a.cor.allgenes,], ko = ppargc1a.ko.cor[ppargc1a.cor.allgenes,], gene = ppargc1a.cor.allgenes)

ppargc1a.cor.all.df <- melt(ppargc1a.cor.all, id.vars = "gene", measure.vars = c("wt", "ko"))
colnames(ppargc1a.cor.all.df) <- c("gene", "condition", "cor")

# select a significant value 0.5,0.7
ppargc1a.wt.genes <- rownames(ppargc1a.wt.cor)[abs(ppargc1a.wt.cor$orginal)>0.5]
ppargc1a.wt.network.both <- data.frame( gene = ppargc1a.wt.genes, wt = ppargc1a.wt.cor[ppargc1a.wt.genes,], ko = ppargc1a.ko.cor[ppargc1a.wt.genes,])
ppargc1a.wt.cor.df <- melt(ppargc1a.wt.network.both, id.vars = "gene", measure.vars = c("wt", "ko"))
colnames(ppargc1a.wt.cor.df) <- c("gene", "condition", "cor")

write.xlsx(ppargc1a.wt.network.both, file = paste0(path_output_fb, "z-score_correlation_wt.0.7.xlsx"), row.names = F, col.names = T, append = F)
# select a significant value 0.5,0.7
ppargc1a.ko.genes <- rownames(ppargc1a.ko.cor)[abs(ppargc1a.ko.cor$orginal)>0.5]
ppargc1a.ko.network.both <- data.frame(wt = ppargc1a.wt.cor[ppargc1a.ko.genes,], ko = ppargc1a.ko.cor[ppargc1a.ko.genes,], gene = ppargc1a.ko.genes)
ppargc1a.ko.cor.df <- melt(ppargc1a.ko.network.both, id.vars = "gene", measure.vars = c("wt", "ko"))
colnames(ppargc1a.ko.cor.df) <- c("gene", "condition", "cor")

################################################

# density plot

library(wesanderson)
p1 <- ggplot(ppargc1a.wt.cor.df, aes(x=cor, color=condition, fill=condition)) + 
   scale_fill_manual(values = wes_palette("Royal1", n = 2)) +
# geom_histogram(aes(y=..density..), colour="black", fill="white")+
 geom_density(alpha=.6, size=1, color="black" ) + theme_bw()
pdf_and_png(p1, paste0(path_output_fb, "Ppargc1a_wt_cor_0.7"), h_to_w_ratio = 0.5)

exp <- AverageExpression(object = fb12, features = ppargc1a.wt.genes, group.by = "KO")
ppargc1a.wt.sig.exp <- data.frame(exp$RNA)
ppargc1a.wt.sig.exp$gene <- rownames(ppargc1a.wt.sig.exp)
ppargc1a.wt.sig.exp.df <- melt(ppargc1a.wt.sig.exp, id.vars = "gene", measure.vars = c("KO", "WT"))
colnames(ppargc1a.wt.sig.exp.df) <- c("gene", "condition", "average.exp")

ppargc1a.wt.sig.exp.df$logExp <- log2(ppargc1a.wt.sig.exp.df$average.exp + 1)
  
p2 <- ggplot(ppargc1a.wt.sig.exp.df, aes(x=condition, y = logExp, fill=condition)) +
   scale_fill_manual(values = wes_palette("Royal1", n = 2)) +
# geom_histogram(aes(y=..density..), colour="black", fill="white")+
 geom_violin(alpha=.6, size=1, color="black" ) + theme_bw()
pdf_and_png(p2, paste0(path_output_fb, "Ppargc1a_exp_cor_0"), h_to_w_ratio = 0.5)


################################################
load(file="/path_to_bulk_TPM_dataset/fibroblast_unprocessed_alldata.RData")
TPM <- fibroblast_data$expr$TPM
GRannot <- fibroblast_data$gene_annot$DF

genes_annot <- as.data.frame(GRannot)
bicor.wt.ppargc1a.genes.ensemID <- 
  genes_annot[genes_annot$external_gene_name %in% ppargc1a.wt.genes,"gene_id"]

tpm.ppargc1a <- TPM[bicor.wt.ppargc1a.genes.ensemID,]
tpm.ppargc1a.df <- melt(tpm.ppargc1a)
colnames(tpm.ppargc1a.df) <- c("EnsemID", "Sample", "TPM")
tpm.ppargc1a.df$Conditon <- rep(c("HO","HO_T","WT","WT_T"), each=95)
tpm.ppargc1a.df$Group <-  rep(c("HO","WT"), each=95*2)

tpm.ppargc1a.df$logTPM <- log2(tpm.ppargc1a.df$TPM + 1)

tpm.group.ppargc1a <- tpm.ppargc1a.df %>% group_by(EnsemID,Group) %>%  summarise(mean_TPM = mean(TPM, na.rm = TRUE))
tpm.group.ppargc1a$mean_TPM

p <- ggplot(tpm.group.ppargc1a, aes(x=mean_TPM, color=Group, fill=Group)) + 
   scale_fill_manual(values = rev(wes_palette("Royal1", n = 2))) +
# geom_histogram(aes(y=..density..), colour="black", fill="white")+
 geom_density(alpha=.6, size=1, color="black" ) + 
  scale_x_continuous(limits = c(0, 1000)) + 
  theme_bw()

pdf_and_png(p, file_name = paste0(path_output_fb, "col1a1_bicor0.3_p10-4_bulk_density"), width = 10, h_to_w_ratio = 0.4)

################################################
library(clusterProfiler)

resource_path <- "/resource_path/" # path to gmt files
gmtReactome <- read.gmt(paste0(resource_path,"Mouse_Human_Reactome_January_01_2023_symbol.gmt"))
gmtGOBP <- read.gmt(paste0(resource_path,"MOUSE_GO_bp_no_GO_iea_symbol.gmt")）
gmtWiki <- read.gmt(paste0(resource_path,"Mouse_WikiPathways_July_05_2022_symbol.gmt"))
gmtHallmark <- read.gmt(paste0(resource_path,"HALLMARK.Mouse_Human_MSigdb_March_01_2020_symbol.gmt")) 

gsea.reactome.cor <- GSEA(ppargc1a.cor, TERM2GENE = gmtReactome, verbose=T, eps = 0, minGSSize=5, pvalueCutoff = 0.2)
gsea.gobp.cor <- GSEA(ppargc1a.cor, TERM2GENE = gmtGOBP, verbose=T, eps = 0, minGSSize=5, pvalueCutoff = 0.2)

gsea.reactome.cor <- enrichplot::pairwise_termsim(gsea.reactome.cor)
gsea.reactome.cor@result$Description <- sapply(gsea.reactome.cor@result$ID,function(x){unlist(strsplit(x,"%"))[1]})
p1 <- treeplot(gsea.reactome.cor)

saveRDS(gsea.reactome.cor, paste0(path_output, "gsea.reactome.cor.RData"))

ppargc1a.cor.pos5 <- names(ppargc1a.cor[ppargc1a.cor>0.5])
ppargc1a.cor.neg5 <- names(ppargc1a.cor[ppargc1a.cor< -0.5])


####
ort.reactome.pos5 <- enricher(ppargc1a.cor.pos5, pvalueCutoff = 0.01,TERM2GENE = gmtReactome)
ort.reactome.pos5 <- enrichplot::pairwise_termsim(ort.reactome.pos5)
colnames(ort.reactome.pos5@termsim) <- rownames(ort.reactome.pos5@termsim) <- sapply(colnames(ort.reactome.pos5@termsim) ,function(x){unlist(strsplit(x,"%"))[1]})
p <- enrichplot::treeplot(ort.reactome.pos5,  cluster.params = list(n = 3))
pdf_and_png(p, file_name = paste0(path_output_fb, "treeplot.ort.reactome.pos5.q0.1"), h_to_w_ratio = 0.4, width = 14)


ort.reactome.pos5@result$Description <- sapply(ort.reactome.pos5@result$ID,function(x){unlist(strsplit(x,"%"))[1]})
ort.reactome.pos5@result$ID <- sapply(ort.reactome.pos5@result$ID,function(x){unlist(strsplit(x,"%"))[3]})

ort.reactome.pos5.output <- ort.reactome.pos5@result %>% filter(p.adjust<0.01)

write.xlsx(ort.reactome.pos5.output, file = paste0(path_output_fb, "top.reactome.output.q0.1.xlsx"), row.names = F, col.names = T, sheetName = "postive.cor.0.5", append = T)

###
ort.gobp.pos5 <- enricher(ppargc1a.cor.pos5, pvalueCutoff = 0.05,TERM2GENE = gmtGOBP)
ort.gobp.pos5@result$Description <- sapply(ort.gobp.pos5@result$ID,function(x){unlist(strsplit(x,"%"))[1]})
ort.gobp.pos5@result$ID <- sapply(ort.gobp.pos5@result$ID,function(x){unlist(strsplit(x,"%"))[1]})

ort.gobp.pos5 <- enrichplot::pairwise_termsim(ort.gobp.pos5)
p <- enrichplot::treeplot(ort.gobp.pos5,  cluster.params = list(n = 3))

###
ort.hallmark.pos5 <- enricher(ppargc1a.cor.pos5, pvalueCutoff = 0.2,TERM2GENE = gmtHallmark)
ort.hallmark.pos5@result$Description <- sapply(ort.hallmark.pos5@result$ID,function(x){unlist(strsplit(x,"%"))[1]})
ort.hallmark.pos5@result$ID <- sapply(ort.hallmark.pos5@result$ID,function(x){unlist(strsplit(x,"%"))[1]})

ort.hallmark.pos5 <- enrichplot::pairwise_termsim(ort.hallmark.pos5)
p <- enrichplot::treeplot(ort.hallmark.pos5,  cluster.params = list(n = 3))


####
ort.reactome.neg5 <- enricher(ppargc1a.cor.neg5, pvalueCutoff = 0.01,TERM2GENE = gmtReactome)
ort.reactome.neg5 <- enrichplot::pairwise_termsim(ort.reactome.neg5)
colnames(ort.reactome.neg5@termsim) <- rownames(ort.reactome.neg5@termsim) <- sapply(colnames(ort.reactome.neg5@termsim) ,function(x){unlist(strsplit(x,"%"))[1]})
p <- enrichplot::treeplot(ort.reactome.neg5,  cluster.params = list(n = 3))
pdf_and_png(p, file_name = paste0(path_output_fb, "treeplot.ort.reactome.neg5.q0.1"), h_to_w_ratio = 0.8, width = 14)

ort.reactome.neg5@result$Description <- sapply(ort.reactome.neg5@result$ID,function(x){unlist(strsplit(x,"%"))[1]})
ort.reactome.neg5@result$ID <- sapply(ort.reactome.neg5@result$ID,function(x){unlist(strsplit(x,"%"))[3]})

ort.reactome.neg5.output <- ort.reactome.neg5@result %>% filter(p.adjust<0.01)

write.xlsx(ort.reactome.neg5.output, file = paste0(path_output_fb, "top.reactome.output.q0.1.xlsx"), row.names = F, col.names = T, sheetName = "negative.cor.0.5", append = T)

ort.reactome.neg5@result$Description <- sapply(ort.reactome.neg5@result$ID,function(x){unlist(strsplit(x,"%"))[1]})
ort.reactome.neg5@result$ID <- sapply(ort.reactome.neg5@result$ID,function(x){unlist(strsplit(x,"%"))[1]})



correlation of the targeted pathways (genes from pathways tested significant by ORT Reactome)
pos.genes1 <- unique(strsplit(paste(ort.reactome.pos5.output$geneID[c(1,3,5,6)], collapse = "/"), split = "/")[[1]])
pos.genes2 <- unique(strsplit(paste(ort.reactome.pos5.output$geneID[c(2,7,8,9)], collapse = "/"), split = "/")[[1]])
pos.genes3 <- unique(strsplit(paste(ort.reactome.pos5.output$geneID[c(4,10,11)], collapse = "/"), split = "/")[[1]])
neg.genes <- unique(strsplit(paste(ort.reactome.neg5.output$geneID, collapse = "/"), split = "/")[[1]])

pos.genes.network <- correlations.all[pos.genes,pos.genes]
p <- pheatmap(pos.genes.network)
cluster = cutree(p$tree_row, k = 2)
names(cluster)[cluster==1]


