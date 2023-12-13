#!/usr/local/bin/Rscript
# R 4.1.1
#---------------------------------------------------------------------------------
# Mouse Marginal zone B-cell with/without Tfh. Ensemble GRCm39
# CAD_mnus421_0002
# Link to publication
# TO ADD ONCE AVAILABLE
#
# Script available from:
# https://github.com/CAD-BFX/CAD_mnus421_0002
#
#
# Analysis Performed by Xiaohui Zhao
# Department of medicine, University of Cambridge, UK
# Copyright Xiaohui Zhao (xz289@cam.ac.uk)
#
#---------------------------------------------------------------------------------
# License Information
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------------

message("+--------------------------------------------------------------------+")
message("+------        Basic settings and libraries           ---------------+")
message("+--------------------------------------------------------------------+")

suppressPackageStartupMessages({
  library("tidyr")
  library("dplyr")
  library("methods")
  library("utils")
  library("ggplot2")
  library("ggrepel")
  library("cowplot")
  library("Matrix")
  library("matrixStats")
  library("useful")
  library("reshape")
  library("reshape2")
  library("DESeq2")
  library("biomaRt")
  library("ggforce")
  library("pheatmap")
  library('RColorBrewer')
  library("scales")
  library("ggbeeswarm")
  library("BiocParallel")
  library("ggalt")
  library("ComplexHeatmap")
  library("apeglm")
  library("openxlsx")
})

NUMCORES      <- 3
register(MulticoreParam(NUMCORES))

Project <-"CAD_mnus421_0002"
baseDir <- "/Users/xz289/Documents/Meri/CAD_mnus421_0002"
setwd(baseDir)

TOPNUM       <- 2000
l2fc         <- 1
significance <- 0.05
elementTextSize <- 6

message("+-------------------------------------------------------------------------------+")
message("+                               Use ensEMBL Annotations                         +")
message("+-------------------------------------------------------------------------------+")

ensembl    =  useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
listEnsembl()
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description','chromosome_name'),  mart = ensembl, useCache = FALSE) 
ensEMBL2id$description <- gsub("..Source.*", "", ensEMBL2id$description)

head(ensEMBL2id)
nrow(ensEMBL2id) ## 55416

save(ensEMBL2id, file = "./DESeq_Analysis/Ensembl_mmusculus_ID_Name_Des_Chr_GRCm39.RData")

message("+-------------------------------------------------------------------------------+")
message("+    DESeq2 Analaysis---->PCA plot                                              +")
message("+-------------------------------------------------------------------------------+")

load("./DESeq_Analysis/deseq2.dds.RData")  ## 55359
cts                   <- assay(dds)
samT                  <- read.csv("./DESeq_Analysis/CAD_mnus421_0002_SampleTable.csv", header=T)
samT$sample           <- paste0(samT$Condition, "_", samT$Replicates)
samT                  <- samT[order(samT$sample),]
## set control group WT level
samT$Condition        <- as.factor(samT$Condition )
samT$Condition        <- relevel(samT$Condition , "WT")

message("+----          DESeq2 model with batch sex and condition in the design                ----+")

dds.cond.batch      <- DESeqDataSetFromMatrix(countData = cts,
                                              colData = samT,
                                              design= ~ Sex+Condition)
dds.cond.batch      <- estimateSizeFactors(dds.cond.batch)
sizeFactors(dds.cond.batch)

dds.cond.batch      <- DESeq(dds.cond.batch, parallel=TRUE)
vsd.cond.batch      <- vst(dds.cond.batch,     blind=F)

mat                 <- assay(vsd.cond.batch)
mat                 <- limma::removeBatchEffect(mat, vsd.cond.batch$Sex)
assay(vsd.cond.batch)  <- mat
counts_batch_corrected <- assay(vsd.cond.batch)
resultsNames(dds.cond.batch)
colData(vsd.cond.batch)

message("+----                     DESeq2 PCA plot                                   -------+")

customPCA     <- function(sampleTBL, RLD, TOPNUM, model, ensEMBL2id, matcols) {
  RLD                 <- as.data.frame(RLD)             
  RLD$ensembl_gene_id <- rownames(RLD)
  RLD.mer             <- merge(RLD, ensEMBL2id, by="ensembl_gene_id") 
  RLD.mer             <- RLD.mer[-which(RLD.mer$external_gene_name==""),]## remove no gene names
  RLD.mer             <- RLD.mer[-(which(duplicated(RLD.mer$external_gene_name)==T)),] # remove duplicated genes
  RLD.new             <- RLD.mer[,matcols]
  rownames(RLD.new)   <- RLD.mer$external_gene_name
  rv     <- rowVars(as.matrix(RLD.new))
  select <- order(rv, decreasing = TRUE)[seq_len(min(TOPNUM, length(rv)))]
  pca    <- prcomp(t(RLD.new[select, ]))
  
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
  
  scores    <- data.frame(sampleName=sampleTBL$sample, pca$x, 
                          condition=sampleTBL$Condition, sex=sampleTBL$Sex)
  
  
  plt.pca.new <- ggplot(scores, aes(x = PC1, y = PC2, colour=condition, fill=condition, shape=sex, label=sampleName) ) +
    geom_text_repel(aes(label=sampleName), show.legend = FALSE, size=4) +
    geom_point(size=4, alpha=0.75) +
    scale_color_manual(values=c("blue", "red")) +
    scale_shape_manual(values = c(15,16)) +
    xlab(pc1lab) +
    ylab(pc2lab) + 
    coord_fixed() +
    theme_bw() +
    theme(text = element_text(size=12), 
          plot.subtitle=element_text(size=8,  face="italic", color="darkgrey"),
          legend.title = element_blank(),
          legend.key = element_rect(fill = "white",colour = "white"),
          legend.position='right', 
          aspect.ratio=1,
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"))
  
  loadings                 <- as.data.frame(pca$rotation)
  
  
  pca.1         <- loadings[order(loadings$PC1,decreasing=TRUE), ]
  pca.1.25      <- pca.1[c(1:25),]
  pca.1.25.plot <- ggplot(data=pca.1.25, aes(x=factor(rownames(pca.1.25),levels=unique(rownames(pca.1.25))), y=PC1)) + 
    geom_point(size = 3 ) + xlab("") + 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, size=8))
  
  pca.2         <- loadings[order(loadings$PC2,decreasing=TRUE), ]
  pca.2.25      <- pca.2[c(1:25),]
  pca.2.25.plot <- ggplot(data=pca.2.25, aes(x=factor(rownames(pca.2.25),levels=unique(rownames(pca.2.25))), y=PC2)) + 
    geom_point(size = 3 ) + xlab("") + 
    theme()+ 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, size=8))
  
  return(list(plt.pca.new,pca.1.25.plot, pca.2.25.plot))
  
}



sampleTBL.batch      <- as.data.frame(colData(dds.cond.batch))
RLD.con.batch        <- assay(vsd.cond.batch)
models               <- c("vsd.cond.batch") 

pca.con.batch       <- customPCA(sampleTBL, RLD.con.batch, TOPNUM, models[1], ensEMBL2id, matcols)

pdf("./DESeq_Analysis/Fig1_PCA_batch_sex_nofilter.pdf")
print(pca.con.batch[[1]])
dev.off()
pdf("./DESeq_Analysis/Fig1_PCA_nofilter_batch_sex_PC1_PC2.pdf", width=7, height= 10)
plot_grid(pca.con.batch[[2]],pca.con.batch[[3]], nrow=2)
dev.off()

message("+---- Housekeep genes individual Plot checking the outliers         -----------+")

res.cond.batch      <- lfcShrink(dds=dds.cond.batch, coef="Condition_NT_vs_WT", type="apeglm", parallel=TRUE)

nrow(subset(res.cond.batch, padj <= significance & abs(log2FoldChange) >= l2fc) ) ## 17
head(subset(res.cond.batch, padj <= significance & abs(log2FoldChange) >= l2fc) )
res.cond.batch.dat  <- as.data.frame(res.cond.batch)
res.cond.batch.dat$ensembl_gene_id <- rownames(res.cond.batch.dat)
res.cond.batch.merE <- merge(res.cond.batch.dat, ensEMBL2id, by = "ensembl_gene_id")
res.cond.batch.merE <- subset(res.cond.batch.merE, external_gene_name!="") ## 52324
res.cond.batch.merE <- res.cond.batch.merE[-which(duplicated(res.cond.batch.merE$external_gene_name)==T),] ## 52324
res.cond.batch.merE.rm <- subset(res.cond.batch.merE, !is.na(padj)) ## 18010

sum(res.cond.batch.merE.rm$padj<0.05 & abs(res.cond.batch.merE.rm$log2FoldChange)>=1)
normCounts.batch    <- as.data.frame(assay(vsd.cond.batch))
normCounts.batch$ensembl_gene_id <- rownames(normCounts.batch)

mzbmarkers  <- c("Notch2", "CD21", "Pdl1", "B220", "Cd1d","36B4")
res.cond.batch.merEC.rm <- merge(res.cond.batch.merE.rm, normCounts.batch, by = "ensembl_gene_id")
write.csv(res.cond.batch.merEC.rm , file = "./DESeq_Analysis/CAD_mnus421_0002-Cond_batchrm_Res_summaryTable.csv", row.names=F)
test <- subset(res.cond.batch.merE.rm, padj<0.05& abs(log2FoldChange)>=1)

save(res.cond.batch, file = "./DESeq_Analysis/Cond_batchrm_dds_res.RData")

## DEGs of gender
res.cond.batch.sex      <- lfcShrink(dds=dds.cond.batch, coef="Sex_Male_vs_Female", type="apeglm", parallel=TRUE)

nrow(subset(res.cond.batch.sex, padj <= significance & abs(log2FoldChange) >= l2fc) ) ## 865
head(subset(res.cond.batch.sex, padj <= significance & abs(log2FoldChange) >= l2fc) )
res.cond.batch.sex.dat  <- as.data.frame(res.cond.batch.sex) ## 55359
res.cond.batch.sex.dat$ensembl_gene_id <- rownames(res.cond.batch.sex.dat)
res.cond.batch.sex.merE <- merge(res.cond.batch.sex.dat, ensEMBL2id, by = "ensembl_gene_id")
res.cond.batch.sex.merE <- subset(res.cond.batch.sex.merE, external_gene_name!="") 
res.cond.batch.sex.merE <- res.cond.batch.sex.merE[-which(duplicated(res.cond.batch.sex.merE$external_gene_name)==T),] ## 52324
res.cond.batch.sex.merE.rm <- subset(res.cond.batch.sex.merE, !is.na(padj)) ## 14265
test.sex <- subset(res.cond.batch.sex.merE.rm, padj<0.05& abs(log2FoldChange)>=1)

message("+---clusterProfiler to perform the GSEA analysis and pathway analysis, with pval < 0.05, regardless fold changes----------+")

library(clusterProfiler)
library(ReactomePA)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
load("/Users/xz289/Documents/Minoru/Project_mt709_0001/Data/References_Data/GRCm39/Ensembl_GRCm39_entrezID_ensembID_exterName.RData") ##
nrow(ensEMBL2id.GO)
degs  <- read.csv("./DESeq_Analysis/CAD_mnus421_0002-Cond_batchrm_Res_summaryTable.csv", header = T)
degsp <- subset(degs, pvalue <= 0.05) ## 1058
upgenes <- subset(degsp, log2FoldChange > 0) ## 515
dwgenes <- subset(degsp, log2FoldChange < 0) ## 543

## prepare genes input for clusterProfiler

GO_Data_conversion <- function(resdat, ensGO, pcut=0.05, selCol.names, outfile){
  resdat <- subset(resdat, !is.na(padj) & !is.na(log2FoldChange))
  resann <- merge(resdat, ensGO, by = "external_gene_name", all.x = T)
  resann.sig <- subset(resann, pvalue <= pcut)
  resdf <- resann.sig[, selCol.names]
  colnames(resdf)  <- c("ENTREZID", "ENSEMBL", "SYMBOL", "L2FC")
  resdf <- resdf[order(-resdf$L2FC),]
  write.csv(resdf, file = outfile, row.names = F)
  gc()
}
ensGO <- ensEMBL2id.GO
sel.column.names <-c("entrezgene_id", "ensembl_gene_id.x", "external_gene_name", "log2FoldChange")
outfile <- "./DESeq_Analysis/ClusterProfiler_resdf_pvalue05_review_03_10_2023.csv"
goinput <- GO_Data_conversion(degs, ensGO, pcut=0.05, selCol.names=sel.column.names, outfile=outfile)

## approach 1, up and down regulated genes all together. Only do BP and KEGG, plus gsea to quantify.
## gene 
## KEGG pathways are not stable, rerun later, 09/03/2032
GO_enrich_fn <- function(input,pval=0.05, ontology, gsemin=10, gsemax=500, outfile){
  ## input is the data, inputInd is the sheetIndex for input data, bkInd is the index for background.
  ## pval=0.05,  ontology (BP, CC, MF or All)
  ## gsemin= 10, gsemax = 500, outfile(saved R object)
  set.seed(1234)
  resdf       <- read.csv(input, header = T)
  geneDat     <- subset(resdf, !is.na(resdf$ENTREZID))
  geneList    <- geneDat$L2FC
  names(geneList) <- geneDat$ENTREZID
  
  enrichGO_Res <- enrichGO(gene = as.character(na.omit(resdf$ENTREZID)), 
                           keyType = "ENTREZID", 
                           OrgDb = org.Mm.eg.db,
                           ont = ontology, 
                           pAdjustMethod = "BH", 
                           pvalueCutoff  = pval,
                           readable      = TRUE)
  
  enrichGSE_Res  <- gseGO(geneList = geneList, OrgDb = org.Mm.eg.db,
                          ont = "BP", minGSSize = 10, maxGSSize = 500, eps=0,
                          pvalueCutoff = 0.05, verbose = FALSE, nPermSimple = 10000)
  if(dim(as.data.frame(enrichGSE_Res))[1]!=0){
    enrichGSE_Res  <- setReadable(enrichGSE_Res, OrgDb = org.Mm.eg.db, keyType="ENTREZID")}
  
  enrichKEGG_Res <-  enrichKEGG(gene = names(geneList),  organism = 'mmu', pvalueCutoff =0.05)
  if(dim(as.data.frame(enrichKEGG_Res))[1]!=0){
    enrichKEGG_Res <- setReadable(enrichKEGG_Res, OrgDb = org.Mm.eg.db, keyType="ENTREZID")}
  
  gseKEGG_Res <- gseKEGG(geneList, organism = "mmu", keyType = "kegg",exponent = 1,
                         minGSSize = 10, maxGSSize = 500, eps = 0, pvalueCutoff = 0.05)
  if(dim(as.data.frame(gseKEGG_Res))[1]!=0){
    gseKEGG_Res <- setReadable(gseKEGG_Res, OrgDb = org.Mm.eg.db, keyType="ENTREZID")}
  
  save(enrichGO_Res, enrichGSE_Res, enrichKEGG_Res, gseKEGG_Res, file = outfile)
  
  gc()
}
gooutfiles   <- paste0("./DESeq_Analysis/GO_Analysis_pval05_all_NTvsWT_BP_KEGG-Oct_2023.RData")

GOallAnalysis <- GO_enrich_fn(outfile, pval=0.05,  "All",  gsemin=10, gsemax=500, gooutfiles)

load(gooutfiles)
goenrich <- as.data.frame(enrichGO_Res)
kenrich <- as.data.frame(enrichKEGG_Res)
kenrich$ONTOLOGY <- "KEGG"
kenrich <- kenrich[,c(10, 1:9)]
allenrich <- rbind(goenrich, kenrich)

write.csv(allenrich, file = "./DESeq_Analysis/EnrichGO_KEGG_pval05_list_Oct_2023.csv", row.names=F)

message("+------------enrich comparison and reactome only with up/down -------------------------+")

## prepare data for comparison  up/dw data frame
mydf <- read.csv(outfile, header = T); 
mydf <- unique(subset(mydf, !is.na(ENTREZID)))
mydf$group <- ifelse(mydf$L2FC > 0, "upregulated", "downregulated")
## the final analysis may have slightly number difference due to the missingness of geneID,Name or others.

write.csv(mydf, file = "./DESeq_Analysis/Input_data_compareCluster_updw_pval05_list_Oct_2023.csv", row.names=F)

## Data conversion for reactome from mouse to human orthology
## get the cmydf geneSymbol to gprofiler2 orthology and convert to human gene name
orthodata <- read.csv("./DESeq_Analysis/gProfiler_mmusculus_hsapiens_03-10-2023_20-56-56.csv", header = T)
colnames(orthodata)[3] <- "ENSEMBL"
orthomydf <- merge(mydf, orthodata, by ="ENSEMBL")
table(orthomydf$group)
#downregulated   upregulated 
#467           427 

orthomydf$EntrezH <- mapIds(org.Hs.eg.db, keys=orthomydf$ortholog_name, column = "ENTREZID", keytype = "SYMBOL")
orthomydf <- orthomydf[order(-orthomydf$L2FC), ]

## enrichGO compare
formula_res1 <- compareCluster(ENTREZID~group, data=mydf, fun="enrichGO", OrgDb='org.Mm.eg.db')
formula_res1 <- setReadable(formula_res1, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
write.csv(as.data.frame(formula_res1), file = "./DESeq_Analysis/EnrichGO_updw_compare_N22_03_10_2023.csv", row.names = F)

## enrich reactome, Nothing
upgList <- orthomydf[orthomydf$L2FC>0, "EntrezH"]
dwgList <- orthomydf[orthomydf$L2FC<0, "EntrezH"]
reactup <- enrichPathway(gene=upgList, pvalueCutoff = 1, pAdjustMethod = "none", readable=TRUE)
reactdw <- enrichPathway(gene=dwgList, pvalueCutoff = 1, pAdjustMethod = "none", readable=TRUE)


message("+------------------------------- individual gene plot -----------------------------------+")

## Individual gene counts plot function
makeGeneCountPlot <- function(DDS, ensEMBL2id, CONDITION, gene2plot, genename) {
  #
  # Plot the log2 normalised read counts for a specified gene
  #
  
  t2            <- plotCounts(DDS, gene=gene2plot, intgroup=c(CONDITION), normalized=TRUE, returnData=TRUE)
  colnames(t2) <- c("count", "condition")
  t2$count <- log2(t2$count)
  
  plt.cont <-  ggplot(t2, aes(x=condition, y=count, fill=condition)) + 
    geom_boxplot(width = 0.5, color=c("red", "blue"), alpha=0.5, outlier.shape=NA) + 
    geom_point(aes(fill=condition),position=position_jitterdodge(),size=0.6, alpha=0.5) +
    scale_fill_manual(name="Condition", values = c("red", "blue")) +
    scale_color_manual(name="Condition", values = c("red", "blue")) +
    theme(text = element_text(size=elementTextSize), legend.position="none") +
    ggtitle(paste0(gene2plot, "-", genename)) + 
    xlab("") + ylab("log2(Normalised count)") +
    theme_classic() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size=elementTextSize,angle = 45, hjust = 1))
  
  
  t2$samples   <- rownames(t2)
  colnames(t2) <- c("count", "condition", "samples")
  t2           <- t2[order(t2$condition),]
  t2$samples2 <- factor(t2$samples, as.character(t2$samples))
  
  plt.ind <- ggplot(t2, aes(x=samples2, y=count, fill=condition, group=condition)) + 
    geom_bar(stat="identity", alpha=0.5) +
    scale_fill_manual(name="Condition", values = c("red", "blue")) +
    xlab("") + ylab("log2(Normalised count)") +
    ggtitle(paste0(gene2plot, "-", genename)) +
    theme_classic() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size=elementTextSize, angle = 45, hjust = 1))
  
  
  print(paste("Created plot for", gene2plot), sep=" ")
  
  return(list(plt.cont, plt.ind))
  
}

housekeepings <- c("Gapdh",  "Ppia", "Ywhaz", "B2m", "Pgk1", "Tbp", "Gusb")
genes2plotd   <- ensEMBL2id[ensEMBL2id$external_gene_name%in%housekeepings,]
ss.sig        <- subset(res.cond, padj <= significance & abs(log2FoldChange) >= l2fc)
ss.sig        <- as.data.frame(ss.sig)
genes2plots   <- ensEMBL2id[ensEMBL2id$ensembl_gene_id%in%rownames(ss.sig),]

pdf("./DESeq_Analysis/Fig_HouseKeep_Gene_Indivi_Group_Plot.pdf", width=14, height = 8)
ind.plt1 <- makeGeneCountPlot(dds.cond, ensEMBL2id, "Condition", genes2plotd$ensembl_gene_id[1], genes2plotd$external_gene_name[1])
ind.plt2 <- makeGeneCountPlot(dds.cond, ensEMBL2id, "Condition", genes2plotd$ensembl_gene_id[2], genes2plotd$external_gene_name[2])
ind.plt3 <- makeGeneCountPlot(dds.cond, ensEMBL2id, "Condition", genes2plotd$ensembl_gene_id[3], genes2plotd$external_gene_name[3])
ind.plt4 <- makeGeneCountPlot(dds.cond, ensEMBL2id, "Condition", genes2plotd$ensembl_gene_id[4], genes2plotd$external_gene_name[4])
ind.plt5 <- makeGeneCountPlot(dds.cond, ensEMBL2id, "Condition", genes2plotd$ensembl_gene_id[5], genes2plotd$external_gene_name[5])
ind.plt6 <- makeGeneCountPlot(dds.cond, ensEMBL2id, "Condition", genes2plotd$ensembl_gene_id[6], genes2plotd$external_gene_name[6])
ind.plt7 <- makeGeneCountPlot(dds.cond, ensEMBL2id, "Condition", genes2plotd$ensembl_gene_id[7], genes2plotd$external_gene_name[7])
ind.plt8 <- makeGeneCountPlot(dds.cond, ensEMBL2id, "Condition", genes2plots$ensembl_gene_id[1], genes2plots$external_gene_name[1])
ind.plt9 <- makeGeneCountPlot(dds.cond, ensEMBL2id, "Condition", genes2plots$ensembl_gene_id[2], genes2plots$external_gene_name[2])
plot_grid(ind.plt1[[1]], ind.plt1[[2]], ind.plt2[[1]], ind.plt2[[2]],ind.plt3[[1]], ind.plt3[[2]],nrow=2, byrow=F)
plot_grid(ind.plt4[[1]], ind.plt4[[2]], ind.plt5[[1]], ind.plt5[[2]],ind.plt6[[1]], ind.plt6[[2]],nrow=2, byrow=F)
plot_grid(ind.plt7[[1]], ind.plt7[[2]], ind.plt8[[1]], ind.plt8[[2]],ind.plt9[[1]], ind.plt9[[2]],nrow=2, byrow=F)
dev.off()

message("+------------------------------- Fig3B, selected DEGs heatmap plot with L2FC -----------------------------------+")
heatmap.dat  <- read.xlsx("Supplementary_Data_1.xlsx", sheet=1)
sigDEGs.dat  <- read.xlsx("Supplementary_Data_1.xlsx", sheet=2)
sigDEGs.datR <- read.xlsx("Supplementary_Data_2.xlsx", sheet=4)

library(dichromat)
summary(sigDEGs.datR$log2FoldChange)
## Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## -1.56823 -0.01415  0.01051  0.05965  0.03493  1.52249
basemeansR           <- as.matrix(sigDEGs.datR[,'log2FoldChange'], ncol=1)
colnames(basemeansR) <- "L2FC"
rownames(basemeansR) <- sigDEGs.datR$external_gene_name
rankBMR   <- basemeansR[order(-basemeansR),]
bk       <- c(seq(-1.6,-0.002,by=0.05)[1:14],seq(0,1.6,0.04)[2:20])
mycols   <- c(colorRampPalette(colors = c("gray28","gray88"))(14),colorRampPalette(colors = c("#FFCE03","#FF6200"))(length(bk)-14)) 
#              heat.colors((length(bk)-14), rev=T))
L2Heat   <- Heatmap(rankBMR, name = "l2FC", col=mycols, width = unit(5, "mm"),
                    show_row_names=T,row_names_gp = gpar(fontsize = 12,  fontface="bold"), cluster_rows = F,
                    row_names_side="right",column_names_gp = gpar(fontsize = 9,  fontface="bold"),
                    heatmap_legend_param = list(at = c(1.522490512,  0.005036496,   -1.568227862), labels = c("Pos1.6", "0.01", "Neg1.6"))) 
pdf("Fig3B-SelMarkers_L2FC_Mean_Heatmap_reduced.pdf", width=3, height= 6)
draw(L2Heat)
dev.off()


message("+------------------------------- Fig3E, selected GESA pathways bar plot from EnrichR ---------------------------+")

require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
convertHumanGeneList <- function(x){
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(humanx)
  return(humanx)
}
pgenes <- list()
for (i in 1:nrow(pathway.dat)){
  subgenes <- convertHumanGeneList(unlist(strsplit(pathway.dat[i, "Genes"] , split=";")))
  pgenes[[i]] <- subgenes
}
library(reshape2) ## melt function
BP_list_up  <- list();
BP_list_tot <- list()
for (i in 1:nrow(pathway.dat)){
  df_tmp  <- sigDEGs.dat[sigDEGs.dat$external_gene_name %in% pgenes[[i]],]
  tmp_up  <- length(subset(df_tmp[,3], df_tmp[,3] > 0))
  tmp_tot <- dim(df_tmp)[1]
  BP_list_up[[i]]  <- tmp_up
  BP_list_tot[[i]] <- tmp_tot
}

pathway.dat$UP    <- as.numeric(unlist(BP_list_up))
pathway.dat$DOWN  <- as.numeric(unlist(BP_list_tot)) - pathway.dat$UP
pathway.dat$DOWN  <- -pathway.dat$DOWN
pathway.dat$Count <- as.numeric(unlist(BP_list_tot))
pathway.dat_molten <- melt(pathway.dat[,c(1:2,9:11)],
                           id.vars=c("Term","Count", "P.value") )
pathway.dat_molten$P.value <- as.numeric(pathway.dat_molten$P.value)
pathway.dat_molten$Tp.value<- -log10(pathway.dat_molten$P.value)
pathway.dat_molten$Tp.value<- round(pathway.dat_molten$Tp.value, digits=4)
pathway.dat_molten$P.value<- round(pathway.dat_molten$P.value, digits=4)
test  <- unlist(lapply(pathway.dat_molten$Term, function(x) gsub("Homo sapiens R-HSA-.*", "", x)))
test1 <- unlist(lapply(test, function(x) gsub("Homo sapiens P.*", "", x)))
pathway.dat_molten$Term <- test1
brks <- seq(-15, 15, 5)
lbls = as.character(c(seq(15, 0, -5), seq(5, 15, 5)))


pdf("Fig3E-EnrichR_selN21_pathways_02_Nov_2021.pdf", width=12, height=8)
barplotSum <- ggplot(pathway.dat_molten, aes(x = reorder(Term, value), y = value, fill = variable)) +   # Fill column
  geom_bar(stat = "identity", width = 0.8) +   # draw the bars
  scale_y_continuous(breaks = brks,   # Breaks
                     labels = lbls,
                     limits=c(-15,15)) + # Labels
  coord_flip() +  
  xlab("") +
  ylab("") +
  guides(fill=FALSE) +
  theme(plot.title = element_text(hjust = .2),
        axis.ticks = element_blank()
  ) +   
  scale_fill_manual(values = c("orange", "gray"))+  # Color palette
  theme_bw() +
  theme_update(axis.title.x = element_text(size=12, face= "bold"),
               axis.text.x = element_text(size=12, face="bold"),
               axis.title.y = element_text(size=12, face= "bold"),
               axis.text.y.left = element_text(size=16, face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),panel.background = element_rect(fill = "white", colour = NA)) 


barplotSum

dev.off()
## ----------------Finish, final editing 13/12/2023--------------------------------##











##------------FINISH---------------------------------------------------------##
0