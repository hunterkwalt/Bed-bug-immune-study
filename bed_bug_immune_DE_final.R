library(tximport)
library(DESeq2)
library(readr)
library("tximportData")
library(GenomicFeatures)
library(ggpubr)
library(ggrepel)

################ DESeq Analysis ################################################

#path to quant files
#directory = "/mnt/md0/cimex_infection/quants/"
directory = "/mnt/md0/cimex_infection/bb_infection_sophie"

#quantfiles list: these are just a list of the sample names
TI <- read.table(file = "/mnt/md0/cimex_infection/TI_samps.txt", header = T)
EN <- read.table(file = "/mnt/md0/cimex_infection/EN_samps.txt", header = T)
PT <- read.table(file = "/mnt/md0/cimex_infection/PT_samps.txt", header = T)

colnames(TI) <- gsub("_quants", "", colnames(TI))
colnames(EN) <- gsub("_quants", "", colnames(EN))
colnames(PT) <- gsub("_quants", "", colnames(PT))

#full list of sample names
allsamps <- read.table(file = "/mnt/md0/cimex_infection/samps.txt", header = F)
colnames(allsamps) <- "samples"

#subset this list.
TI_names <- allsamps$samples[13:18]
EN_names <- allsamps$samples[1:6]
PT_names <- allsamps$samples[7:12]

#make file paths s R can find the files
TI_files <- file.path(directory, colnames(TI), "quant.sf")
names(TI_files) <- paste0(TI_names)

EN_files <- file.path(directory, colnames(EN), "quant.sf")
names(EN_files) <- paste0(EN_names)

PT_files <- file.path(directory, colnames(PT), "quant.sf")
names(PT_files) <- paste0(PT_names)

#annotation file to make txdb
gff = "/mnt/md0/cimex_infection/bb_annotation.gff.gz"

#to save txdb
txdb.filename = "/mnt/md0/cimex_infection/bb_annotation.sqlite"

#make and save txdb 
txdb <- makeTxDbFromGFF(file = gff)
saveDb(txdb, txdb.filename)

#maketx2gene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
tx2gene <- tx2gene[-c(1),]

#savetx2gene
write.csv(tx2gene, "/mnt/md0/cimex_infection/bb_tx2gene.csv", row.names = F)
tx2gene <- read_csv("/mnt/md0/cimex_infection/bb_tx2gene.csv")
head(tx2gene)

#make object TI
txi_TI <- tximport(TI_files, type = "salmon", tx2gene = tx2gene)
names(txi_TI)

head(txi_TI$counts)

#make object EN
txi_EN <- tximport(EN_files, type = "salmon", tx2gene = tx2gene)
names(txi_EN)

head(txi_EN$counts)

#make object PT
txi_PT <- tximport(PT_files, type = "salmon", tx2gene = tx2gene)
names(txi_PT)

head(txi_PT$counts)

#import EN to DEseq
EN_table <- data.frame(condition = factor(rep( c("CT", "PF"), each = 3)))
rownames(EN_table) <- colnames(txi_EN$counts)
ENdds <- DESeqDataSetFromTximport(txi_EN, EN_table, ~condition)

#import TI to DEseq
TI_table <- data.frame(condition = factor(rep( c("BC", "NV"), each = 3)))
rownames(TI_table) <- colnames(txi_TI$counts)
TIdds <- DESeqDataSetFromTximport(txi_TI, TI_table, ~condition)

#import PT to DEseq
PT_table <- data.frame(condition = factor(rep( c("BD", "CT"), each = 3)))
rownames(PT_table) <- colnames(txi_PT$counts)
PTdds <- DESeqDataSetFromTximport(txi_PT, PT_table, ~condition)


#relevel factors - it does alhpabetical by default and I want my treatment group to be the positive numbers for DE genes
TIdds$condition <- factor(TIdds$condition, levels = c("NV", "BC"))
PTdds$condition <- factor(PTdds$condition, levels = c("CT", "BD")) 

#get results
ENdds <- DESeq(ENdds)
TIdds <- DESeq(TIdds)
PTdds <- DESeq(PTdds)

ENresults <- results(ENdds, alpha = 0.05)
TIresults <- results(TIdds, alpha = 0.05)
PTresults <- results(PTdds, alpha = 0.05)

#summarize results
summary(ENresults)
summary(TIresults)
summary(PTresults)

######################################### QC Plots ##################################################

library(pheatmap)

### TI PCA and Cluster
vsd <- vst(TIdds, blind=FALSE)
TIpca <- plotPCA(vsd, intgroup = "condition")
TIpca <- TIpca + ggtitle("TI Bacteria (Bacillus spp.)")

select <- order(rowMeans(counts(TIdds,normalized=TRUE)),
                decreasing=TRUE)[1:1000]

df <- as.data.frame(colData(TIdds))
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=T,
         cluster_cols=T, annotation_col=df)


### EN PCA and Cluster
vsd <- vst(ENdds, blind=FALSE)
ENpca <- plotPCA(vsd, intgroup = "condition")
ENpca <- ENpca + ggtitle("Pseudomonas fluorescens")

select <- order(rowMeans(counts(ENdds,normalized=TRUE)),
                decreasing=TRUE)[1:1000]

df <- as.data.frame(colData(ENdds))
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=T,
         cluster_cols=T, annotation_col=df)


### PT PCA and Cluster
vsd <- vst(PTdds, blind=FALSE)
PTpca <- plotPCA(vsd, intgroup = "condition")
PTpca <- PTpca + ggtitle("Borrelia duttoni")

select <- order(rowMeans(counts(PTdds,normalized=TRUE)),
                decreasing=TRUE)[1:100]

df <- as.data.frame(colData(PTdds))
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=T,
         cluster_cols=T, annotation_col=df)

#export pcas for supp data
#ggsave(ENpca, filename = "/mnt/md0/cimex_infection/figures/updated/EN_pca.svg")
#ggsave(TIpca, filename = "/mnt/md0/cimex_infection/figures/updated/TI_pca.svg")
#ggsave(PTpca, filename = "/mnt/md0/cimex_infection/figures/updated/PT_pca.svg")


######################## Make volcano Plots #################

#full immune genes - this is a list of all bed bug immune genes. from Benoit 2016 and Meraj 2024
#these will be used to idenify differentially expressed immune-related genes
all_immune_genes <- read.table("/mnt/md0/cimex_infection/meraj_data/immune_genes.txt", header = F)
all_immune_genes <- all_immune_genes$V1

#make volcano plot for P. fluorescens. Upregulated genes in Red, downregulated in blue, immune genes labeled in Yellow. insignif = gray
EN_sig_up <- subset(ENresults, padj<= 0.05 & log2FoldChange >= 1)
EN_sig_down <- subset(ENresults, padj<= 0.05 & log2FoldChange <= -1 )
EN_sig <- na.omit(subset(ENresults, padj <= 0.05 & (log2FoldChange >= 1 | log2FoldChange <= -1)))
EN_immune_de <- EN_sig[rownames(EN_sig) %in% all_immune_genes,]
EN_volcano = ggplot(data=data.frame(ENresults), aes(x=log2FoldChange, y=-log10(padj))) + geom_point(size = 2.5, color = "grey") + 
  geom_point(data = data.frame(EN_sig_up), aes(x=log2FoldChange, y=-log10(padj)), color = "red4",size = 2.5) + ggtitle("P. fluorescens Infected vs. Uninfected")+
  geom_point(data = data.frame(EN_sig_down), aes(x=log2FoldChange, y=-log10(padj)), color ="blue",size = 2.5) +
  theme(text = element_text(size = 12)) + geom_point(data = data.frame(EN_immune_de),aes(x=log2FoldChange, y=-log10(padj)), color = "goldenrod",size = 2.5 ) +
  geom_label_repel(data = data.frame(EN_immune_de), aes(label = row.names(EN_immune_de)), min.segment.length = 0, box.padding = 0.5)

#make volcano plot for B. duttoni. Upregulated genes in Red, downregulated in blue, immune genes labeled in Yellow. insignif = gray
PT_sig_up <- subset(PTresults, padj<= 0.05 & log2FoldChange >= 1)
PT_sig_down <- subset(PTresults, padj<= 0.05 & log2FoldChange <= -1 )
PT_sig <- na.omit(subset(PTresults, padj <= 0.05 & (log2FoldChange >= 1 | log2FoldChange <= -1)))
PT_immune_de <- PT_sig[rownames(PT_sig) %in% all_immune_genes,]
PT_volcano = ggplot(data=data.frame(PTresults), aes(x=log2FoldChange, y=-log10(padj))) + geom_point(size = 2.5, color = "grey") + 
  geom_point(data = data.frame(PT_sig_up), aes(x=log2FoldChange, y=-log10(padj)), color ="red4",size = 2.5) + 
  geom_point(data = data.frame(PT_sig_down), aes(x=log2FoldChange, y=-log10(padj)), color ="blue",size = 2.5) +
  ggtitle("B. duttoni Infected vs. Uninfected")+
  theme(text = element_text(size = 12), axis.title.y = element_blank())  + 
  geom_point(data = data.frame(PT_immune_de),aes(x=log2FoldChange, y=-log10(padj)), color = "goldenrod",size = 2.5 ) +
  geom_label_repel(data = data.frame(PT_immune_de), aes(label = row.names(PT_immune_de)), min.segment.length = 0, box.padding = 0.5)


#make volcano plot for TI bacteria. Upregulated genes in Red, downregulated in blue, immune genes labeled in Yellow. insignif = gray
TI_sig_up <- subset(TIresults, padj<= 0.05 & log2FoldChange >= 1)
TI_sig_down <- subset(TIresults, padj<= 0.05 & log2FoldChange <= -1 )
TI_sig <- na.omit(subset(TIresults, padj <= 0.05 & (log2FoldChange >= 1 | log2FoldChange <= -1)))
TI_immune_de <- TI_sig[rownames(TI_sig) %in% all_immune_genes,]
TI_volcano = ggplot(data=data.frame(TIresults), aes(x=log2FoldChange, y=-log10(padj))) + geom_point(size = 2.5, color = "grey") + 
  geom_point(data = data.frame(TI_sig_up), aes(x=log2FoldChange, y=-log10(padj)), color = "red4",size = 2.5) + ggtitle("TI Contaminated vs. Sterile")+
  geom_point(data = data.frame(TI_sig_down), aes(x=log2FoldChange, y=-log10(padj)), color ="blue",size = 2.5) + xlim(-10,10) +
  theme(text = element_text(size = 12), axis.title.y= element_blank()) + geom_point(data = data.frame(TI_immune_de),aes(x=log2FoldChange, y=-log10(padj)), color = "goldenrod",size = 2.5 ) +
  geom_label_repel(data = data.frame(TI_immune_de), aes(label = row.names(TI_immune_de)), min.segment.length = 0, box.padding = 0.5)

#the following commands are to save the plots but I ended up manually exporting to get dimensions correct.
#svg("/mnt/md0/cimex_infection/figures/updated/labeled_volcano.svg")
ggarrange(EN_volcano,TI_volcano, PT_volcano, nrow = 1)
#dev.off()


################ Make Immune Gene Tables #################

#get all immune genes from each result
EN_immune_total <- ENresults[rownames(ENresults) %in% all_immune_genes,]
TI_immune_total <- TIresults[rownames(TIresults) %in% all_immune_genes,]
PT_immune_total <- PTresults[rownames(PTresults) %in% all_immune_genes,]


#write out immune genes
#write.table(EN_immune_total, quote = F, col.names = T, sep = "\t", row.names = T,  file = "/mnt/md0/cimex_infection/immune/EN_immune_total.tsv")
#write.table(TI_immune_total, quote = F, col.names = T, sep = "\t", row.names = T,  file = "/mnt/md0/cimex_infection/immune/TI_immune_total.tsv")
#write.table(PT_immune_total, quote = F, col.names = T, sep = "\t", row.names = T,  file = "/mnt/md0/cimex_infection/immune/PT_immune_total.tsv")


################### GO Analysis ##################

library("topGO")

#file with genes mapped to GO terms
geneID2GO <- readMappings(file = "/mnt/md0/cimex_infection/eggnog/gene_level_annotation/GO_gene_level.tsv")

geneNames <- names(geneID2GO)

#function to get genes of interest.
getGeneUniverse <- function(deseqRes){
  vgs = as.numeric(deseqRes$padj)
  names(vgs) = rownames(deseqRes)
  vgs
}

#get padj < 0.05
sel_pval <- function(allScore){ return(allScore < 0.05)}

#P. fluorescens GO
ent_gu <- getGeneUniverse(ENresults)
ent_gu[is.na(ent_gu)] <- 1
GOdata <- new("topGOdata", ontology = "BP", allGenes = ent_gu, geneSel = sel_pval, 
              annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
allGO <- usedGO(object = GOdata)
GOresult <- runTest(GOdata, statistic = "fisher", algorithm = "classic")
allResEN <- GenTable(GOdata, pval = GOresult, topNodes = length(allGO))
knitr::kable(allResEN[1:10,])


#TI GO
ti_gu <- getGeneUniverse(TIresults)
ti_gu[is.na(ti_gu)] <- 1
GOdata <- new("topGOdata", ontology = "BP", allGenes = ti_gu, geneSel = sel_pval, 
              annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
allGO <- usedGO(object = GOdata)
GOresult <- runTest(GOdata, statistic = "fisher", algorithm = "Classic")
allResTI <- GenTable(GOdata, pval = GOresult, topNodes = length(allGO))
knitr::kable(allResTI[1:10,])



# B duttoni GO
bd_gu <- getGeneUniverse(PTresults)
bd_gu[is.na(bd_gu)] <- 1
GOdata <- new("topGOdata", ontology = "BP", allGenes = bd_gu, geneSel = sel_pval, 
              annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
allGO <- usedGO(object = GOdata)
GOresult <- runTest(GOdata, statistic = "fisher", algorithm = "classic")
allResPT <- GenTable(GOdata, pval = GOresult, topNodes = length(allGO))
knitr::kable(allResPT[1:10,])

#write out full GO TSVs
##write.table(subset(allResEN, pval <= 0.05), quote = F, col.names = T, sep = "\t", row.names = F,  file = "/mnt/md0/cimex_infection/GO/p_flour_go.tsv")
##write.table(subset(allResTI, pval <= 0.05), quote = F, col.names = T, sep = "\t", row.names = F, file = "/mnt/md0/cimex_infection/GO/TI_go.tsv")
##write.table(subset(allResPT, pval <= 0.05), quote = F, col.names = T, sep = "\t", row.names = F, file = "/mnt/md0/cimex_infection/GO/b_dut_go.tsv")

################################ VENN DIAGRAM #######################################3

#get differentially expressed Names
ti_venn_l <- row.names(TI_sig)
en_venn_l <- row.names(EN_sig)
pt_venn_l <- row.names(PT_sig)

intersect(ti_venn_l, en_venn_l)
intersect(ti_venn_l, pt_venn_l)
intersect(en_venn_l, pt_venn_l)


#make venn diagram automatically. I will use these values to make a better one manually using eulerr (below).
library(VennDiagram)
venn.diagram(
  x = list(ti_venn_l, pt_venn_l, en_venn_l),
  category.names = c("STI", "BD", "PF"),
  filename = "/mnt/md0/cimex_infection/test_venn_new.png",
  output = T,
  col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3))
)

#use those numbers to make prettier/proportional venn diagram
library(eulerr)

fit1 <- euler(c("P. fluorescens" = 69, "B. duttoni" = 79, "STI" = 23, "P. fluorescens&B. duttoni" = 3, "B. duttoni&STI" = 1,
                "P. fluorescens&STI" = 10, "P. fluorescens&B. duttoni&STI" = 0))

par(cex.main = 0.5)
proportional_venn <- plot(fit1,shape = "ellipse", quantities = list(cex = 1), labels = list(font = 2, cex = 1.7), fill = c("red3", "dodgerblue", "gray"))
svg("/mnt/md0/cimex_infection/figures/updated/updated_venn.svg")
proportional_venn
dev.off()

###################### Check if more genes than expected at random ##################
#followed this tutorial
### https://rdrr.io/bioc/GeneOverlap/man/GeneOverlap.html


ti_en_union <- length(union(row.names(na.omit(ENresults)), row.names(na.omit(TIresults))))
ti_pt_union <- length(union(row.names(na.omit(PTresults)), row.names(na.omit(TIresults))))
pt_en_union <- length(union(row.names(na.omit(ENresults)), row.names(na.omit(PTresults))))

# number of DE genes in each group
en <- nrow(EN_sig)
ti <- nrow(TI_sig)
pt <- nrow(PT_sig)

#intersection of each group
ti_en_int <- length(intersect(ti_venn_l, en_venn_l))
ti_pt_int <- length(intersect(ti_venn_l, pt_venn_l))
en_pt_int <- length(intersect(en_venn_l, pt_venn_l))

#get union of two subsets of genes
ti_en_sub_union <- length(union(ti_venn_l, en_venn_l))
ti_pt_sub_union <- length(union(ti_venn_l, pt_venn_l))
pt_en_sub_union <- length(union(pt_venn_l, en_venn_l))


#test if there are significantly more overlapping genes than one would expect at random.
fisher.test(matrix(c(ti_en_union-ti_en_sub_union,ti - ti_en_int ,en - ti_en_int, ti_en_int),nrow = 2, ncol = 2), alternative = "greater")
fisher.test(matrix(c(ti_pt_union-ti_pt_sub_union,ti - ti_pt_int ,pt - ti_pt_int, ti_pt_int),nrow = 2, ncol = 2), alternative = "greater")
fisher.test(matrix(c(pt_en_union-pt_en_sub_union,pt - en_pt_int ,en - en_pt_int, en_pt_int),nrow = 2, ncol = 2), alternative = "greater")


##### THIS CODE WAS USED TO MAKE SUPPLEMENTARY FIGS 2-4 #######

########## OVERALL IMMUNE GENE EXPRESSION HEATMAPS #########33

##################### Full Set of Immune Genes ##########################################

library(ggplot2)
library(tidyverse)

#The final immune gene list was used to get the transcript names from the bed bug ref.
#These names were used to get the TPM for each transcript.
setwd("/mnt/md0/cimex_infection/meraj_data/")

#this table is the result of that. It is formatted as our treatment, exp/control group, replicate, transcript name, then TPM
orf <- read.table(file = "/mnt/md0/cimex_infection/meraj_data/full_immune_genes_tpm_final.tsv", header = F)

#make header
colnames(orf) <- c("Treatment", "Bacteria","Replicate",  "Transcript", "TPM")

#this is a tx2gene file so that our final result is shown by gene, not by transcript
txs <- read.table("/mnt/md0/cimex_infection/meraj_data/immune_tx_to_gene.tsv", header = F)
colnames(txs) <- c("Transcript", "Gene")

#make a gene column
genes <- rep(txs$Gene, each = 18)
orf$Gene <- genes

#average transcript expression
av <- orf %>% 
  group_by(Treatment, Bacteria, Gene) %>%
  summarise(avg = mean(TPM), .groups = "keep")

#subset data by experiment
en_sub <- subset(orf, Treatment == "EN")
pt_sub <- subset(orf, Treatment == "PT")
ti_sub <- subset(orf, Treatment == "TI")

#make a log2TPM row for visual clairity
en_sub$log2TPM <- log(en_sub$TPM + 1, 2)
pt_sub$log2TPM <- log(pt_sub$TPM + 1, 2)
ti_sub$log2TPM <- log(ti_sub$TPM + 1, 2)


#add column for each sample for heatplot
en_sub$Sample <- paste(en_sub$Treatment, en_sub$Bacteria, en_sub$Replicate, sep = "_")

#make plots

######### P. fluorescens
en_hm <- ggplot(en_sub, aes(Sample, Gene, fill = log2TPM)) + geom_tile() + theme_minimal() + 
  scale_fill_gradient(low = "white", high = "red") + theme(plot.title = element_text(face = "italic"),
                                                           axis.text.y = element_text(size = 4)) + ggtitle("Pseudomonas fluorescens")

########## TI Bacteria
ti_sub$Sample <- paste(ti_sub$Treatment, ti_sub$Bacteria, ti_sub$Replicate, sep = "_")

ti_hm <- ggplot(ti_sub, aes(Sample, Gene, fill = log2TPM)) + geom_tile() + theme_minimal() + 
  scale_fill_gradient(low = "white", high = "red") + theme(axis.text.y = element_text(size = 4 )) + ggtitle("TI bacteria")

############# Borrelia duttoni
pt_sub$Sample <- paste(pt_sub$Treatment, pt_sub$Bacteria, pt_sub$Replicate, sep = "_")

pt_hm <- ggplot(pt_sub, aes(Sample, Gene, fill = log2TPM)) + geom_tile() + theme_minimal() + 
  scale_fill_gradient(low = "white", high = "red") + theme(plot.title = element_text(face = "italic"),
                                                           axis.text.y  = element_text(size = 5)) + ggtitle("Borellia duttoni")

#show all together
ggarrange(en_hm, ti_hm, pt_hm, ncol = 1)

