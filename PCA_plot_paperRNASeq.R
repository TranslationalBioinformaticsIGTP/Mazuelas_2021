####### Packages Needed
library(DESeq2)
library(tximport)
library(org.Hs.eg.db)
library(yaml)
library(ggplot2)
library (ggbeeswarm)
library(pheatmap)
library(viridis)
library(clusterProfiler)
library(DOSE)
library(gage)
library(apeglm)
library(biomaRt)
library(ggrepel)
library(RColorBrewer)
library(gplots)

#######################  Loading Functions #########################
source(file= "./rna_seq_Functions.R")

############################  Parameters  ##########################
#Model of Analysis
model <- "2D"
samples.group <-"Diff.Day"
condition <-samples.group
goterms <- c("BP", "CC", "MF","KEGG")


# Directories
heatmap.dir <- file.path("Results",model, "DEG_heatmap")
if(!file.exists(heatmap.dir)) dir.create(heatmap.dir)

markers.dir <- file.path("Results", model, "markers")
if(!file.exists(markers.dir)) dir.create(markers.dir)

gsea.dir <- file.path("Results", model, "GSEAGO")
if(!file.exists(gsea.dir)) dir.create(gsea.dir)
for(term in goterms){
  if(!file.exists(file.path(gsea.dir, term))) dir.create(file.path(gsea.dir, term))
}
pca.dir <-  file.path("Results", model, "PCAplot")
if(!file.exists(pca.dir)) dir.create(pca.dir)


# Salmon alignement and quantification parameters
file1.suffix <- "_1.fastq.gz"
file2.suffix <- "_2.fastq.gz"
fastqdir <- "/imppc/labs/eslab/mmagallon/Projects/RNA-Seq-timecourse.2/Data"
# transcript.index <- "/imppc/labs/eslab/mmagallon/Projects/RNA-Seq-MPNSTcellLinesVsPNFCellLines/Results/Salmon/salmon_indexes_UCSC_hg38"
output.suffix <- "_quant"
output.quants <- "./Results/Salmon"
# threads <- 8

# Tximport paramenters
orgdb <- org.Hs.eg.db
org.columns <- "SYMBOL"
org.keytype <- "REFSEQ"

# DESeq2 parameters
filt.min.reads <- 5
filt.min.samples <- 1

# color.plate <- viridis(n = 256, alpha = 1, begin = 0, end = 1, option = "D")
color.plate <-  bluered(75)
pvalue <- 0.05
stages <- c("PSC", "NC", "day7", "day14", "day30")
gene.markers <- c() 
for(i in seq_len(length(stages))){
  gr <- stages[i]
  markers<- read.table(file = file.path(markers.dir, paste0("up_", gr,".txt")), header = FALSE,stringsAsFactors = FALSE)[,1]
  names(markers) <- rep(gr, length(markers))
  gene.markers <- c(gene.markers,markers)
  
}

#Enriched terms 2D NF1-/-
cells <- c("3MM", "5MM", "D12")
enrichment.2DNF1 <- list()

for (i in seq_len(length(cells))){
  cell <- cells[i]
  enrichment.2DNF1[[cell]] <- read.table(file = file.path(gsea.dir,"BP", paste0("Enrichment_", cell, "_BPmarkers_up_day30.csv")), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
}


####################################  Loading the data   ######################
sample.data <- read.table(file = "/imppc/labs/eslab/mmagallon/Projects/RNA-Seq-timecourse.2/Sample.Info.AllSamples.3.csv", header = TRUE, sep = "\t", stringsAsFactors = FALSE,comment.char = "")
# sample.data <- sample.data[sample.data$Genotype %in% c("PP"),]
sample.data <- sample.data[which(sample.data$Model == "2D"),]

# sample.data <- sample.data[sample.data$Genotype %in% c("PP","MM_SC"),]
# sample.data <- sample.data[-which(sample.data$Sample.Name == "SC_6PNF"),]

# Getting the file names of the analysis
file.names <- sample.data$File.Name
sample.names <- sample.data$graph.Names

# mycolors <- unique(sample.data$colors.p1[sample.data$Diff.Day %in% c("PSC", "NC", "day7", "day14", "day30")])
# names(mycolors) <- c("PSC", "NC", "day7", "day14", "day30")
# # names(mycolors) <- c("PSC", "NC", "day7", "day14", "day30","SC")

##################              Tximport              #####################
salmonquants.fl <-file.path(output.quants, paste0(file.names,"_quant"), "quant.sf")
names(salmonquants.fl)<- file.names
txi.salmon <- importQuantsData(quant.files = salmonquants.fl, orgdb = orgdb,
                               orgdb.keytype = orgdb.keytype, org.columns = org.columns,
                               tximpot.type = "salmon")
#Change Sample name
colnames(txi.salmon$abundance) <- sample.names
colnames(txi.salmon$counts) <- sample.names
colnames(txi.salmon$length) <- sample.names

# txi.salmon <- selectDataFromTximport(tximport = txi.salmon, sample_names = sample.data$Sample.Name[sample.data$Genotype == "PP"])


##################          Quality Control       ######################

#Normalization
filtered.dds <- getFilteredDDS(tximport = txi.salmon, samples_group = samples.group, samples_df = sample.data)

# dds.deg <- DESeq(object = dds, test = "Wald")
pca.pseudocounts <- rlog(filtered.dds)
# assay(filtered.dds)
# pca.pseudocounts <- pca.pseudocounts[rownames(pca.pseudocounts) %in% gene.markers[!names(gene.markers)%in% "PSC"],]

# cell.type <- sample.data$Spheresvs2D
colData(pca.pseudocounts)$Diff.Day <- factor(as.character(colData(pca.pseudocounts)$Diff.Day), levels=c("PSC","NC","day7","day14","day30","SC"))
Diff.Day <- colData(pca.pseudocounts)$Diff.Day
#WT
Diff.Days <- factor(c("PSC",
                      "NC","NC","NC",
                      "7d","7d","7d",
                      "14d","14d","14d",
                      "30d","30d","30d"), levels = c("PSC","NC","7d","14d","30d"))
colData(pca.pseudocounts)$Diff.Day <-Diff.Days
Diff.Day <- colData(pca.pseudocounts)$Diff.Day
#NF1 and WT
Diff.Days <- factor(c("PSC","PSC","PSC","PSC",
                      "NC","NC","NC","NC","NC","NC",
                      "7d","7d","7d","7d","7d","7d",
                      "14d","14d","14d","14d","14d","14d",
                      "30d","30d","30d","30d","30d","30d"), levels = c("PSC","NC","7d","14d","30d"))
colData(pca.pseudocounts)$Diff.Day <-Diff.Days
Diff.Day <- colData(pca.pseudocounts)$Diff.Day



colData(pca.pseudocounts)$Cell.Type[sample.data$Sphere.Type %in% "Heterotypic" ] <- "Sphere"
colData(pca.pseudocounts)$Cell.Type[sample.data$Sphere.Type %in% "Homotypic" ] <- "Sphere"

colData(pca.pseudocounts)$Cell.Type <- factor(as.character(colData(pca.pseudocounts)$Cell.Type), levels=unique(colData(pca.pseudocounts)$Cell.Type))
Samples <- colData(pca.pseudocounts)$Cell.Type
Samples <- as.character(Samples)
Samples[Samples == "NF1"] <- "NF1(-/-)"
Samples <- factor(Samples,levels = unique(Samples))
colData(pca.pseudocounts)$Cell.Type <- Samples
# mycolors <- c(PSC = "orange",SC = "skyblue", NC="#F8766D", day7 = "#00BF7D", day14 = "#A3A500", day30 = "purple1")
mycolors <- unique(sample.data$colors.p1)
names(mycolors) <- unique(colData(pca.pseudocounts)$Diff.Day)

# plotPCA(pca.pseudocounts, intgroup=c("Diff.Day"))


pcaData <- plotPCA(pca.pseudocounts, intgroup=c("Diff.Day","Cell.Type", "Model"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))


png(filename=file.path(pca.dir,paste0(model, "_FiPS_MM_NOlabel_newshape_pca.p3.png")), width = 1500, height = 800)
ggplot(pcaData, aes(PC1, PC2, shape = Cell.Type, color = Diff.Day))+
  geom_point(fill = sample.data$colors.p1, size = 15, stroke =3)+ 
  # scale_fill_manual(name = "Stages", values = mycolors)+
  scale_color_manual(name = "Samples", values = mycolors)+
  # scale_shape_manual(name = "Genotype", values = c(23))+
  scale_shape_manual(name = "Genotype", values = c(23, 9))+
 
theme(legend.title = element_text(size=50), 
      legend.text = element_text( size=50), 
      axis.text=element_text(size=60),
      axis.title=element_text(size=60),
      plot.title = element_text(size=30),
      panel.background = element_rect(fill="white"), 
      panel.grid.major = element_line(colour = "grey50"), plot.margin =margin(30,10,10,10))+
  guides(colour = guide_legend(override.aes = list(shape = 15, size = 15)), 
         shape = guide_legend( override.aes = list(size = 8, color ="#7C8e8F", fill = "#7C8e8F" )))+
  xlab(paste0("PC1: ",percentVar[1],"% variance") ) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

dev.off()





       