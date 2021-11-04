##########################################################
#     Heatmap Representation Mazuelas et al. 2021        #
##########################################################
####### Packages Needed
library(DESeq2)
library(tximport)
library(org.Hs.eg.db)
library(yaml)
library(ggplot2)
library(ggbeeswarm)
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
source(file = "./rna_seq_Functions.R")

###### Parameters ######

#All samples parameters
params <- yaml.load_file("./Parameters/AllSamples_pipeline_parameters.yaml")

#Model of Analysis
model <- params$model
samples.group <-params$samples.group
stages <- params$stages

# loading the sample data information 
sample.data <- read.table(file = params$sample.data.file , header = T, sep = "\t", stringsAsFactors = FALSE, comment.char = "")
sample.data <- sample.data[sample.data$Genotype %in% c(params$genotype),]
sample.data <- sample.data[!sample.data$graph.Names %in% c( "", "FB_E"),]
sample.data$graph.Names <- gsub("FB","Fb",sample.data$graph.Names)

# Getting the file names of the analysis
file.names <- sample.data$File.Name
sample.names <- sample.data$graph.Names
# Salmon alignement and quantification parameters
file1.suffix <- params$file1.suffix
file2.suffix <- params$file2.suffix
fastqdir <- params$fastqdir
transcript.index <- params$transcript.index
output.suffix <- params$output.suffix
output.quants <- params$output.quants
# threads <- 8

# Tximport paramenters
orgdb <- org.Hs.eg.db
org.columns <- params$org.columns
org.keytype <- params$org.keytype

# DESeq2 parameters
filt.min.reads <- params$filt.min.reads
filt.min.samples <-params$filt.min.samples
pvalue <- params$pvalue

#Color heatmap
color.plate <-  bluered(75)

# Genes FiPS pipeline
gene.markers <- c() 
stages <- c("PSC", "NC", "day7", "day14", "day30")
for(i in seq_len(length(stages))){
  gr <- stages[i]
  mks<- read.table(file = file.path('./Results/2D/markers', paste0("up_", gr,".txt")), header = FALSE,stringsAsFactors = FALSE)[,1]
  # mks <- markers[[gr]]
  names(mks) <- rep(gr, length(mks))
  gene.markers <- c(gene.markers,mks)
  
}
# #Deleting PSC markers
# gene.markers <- gene.markers[names(gene.markers) != "PSC"]

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


########### Procesing the data for the heatmap representation ########
filtered.dds <- getFilteredDDS(tximport = txi.salmon, samples_group = samples.group, samples_df = sample.data,filter_min_reads = filt.min.reads,filter_min_samples = 2)
dds.rlog <- rlog(filtered.dds)

#### Heatmap Annotation
# All Samples
Diff.Days <- factor(c("PSC","PSC","PSC","PSC",
                      "NC","NC","NC","NC","NC","NC",
                      "7d","7d","7d","7d","7d","7d",
                      "14d","14d","14d","14d","14d","14d",
                      "30d","30d","30d","30d","30d","30d","WT_Heterotypic","WT_Heterotypic","WT_Heterotypic",
                      "NF1_Heterotypic","NF1_Heterotypic","NF1_Heterotypic",
                      "NF1_Homotypic","NF1_Homotypic","NF1_Homotypic","Fb","Fb","PNF SC","PNF SC","PNF SC"),
                    levels = c("PSC","NC","7d","14d","30d","WT_Heterotypic","NF1_Heterotypic","NF1_Homotypic","Fb","PNF SC"))
colData(dds.rlog)$Diff.Day <-Diff.Days
colnames(dds.rlog) <- gsub("FB","Fb",colnames(dds.rlog))

#Stage annotation
annot.col <- data.matrix( colnames(dds.rlog))
samples.group <- "Diff.Day"
colnames(annot.col)<- "Samples"
rownames(annot.col) <- annot.col[,1]
annot.col[,1] <- as.character(colData(dds.rlog)[,samples.group])
annotcol <- data.frame(annot.col)
sample_grououp <- unique(sample.data$colors.p1)
names(sample_grououp) <- unique(as.character(colData(dds.rlog)[,samples.group]))
sample_grououp <- sample_grououp[-1]

sample_grououp[5:9] <- "white"
names(sample_grououp)[5:9] <- ""
sample_grououp <- sample_grououp[-c(6:9)]

#Genotype Annotation
gno <- colData(dds.rlog)[colnames(dds.rlog),]

gno <- gno$Cell.Type

annotcol$Genotype <- gno
annotcol$Genotype[grepl("NF1",annotcol$Genotype)] <- "NF1(-/-)"
annotcol$Genotype[grepl("Fb",annotcol$Genotype)] <- "NF1(+/-)"

annotcol <- annotcol[,c(2,1)]
colnames(annotcol) <- c('Genotype   ',"Samples   ")

# Sample groups colors
mycolors.col <- c("grey44","grey70","grey89")
names(mycolors.col) <- c("WT", "NF1(+/-)","NF1(-/-)" )
mycolors <- list("Samples   " =sample_grououp, "Genotype   " = mycolors.col)

# No PSC
annotcol <- annotcol[-which(grepl("PSC",rownames(annotcol))),]
annotcol[,"Samples   "][25:nrow(annotcol)] <- ""

col.selected <- c(1,3,4,2,5:7,9,10,8,11:13,15,16,14,17:19,21,22,20,23:25,27,28,26,29:31,33,34,32,36,37,35,40:42,38:39) #Fips+Spheres+SC+FB

####### Data counts to plot #####
count.data <- assay(dds.rlog)[gene.markers,] # selection of FiPS markers
colnames(count.data)

# Data normalized to plot in the heatmap
data_subset_norm <- t(apply(count.data, 1, cal_z_score))
data_subset_norm <- data_subset_norm[,col.selected]
colnames(data_subset_norm)

#without iPSC
data_subset_norm <- data_subset_norm[,!grepl(x=colnames(data_subset_norm),pattern = "PSC")]

# Adjusting the zscore color associated with the different heatmaps
zscore <- data_subset_norm
heat.col <- colorRampPalette(c("blue","white","red"))(75)
myBreaks <- c(seq(min(zscore), 0, length.out=ceiling(75/2) + 1), 
              seq(max(zscore)/75, max(zscore), length.out=floor(75/2)))


# Repersentation of WT+SC (F1)
data_subset_norm <- data_subset_norm[,c(1:3,7:9,13:15,19:21,34:36)]
annotcol.f <- annotcol[c(1:3,7:9,13:15,19:21,36:38),]
annotcol <- annotcol.f
annotcol[grepl("^SC",rownames(annotcol)),2] <- "PNF SC"
annotcol[grepl("^SC",rownames(annotcol)),1] <- "NF1(-/-)"

sample.g <- sample_grououp
sample.g[length(sample.g)]<-"#F0CE77"
names(sample.g)[length(sample.g)] <- "PNF SC"
colnames(annotcol) <- c("Genotype   ","Samples   ")
mycolors <- list("Samples   " =sample.g, "Genotype   " = mycolors.col[-2])

# Repersentation of WT+NF1 2D (F2)
# data_subset_norm <- data_subset_norm[,c(1:24)]
# colnames(data_subset_norm)
# annot.nf <- annotcol[c(1:24),]
# sg <- sample_grououp
# sg <- sg[-length(sg)]
# gn <- mycolors.col
# gn <- gn[-2]
# mycolors <- list("Samples   " =sample_grououp, "Genotype   " = mycolors.col)
# annotcol <- annot.nf

# Representation of WT+spheres+sc+fb (F3)
# data_subset_norm <- data_subset_norm[,c(1:3,7:9,13:15,19:21,25:38)]
# annotcol[grepl("^SC",rownames(annotcol)),2] <- "PNF SC"
# annotcol[grepl("^Fb",rownames(annotcol)),2] <- "PNF Fb"
# sample.g <- sample_grououp
# names(sample.g)[length(sample.g)] <- "PNF SC"
# sample.g <- c(sample.g,"")
# sample.g[length(sample.g)]<-"white"
# names(sample.g)[length(sample.g)] <- "PNF Fb"
# sample.g <- c(sample.g,"")
# sample.g[length(sample.g)]<-"white"
# names(sample.g)[length(sample.g)] <- ""
# colnames(annotcol) <- c("Genotype   ","Samples   ")
# mycolors <- list("Samples   " =sample.g, "Genotype   " = mycolors.col)

colnames(data_subset_norm)

#heatmap of stage specific markers

# png(filename = file.path(file.path(heatmap.dir, paste0("FiPS_AllSpheres_supervised_","Markers_NoLegend_", model, ".colBar.2.png"))), width = 1000, height = 800)
png(filename = file.path(file.path(heatmap.dir, paste0("FiPS_supervised_","Markers_NoLegend_", model, ".png"))), width = 1000, height = 800)
# png(filename = file.path(file.path(heatmap.dir, paste0("FiPS_Markers_AllSamples_supervised_", model, "_3verticallLegend_colBar.2.png"))), width = 1000, height = 800)
pheatmap(data_subset_norm,
         # annotation_row = annotation.row,
         annotation_col = annotcol,
         annotation_colors = mycolors,
         cluster_cols = F,
         cluster_rows =F,
         fontsize = 20,
         # color = bluered(75),
         color =heat.col,
         fontsize_col = 25,
         show_rownames = FALSE,#gaps_col = c(12,15,18),
         # annotation_legend = TRUE,
         legend = TRUE,
         legend_breaks = c(round(min(zscore)),round(max(zscore))),
         breaks =myBreaks ,
         legend_labels = c("Down","Up"),
         annotation_names_row = FALSE,
         annotation_names_col = FALSE,
         lwd =3,
         treeheight_row = 40,
         treeheight_col = 40,
         angle_col = 90, 
         margins=c(50,50,50,50),
         # gaps_col = c(12)
         # gaps_col = c(12,15,18,21,24)
         # gaps_col= c(24,27,30,33)
         
         # main = paste0("Clustering all DEG in all samples", height," height ", nrow(annotation_row), " genes"))
         # main = c("Stage Specific Markers in FiPS")
)
dev.off()
