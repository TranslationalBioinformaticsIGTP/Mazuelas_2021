##################################################################
#                        Expression Plots                        #
##################################################################
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
library(bezier)
#######################  Loading Functions #########################
source(file= "/imppc/labs/eslab/mmagallon/Projects/RNA-Seq-timecourse.2/Analysis/rna_seq_Functions.R")
source(file = "/imppc/labs/eslab/mmagallon/Projects/RNA-Seq-timecourse.2/Analysis/normalization_ExpPlots_Funtions.R")

bezier.interpolate <- function(m, cp.dist.type="equal", cp.dist=0.5, npoints=100) {
  
  if(cp.dist.type == "equal") {
    cp1 <- m
    cp1[,1] <- cp1[,1] + cp.dist
    cp1 <- cp1[-nrow(cp1),]
    
    cp2 <- m
    cp2[,1] <- cp2[,1] - cp.dist
    cp2 <- cp2[-1,]
  } else {
    x.dist <- c(m[,1],0)-c(0, m[,1])
    x.dist <- x.dist[2:(nrow(m))]
    
    cp1 <- m
    cp1[,1] <- cp1[,1] + c(x.dist*cp.dist,0)
    cp1 <- cp1[-nrow(cp1),]
    
    cp2 <- m
    cp2[,1] <- cp2[,1] - c(0,x.dist*cp.dist)
    cp2 <- cp2[-1,]
  }
  
  l <- list(my.points, cp1, cp2)
  bz.matrix <- do.call(rbind, l)[order(sequence(sapply(l, nrow))), ]
  
  t <- seq(0, nrow(m)-1, length=100)
  
  return(bezier(t=t, p=bz.matrix, deg = 3))
}

############################  Parameters  ##########################
#Model of Analysis
model <- "2D"
samples.group <-"Diff.Day"
condition <-samples.group

#2D and 3D
model <- c("2D","3D")
samples.group <- "Spheresvs2D"
# AllSamples
# model <- c("2D","3D", "Cells")
# samples.group <- "Model"

stages <- c("PSC", "NC", "day7", "day14", "day30")
goterms <- c("BP", "CC", "MF","KEGG")

# Directories
if(length(model) ==2) model <- "2Dvs3D"
if(length(model) ==3) model<- "AllSamples"
if(model == "Cells") model <- "FbvsSC"
pca.dir <-  file.path("Results", model, "PCAplot")
if(!file.exists(pca.dir)) dir.create(pca.dir)

distance.dir <- file.path("Results", model, "DistancePlots")
if(!file.exists(distance.dir)) dir.create(distance.dir)

heatmap.QC.dir <- file.path("Results", "heatmap.QC")
if(!file.exists(heatmap.QC.dir)) dir.create(heatmap.QC.dir)

heatmap.dir <- file.path("Results",model, "DEG_heatmap")
if(!file.exists(heatmap.dir)) dir.create(heatmap.dir)

DEG.dir <- file.path("Results", model, "DESeqResults")
if(!file.exists(DEG.dir)) dir.create(DEG.dir)

selected.genes.dir <- file.path("Results", model, "DEGgenes")
if(!file.exists(selected.genes.dir)) dir.create(selected.genes.dir)

markers.dir <- file.path("Results", model, "markers")
if(!file.exists(markers.dir)) dir.create(markers.dir)

gsea.dir <- file.path("Results", model, "GSEAGO")
if(!file.exists(gsea.dir)) dir.create(gsea.dir)
for(term in goterms){
  if(!file.exists(file.path(gsea.dir, term))) dir.create(file.path(gsea.dir, term))
}

revigo.dir <- file.path("Results", model, "REViGO")
if(!file.exists(revigo.dir)) dir.create(revigo.dir)

annot.dir <- file.path("Results",model,  "Annotation")
if(!file.exists(annot.dir)) dir.create(annot.dir)

gene.expression.dir <- file.path("Results",model,  "GeneExpressionPlot")
if(!file.exists(gene.expression.dir)) dir.create(gene.expression.dir)

# Salmon alignement and quantification parameters
# salmonDir <- "/soft/bio/salmon-1.1.0/bin/salmon"
file1.suffix <- "_1.fastq.gz"
file2.suffix <- "_2.fastq.gz"
fastqdir <- "/imppc/labs/eslab/mmagallon/Projects/RNA-Seq-timecourse.2/Data"
transcript.index <- "/imppc/labs/eslab/mmagallon/Projects/RNA-Seq-MPNSTcellLinesVsPNFCellLines/Results/Salmon/salmon_indexes_UCSC_hg38"
output.suffix <- "_quant"
output.quants <- "/imppc/labs/eslab/mmagallon/Projects/RNA-Seq-timecourse.2/Analysis/RNA-Seq_2D_3D_ipsDifferentiation/Results/Salmon"
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
###### loading the sample data information 
sample.data <- read.table(file = "/imppc/labs/eslab/mmagallon/Projects/RNA-Seq-timecourse.2/Sample.Info.AllSamples.3.csv", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "")
sample.data <- sample.data[sample.data$Model %in% c("2D"),]
sample.data <- sample.data[sample.data$Model %in% c("2D", "3D"),]
sample.data <- sample.data[sample.data$Genotype %in% "PP",]
# Delete Fibroblast from sample.data
# sample.data <-  sample.data[!sample.data$Cell.Type %in% c("Fb"),]
# No iPSC
# sample.data <- sample.data[!sample.data$Diff.Day %in% "PSC",]

file.names <- sample.data$File.Name
sample.names <- sample.data$graph.Names

# Loading the FiPS stage specific markers
gene.markers <- c() 
for(i in seq_len(length(stages))){
  gr <- stages[i]
  markers<- read.table(file = file.path('/imppc/labs/eslab/mmagallon/Projects/RNA-Seq-timecourse.2/Analysis/RNA-Seq_2D_3D_ipsDifferentiation/Results/2D/markers', paste0("up_", gr,".txt")), header = FALSE,stringsAsFactors = FALSE)[,1]
  names(markers) <- rep(gr, length(markers))
  gene.markers <- c(gene.markers,markers)
  
}

#FiPS stage markers in scRNAseq

df <- read.table(file.path("/imppc/labs/eslab/mmagallon/Projects/scRNA-seq_PNFs/results/PNF23/MarkersRoadmap/PNF23_RoadmapMarkersScoring.csv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

common.singcell.markers <- read.table("/imppc/labs/eslab/mmagallon/Projects/scRNA-seq_PNFs/results/commonRoadmapMarkersinPNFs.csv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
genes <- common.singcell.markers$Genes
# genes<-df$Gene[!df$Gene%in%genes]


genes <- c("SOX10","CDH19", "TWIST1","MEDAG","MPCL1")
# Lines as gene expression
# 
# gene.color <- c("#388D36", "#388D36","#EF171A","#EF171A","#650A61","#EF171A","#650A61","#EE9D3A","#650A61","#EF171A")
# genes <- c("NGFR","SOX10", "CDH19", "ITGA4","PLP1", "GAP43", "S100B", "EGR2", "PMP22","MPZ")
# names(gene.color) <- genes
# color.genes <- list( gene.color.red = gene.color [which(gene.color =="#EF171A")],
#                      gene.color.green = gene.color[which(gene.color =="#388D36")],
#                      gene.color.purple = gene.color[which(gene.color =="#650A61")],
#                      gene.color.orange = gene.color[which(gene.color =="#EE9D3A")])


# MUSCLE AND CHODROCYTE GENES
# genes <- c("MYOD1", "MYH3", "MYF5", "MYH2", "PAX7", "PAX3")
# 
# genes <- c( "ACAN", "FMOD", "BMP2", "CHAD", "HAPLN1", "COMP")
# 
# 
# gene.color <-   c("#808080","#FF8000","#0F99B2","#610051")
# names(gene.color) <- c( "WT", "NF1_A", "NF1_B", "NF1_C")


#Genes to plot by sample group 
Muscle <- c("MYOD1", "MYH3", "MYF5", "MYH2", "PAX7")
names(Muscle) <- rep("darkred",length(Muscle))
# markers.gr <- "Muscle"
# genes <- c( "ACAN", "FMOD", "BMP2", "CHAD", "HAPLN1", "COMP")
Chondrocyte <- c( "ACAN", "FMOD",  "CHAD", "HAPLN1", "COMP")
names(Chondrocyte) <- rep("darkolivegreen",length(Chondrocyte))
# markers.gr <- "Chondrocyte"

# genes <- c("SOX10", "CDH19", "PLP1", "GAP43", "S100B")
SchwannCell <- c("MPZ", "PLP1", "S100B")
# genes <- c("GAP43", "GJC3", "GFAP","POU3F1")

names(SchwannCell) <- rep("cadetblue4",length(SchwannCell))
# markers.gr <- "SchwannCell"

# genes <- c("SOX10", "CDH19", "PLP1", "GAP43", "S100B")
NewSchwannCell <- c("GAS7", "GFRA3", "PLPP1")
names(NewSchwannCell) <- rep("cadetblue4",length(NewSchwannCell))
# markers.gr <- "NewSchwannCell"
markers.gr <- list(Muscle = Muscle, Chondrocyte =Chondrocyte, SchwannCell=SchwannCell,NewSchwannCell=NewSchwannCell)


########Classic and new stage markers#####################
old.new.genes <- read.table(file = "Markers_single_RNA_oldNews.csv", sep = "\t",header = TRUE,  stringsAsFactors = FALSE)
gene.old <- old.new.genes[old.new.genes$type == "old",]
gene.new <- common.singcell.markers

old.nc <- gene.old[gene.old$stage == "NC","genes"]
old.nc <- c(old.nc,"ERBB3")
names(old.nc) <- rep(unique(sample.data$colors.p1[sample.data$Diff.Day == "NC"]),length(old.nc))

old.d7 <- gene.old[gene.old$stage == "day7","genes"]
old.d7 <- c(old.d7,"DHH")
names(old.d7) <- rep(unique(sample.data$colors.p1[sample.data$Diff.Day == "day7"]),length(old.d7))
old.d7 <- old.d7[c(1,3,5)]

old.d14 <- gene.old[gene.old$stage == "day14","genes"]
names(old.d14) <- rep(unique(sample.data$colors.p1[sample.data$Diff.Day == "day14"]),length(old.d14))

old.d30 <- gene.old[gene.old$stage == "day30","genes"]
old.d30 <- c(old.d30,"MPZ")
names(old.d30) <- rep(unique(sample.data$colors.p1[sample.data$Diff.Day == "day30"]),length(old.d30))
old.d30 <- old.d30[c(2,4,6)]
color.genes <- list(NC=old.nc, day7 = old.d7, day14 = old.d14, day30 = old.d30)

genes <- c(old.nc,old.d7,old.d14,old.d30)
gene.color <- c(names(old.nc),names(old.d7),names(old.d14),names(old.d30))
names(gene.color) <- genes


new.nc <- gene.new[gene.new$Stage == "NC","Genes"]
names(new.nc) <- rep(unique(sample.data$colors.p1[sample.data$Diff.Day == "NC"]),length(new.nc))
new.nc <- new.nc[c(1,5,6)]
new.nc[c(1,2,3)] <- c("SSUH2","SCML1","STMN1")

new.d7 <- gene.new[gene.new$Stage == "day7","Genes"]
names(new.d7) <- rep(unique(sample.data$colors.p1[sample.data$Diff.Day == "day7"]),length(new.d7))
new.d7 <- new.d7[new.d7 %in%c("EHBP1","RCAN1","FRMD6")]

new.d14 <- gene.new[gene.new$Stage == "day14","Genes"]
names(new.d14) <- rep(unique(sample.data$colors.p1[sample.data$Diff.Day == "day14"]),length(new.d14))

new.d30 <- gene.new[gene.new$Stage == "day30","Genes"]
names(new.d30) <- rep(unique(sample.data$colors.p1[sample.data$Diff.Day == "day30"]),length(new.d30))
new.d30 <- new.d30[new.d30 %in%c("GAS7","GFRA3","PLPP1")]

color.genes <- list(NC=new.nc, day7 = new.d7, day14 = new.d14, day30 = new.d30)
genes <- c(new.nc,new.d7,new.d14,new.d30)
gene.color <- c(names(new.nc),names(new.d7),names(new.d14),names(new.d30))
names(gene.color) <- genes

genes <-genes[genes %in% c("HDAC9","LMO4","PPP1R10","SCML1","SNAI1","SSUH2","STMN1","TOB1","GFRA1",
                           "BIRC2", "EHBP1","COTL1","FRMD6","GPAT3","RCAN1","SERPINE2","USP53",
                           "ABCA8","CCL2","COL6A1","DAG1","ENDOD1","FXYD1","FXYD3","GAS7","GFRA3","KLF9","LGI4","MIA","PLP1","PLPP1","S100B")]

table(genes %in% c("HDAC9","LMO4","PPP1R10","SCML1","SNAI1","SSUH2","STMN1","TOB1","GFRA1",
                   "BIRC2", "EHBP1","COTL1","FRMD6","GPAT3","RCAN1","SERPINE2","USP53",
                   "ABCA8","CCL2","COL6A1","DAG1","ENDOD1","FXYD1","FXYD3","GAS7","GFRA3","KLF9","LGI4","MIA","PLP1","PLPP1","S100B"))

genes <- c("HDAC9","LMO4","PPP1R10","SCML1","SNAI1","SSUH2","STMN1","TOB1","GFRA1",
          "BIRC2", "EHBP1","COTL1","FRMD6","GPAT3","RCAN1","SERPINE2","USP53",
          "ABCA8","CCL2","COL6A1","DAG1","ENDOD1","FXYD1","FXYD3","GAS7","GFRA3","KLF9","LGI4","MIA","PLP1","PLPP1","S100B")

genes <- c("TFAP2A","TWIST1","SNAI2","SNAI1")
# genes <- c("EHBP1","GAS7","FXYD1","GFRA3")
# new.d30
# gene.markers.nc <- c("SOX10","NGFR","SSUH2","PTX3","ARHGAP5","SOX6","SCD","CIT","GFRA1","RXRG")
# names(gene.markers.nc) <- rep(unique(sample.data$colors.p1[sample.data$Diff.Day == "NC"]),length(gene.markers.nc))
# 
# gene.markers.nc <- c("SOX10","NGFR","SSUH2","PTX3","ARHGAP5","SOX6","SCD","CIT","GFRA1","RXRG")
# names(gene.markers.nc) <- rep(unique(sample.data$colors.p1[sample.data$Diff.Day == "NC"]),length(gene.markers.nc))
# # gene.markers.d7 <- gene.markers[names(gene.markers) == "day7"]
# gene.markers.d7 <- c("CDH19","CADM1", "EHBP1","SPP1", "PMEPA1", "ITGB8")
# names(gene.markers.d7) <- rep(unique(sample.data$colors.p1[sample.data$Diff.Day == "day7"]),length(gene.markers.d7))
# 
# # gene.markers.d14 <- gene.markers[names(gene.markers) == "day14"]
# gene.markers.d14 <- c("PRSS23")
# names(gene.markers.d14) <- rep(unique(sample.data$colors.p1[sample.data$Diff.Day == "day14"]),length(gene.markers.d14))
# 
# # gene.markers.d30<- gene.markers[names(gene.markers) == "day30"]
# gene.markers.d30 <- c("S100B","PLP1","ERBB3","DAG1","FXYD1","GFRA3","LGI4","PDLIM4","MIA","FXYD3","GAS7")
# names(gene.markers.d30) <- rep(unique(sample.data$colors.p1[sample.data$Diff.Day == "day30"]),length(gene.markers.d30))

# genes <- c(gene.markers.nc,gene.markers.d7,gene.markers.d14,gene.markers.d30)
# 
# gene.color <- brewer.pal(4,"Blues")
# names(gene.color) <- c("NC","day7","day14","day30")





############################## Preparing the data ##############################

##################  Tximport  #####################
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
# genes <- c("TWIST1","CD34")

##########################  Preparing the data  ################################
#Normalization
filtered.dds <- getFilteredDDS(tximport = txi.salmon, samples_group = samples.group, samples_df = sample.data)

# dds.deg <- DESeq(object = dds, test = "Wald")
# dds.rlog <- rlog(filtered.dds)

# Sample Groups for sample normalization
group <- c("WT","NF1_A","NF1_B", "NF1_C")

# group <- c("WT")
# group <- c("WT","NF1_A","NF1_B", "NF1_C","WT_Het","NF1_A_Het","NF1_B_Het", "NF1_C_Het","NF1_A_Hom","NF1_B_Hom", "NF1_C_Hom")
genes <- unlist(markers.gr)
norm.data <- getNormalizedFiltExp(filtered.dds = filtered.dds,
                                  genes = genes, 
                                  sample.data = sample.data,
                                  norm.data.column = "Diff.Day",
                                  sample.groups = group)
m.list  <- getMeanNomrFiltExp(norm.data)
desv.list <- getSdNomrFiltExp(norm.data = norm.data)

norm.data.list <- split( norm.data , f = norm.data$sample.group)

# g <- 1
# color.genes$mix <- c(color.genes$gene.color.green["NGFR"],color.genes$gene.color.purple["S100B"])
######################## Gene expression Plot #############################
#  Line by group of genes
g="WT"
for( g in seq_len(length(group))){
  f <- norm.data.list[[g]]
  f$data.group <- factor(f$data.group, levels = c("PSC","NC","day7","day14", "day30"))
  m <- m.list[[g]]
  # i =2
  for(i in seq_len(length(color.genes))){
    # gn <- names(color.genes[[i]])
    sample.group <- names(color.genes)[i]
    gn <- color.genes[[sample.group]]
    gn <- gn[gn %in% colnames(f)]
    color <- names(gn)
    
    # color <- gene.color[gns]
    
    # png(filename = file.path(gene.expression.dir, paste0(names(color.genes)[i], "_Expression_Grupal_boxplot.png")), width = 1000, height = 800)
    # png(filename = file.path(gene.expression.dir, paste0(stage, "_new_Expression_Grupal_SmoothLine_nogenename.png")), width = 1000, height = 800)
    # pdf(file = file.path(gene.expression.dir, paste0(stage, "_old_Expression_Grupal_SmoothLine.pdf")), width = 1000, height = 800)
    svg(filename = file.path(gene.expression.dir, paste0(sample.group, "_NEW.2_Expression_Grupal_SmoothLine_nogenenameNoBold.test.svg")), width = 7, height = 5)
    par(oma=c( 0,0,0,0), mar =c(5,5,2,0)+0.1, mgp=c(5,1,0))
    
    plot(0,0, type= "n",
         xlim = c(-10, 35),
         # ylim = c(0,1+ max(colMaxs(as.matrix(desv[,-7]),na.rm = T))),
         ylim = c(0,1),# Raw Data
         xlab = "", 
         ylab = "",
         cex.lab= 4,
         cex.axis =4,
         xaxt ="n",
         yaxt = "n",
         bty="n", font =2)
         # main = paste0(  gr," ",stage, " markers"), cex.main = 4)
    
    axis(side = 1, lwd = 3, at= c(-15,-10,0,7,16,30,35),
         labels =FALSE,
         cex.axis=3)
    
    text(cex=3, x=c(-15,-10,0,7,16,30), labels =  c("","PSC","NC","7d"," 14d ","30d"),
         y=-0.2, xpd=TRUE, srt=0 )
    #y=-0.32
    axis(side = 2, at = c(0,0.5,1), cex.axis = 3,lwd = 3, labels =TRUE,las =1)
    
    for(j in seq_len(length(gn))){
      gene <- gn[j]
      
      my.points <- data.frame(c(-10,0,7,16,30),m[,gene])
      bz.points <- bezier.interpolate(my.points, cp.dist.type = "prop", cp.dist=0.5)
      lines(bz.points, type = "l",col = unique(color),lwd =5)
      
      text(x = 33.5, y = m[5,gene], labels = gene, cex=1.5)
      segments(x0 =c(-10,0,7,16,30) +0.05, x1 = c(-10,0,7,16,30) - 0.05, y0 = m[,gene], y1 = m[,gene], col = "black")
      # arrows( x0 = c(-10,0,7,16,30), m[,gene]-desv[,gene],
      #         x1 =c(-10,0,7,16,30), m[,gene]+desv[,gene], length = 0.05, angle=90, code= 3, lwd =2)


    }
   
    dev.off()
  }

}


### Representation of smooth lines of specific genes by group
i=2
# markers.gr <-"mesenquimal"
# names(genes)<-c("grey","grey")
m=1
for(m in seq_len(length(markers.gr))){
  m.gr <- markers.gr[[m]]
  m.nm <- names(markers.gr)[m]
  i=1
  for(i in seq_len(length(group))){
    gr <- group[i]
    f <- norm.data.list[[gr]]
    f$data.group <- factor(f$data.group, levels = c("PSC","NC","day7","day14", "day30"))
    # m <- m.list[[gr]]
    # mg <- m[grepl(gr,rownames(m)),]
    m <- m.list[[gr]]
    m <- m[,m.gr]
    desv<- desv.list[[gr]]
    desv <- desv[m.gr]
    # desv<- do.call(rbind,desv.list)
    # mgene <- max(m[,-ncol(m)],na.rm = TRUE)
    mgene <- max(m,na.rm = TRUE)
    
    # png(filename = file.path(gene.expression.dir, paste0(gr,"_",markers.gr,"_Expression_Grupal_2D_lines_Bezier_prop.png")), width = 1000, height = 800)
    svg(filename = file.path(gene.expression.dir, paste0(gr,"_",m.nm,"_Expression_Grupal_2D_lines_Bezier_prop_NoBold.svg")), width = 7, height = 5)
    
    par(oma=c( 0,0,0,0), mar =c(5,5,2,0)+0.1, mgp=c(5,1,0))
    
    plot(0,0, type= "n",
         xlim = c(-10, 35),
         # ylim = c(0,1+ max(colMaxs(as.matrix(desv[,-7]),na.rm = T))),
         ylim = c(0,1),# Raw Data
         xlab = "", ylab = "",
         cex.lab= 3.5,
         cex.axis = 3.5,
         xaxt ="n",
         yaxt = "n",
         bty="n",
         # main = paste0(  gr," ",m.nm, " markers"),
         cex.main =1, 
         font=2)
    
    axis(side = 1, lwd = 3, at= c(-15,-10,0,7,16,30,35),
         labels =FALSE,
         cex.axis=3)
    
    text(cex=3, x=c(-15,-10,0,7,16,30), labels =  c("","PSC","NC","7d"," 14d ","30d"),
         y=-0.2, xpd=TRUE, srt=0)
    
    axis(side = 2, at = c(0,0.5,1), cex.axis = 3,lwd = 3, labels =TRUE,las =1)
    
    for(j in seq_len(length(m.gr))){
      g <- m.gr[j]
      color <-names(m.gr)[j]
      my.points <- data.frame(c(-10, 0, 7, 16, 30),m[,g])
      bz.points <- bezier.interpolate(my.points, cp.dist.type = "prop", cp.dist=0.5)
      lines(bz.points, type = "l",col = unique(color),lwd =5)
      
      
      # lines(x =c(-10, 0, 7, 14, 30), y =m[,g], col = unique(color),lwd =3)
      segments(x0 =c(-10, 0, 7, 16, 30) +0.05, x1 = c(-10, 0, 7, 16, 30) - 0.05, y0 = m[,g], y1 =m[,g], col = "black")
      # arrows( x0 = c(-10, 0, 7, 14, 30),  m[,g]-desv[,g],
      # x1 = c(-10, 0, 7, 14, 30), m[,g]+desv[,g], length = 0.05, angle=90, code= 3, lwd =2)
      
      
    }
    dev.off()
  }
  
}


# Representation of mesenchymal markers

color.genes <- c("red","blue", "darkgreen", "orange")
names(color.genes) <- group

for( g in seq_len(length(group))){
  f <- norm.data.list[[g]]
  f$data.group <- factor(f$data.group, levels = c("PSC","NC","day7","day14", "day30"))
  m <- m.list[[g]]
  color <- color.genes[g]
  sample.group <- group[g]
  
  svg(filename = file.path(gene.expression.dir, paste0(sample.group, "_NEW.2_Expression_Grupal_SmoothLine_nogenenameNoBold.test.svg")), width = 7, height = 5)
  par(oma=c( 0,0,0,0), mar =c(5,5,2,0)+0.1, mgp=c(5,1,0))
  
  plot(0,0, type= "n",
       xlim = c(-10, 35),
       # ylim = c(0,1+ max(colMaxs(as.matrix(desv[,-7]),na.rm = T))),
       ylim = c(0,1),# Raw Data
       xlab = "", 
       ylab = "",
       cex.lab= 4,
       cex.axis =4,
       xaxt ="n",
       yaxt = "n",
       bty="n", font =2)
  # main = paste0(  gr," ",stage, " markers"), cex.main = 4)
  
  axis(side = 1, lwd = 3, at= c(-15,-10,0,7,16,30,35),
       labels =FALSE,
       cex.axis=3)
  
  text(cex=3, x=c(-15,-10,0,7,16,30), labels =  c("","PSC","NC","7d"," 14d ","30d"),
       y=-0.2, xpd=TRUE, srt=0 )
  #y=-0.32
  axis(side = 2, at = c(0,0.5,1), cex.axis = 3,lwd = 3, labels =TRUE,las =1)
  
  
  for(i in seq_len(length(genes))){
    gene <- genes[i]
    my.points <- data.frame(c(-10,0,7,16,30),m[,gene])
    bz.points <- bezier.interpolate(my.points, cp.dist.type = "prop", cp.dist=0.5)
    lines(bz.points, type = "l",col = unique(color),lwd =5)
    
    text(x = 33.5, y = m[5,gene], labels = gene, cex=1.5)
    segments(x0 =c(-10,0,7,16,30) +0.05, x1 = c(-10,0,7,16,30) - 0.05, y0 = m[,gene], y1 = m[,gene], col = "black")
    # arrows( x0 = c(-10,0,7,16,30), m[,gene]-desv[,gene],
    #         x1 =c(-10,0,7,16,30), m[,gene]+desv[,gene], length = 0.05, angle=90, code= 3, lwd =2)
    
    
  }
  dev.off()
}






colors <- c("red","blue", "darkgreen", "orange")
########## Boxplot and Line expression plot individual genes  ######
i =2
groups <-c("WT","NF1_A","NF1_B","NF1_C")
# groups <- "WT"
# markers.gr <- "Selected_markers"
m =1
k=1
i=1
for(m in seq_len(length(markers.gr))){
  m.gr <- markers.gr[[m]]
  m.nm <- names(markers.gr)[m]
  
  for(k in seq_len(length(groups))){
    grs <- groups[k]
    for(i in seq_len(length(m.gr))){
      # g <- names(gene.color)[i]
      g <- m.gr[i]
      if(!g%in% colnames(m.list[[grs]])) next
      
      color <- names(g)
      
      # png(filename = file.path(gene.expression.dir, paste0(g, "_",stg,"_",group,"_",markers.gr,"_FiPS_Expression_individual_SmoothLine.png")), width = 1000, height = 800)
      
      svg(filename = file.path(gene.expression.dir, paste0(g, "_",grs,"_",m.nm,"_FiPS_Expression_individual_SmoothLine_SUPP.svg")), width = 7, height = 5)
      
      par(oma=c( 0,0,0,0), mar =c(5,5,2,0)+0.1, mgp=c(5,1,0))
      
      plot(0,0, type= "n",
           xlim = c(-10, 35),
           # ylim = c(0,1+ max(colMaxs(as.matrix(desv[,-7]),na.rm = T))),
           ylim = c(0,1),# Raw Data
           xlab = "", ylab = "",
           cex.lab= 3.5,
           cex.axis = 3.5,
           xaxt ="n",
           yaxt = "n",
           bty="n",
           # main = paste0(  gr," ",m.nm, " markers"),
           cex.main =1, 
           font=2)
      
      axis(side = 1, lwd = 3, at= c(-15,-10,0,7,16,30,35),
           labels =FALSE,
           cex.axis=3)
      
      text(cex=3, x=c(-15,-10,0,7,16,30), labels =  c("","PSC","NC","7d"," 14d ","30d"),
           y=-0.2, xpd=TRUE, srt=0)
      
      axis(side = 2, at = c(0,0.5,1), cex.axis = 3,lwd = 3, labels =TRUE,las =1)
      
      my.points <- data.frame(c(-10, 0, 7, 16, 30),m.list[[grs]][,g])
      bz.points <- bezier.interpolate(my.points, cp.dist.type = "prop", cp.dist=0.5)
      lines(bz.points, type = "l",col = unique(color),lwd = 5)
      
      
      # lines(x =c(-10, 0, 7, 14, 30), y =m.list[[gr]][,g], col = unique(color),lwd =3)
      segments(x0 =c(-10, 0, 7, 16, 30) +0.05, x1 = c(-10, 0, 7, 16, 30) - 0.05, y0 = m.list[[grs]][,g], y1 = m.list[[grs]][,g], col = "black")
      arrows( x0 = c(-10, 0, 7, 16, 30), m.list[[grs]][,g]-desv.list[[grs]][,g],
              x1 = c(-10, 0, 7, 16, 30), m.list[[grs]][,g]+desv.list[[grs]][,g], length = 0.05, angle=90, code= 3, lwd =3)
      
      dev.off()
      
     
      
    }
  }
  
}



#Figure3 paper
groups <- c("WT_Het","NF1_A","NF1_B","NF1_C")

#colors
#Colors
colors <- group
names(colors)[grepl("WT", x = group)] <-"#808080"
# "#FF8000"
names(colors)[grepl("NF1_A", x = colors)] <- "#A35200"
names(colors)[grepl("NF1_A_Hom", x = colors)] <- karyoploteR::lighter("#A35200")
names(colors)[grepl("NF1_A_Het", x = colors)] <- karyoploteR::darker("#A35200")


names(colors)[grepl("NF1_B", x = colors)] <- "#0F99B2"
names(colors)[grepl("NF1_B_Hom", x = colors)] <- karyoploteR::lighter("#0F99B2")
names(colors)[grepl("NF1_B_Het", x = colors)] <- karyoploteR::darker("#0F99B2")


names(colors)[grepl("NF1_C", x = colors)] <- "#610051"
names(colors)[grepl("NF1_C_Hom", x = colors)] <- karyoploteR::lighter("#610051")
names(colors)[grepl("NF1_C_Het", x = colors)] <-  karyoploteR::darker("#610051")
j=3
i=1

for(j in seq_len(length(groups))){
  gr.wt <- "WT"
  gr.nf <- groups[j]
  fwt <- norm.data.list[[gr.wt]]
  fnf <-f <- norm.data.list[grepl(gr.nf,names(norm.data.list))]
  
  het<- names(fnf)[grepl("_Het",names(fnf))]
  hom <-names(fnf)[grepl("_Hom",names(fnf))] 
  
  f2ds <- f[[gr.nf]]
  f3dhoms<-c()
  if(length(hom)!=0){
    f3dhoms <- f[[hom]]
  } 
  f3dhets <- f[[het]]
  
  for(i in seq_len(length(genes))){
    # g <- names(gene.color)[i]
    #gene selection
    g <- genes[i]
    if(!g %in% colnames(fwt))next
    f2d <- f2ds[g]
    color2d <- names(colors)[colors == gr.wt]
    if(length(f3dhoms)>0){
      f3dhom <- f3dhoms[g]
      color3dhom <- names(colors)[colors == hom]
      
    }
    
    f3dhet <- f3dhets[g]
    color3dhet <- names(colors)[colors == het]
    
      # color <- gene.color[i]
    # png(filename = file.path(gene.expression.dir, paste0(g, "_", gr.nf,"_FiPS_Expression_individual_SmoothLine.2.png")), width = 1000, height = 800)
    svg(filename = file.path(gene.expression.dir, paste0(g, "_", gr.nf,"_FiPS_Expression_individual_SmoothLine.2.noBold.svg")), width = 7, height =5)
    
    par(oma=c( 0,0,0,0), mar =c(5,5,2,0)+0.1, mgp=c(5,1,0))
    
    plot(0,0, type= "n",
         xlim = c(-10, 35),
         # ylim = c(0,1+max(as.numeric(desv.list[["WT"]][,g]),na.rm = TRUE)),
         ylim = c(0,1),
         
         xlab = "", ylab = "",
         # cex.lab= 5,
         # cex.axis = 5,
         xaxt ="n",yaxt ="n",
         bty="n", main = paste0(g), cex.main = 2, font=2)
    
    axis(side = 1, lwd = 5, at= c(-15,-10,0,7,16,30),
         labels =FALSE,
         cex.axis=3,font = 2)
    
    text(cex=3, x=c(-15,-10,0,7,16,30), labels =  c("","PSC","NC","7d"," 14d ","30d"),
         y=-0.2, xpd=TRUE, srt=0)
    
    axis(side = 2, at = c(0,0.5,1), cex.axis = 3,lwd =5, labels =TRUE,las =1)
    
    my.points <- data.frame(c(-10, 0, 7, 16, 30), m.list[[gr.wt]][,g])
    bz.points <- bezier.interpolate(my.points, cp.dist.type = "prop", cp.dist=0.5)
    lines(bz.points, type = "l", col = color2d ,lwd =15)
    
    if(gr.nf != "WT_Het"){
      my.points <- data.frame(c(-10, 0, 7, 16, 30), m.list[[gr.nf]][,g])
      bz.points <- bezier.interpolate(my.points, cp.dist.type = "prop", cp.dist=0.5)
      lines(bz.points, type = "l", col = names(colors)[colors==gr.nf] ,lwd =15)
    }
    
    
    points(x = rep(16,nrow(f3dhet)), y = f3dhet[,g], pch = 21, bg =color3dhet, col = "black", cex=3.5  )
    
    if(length(f3dhoms)>0){
    points(x = rep(16,nrow(f3dhom)),y = f3dhom[,g], pch = 21,  bg = color3dhom, col = "black",cex=3.5  )
    }
    
   dev.off()
  }
  
}












png(filename = file.path(gene.expression.dir, "GreenColor_Markers_RTPCR.2.png"), width = 1000, height = 800)
# par(oma=c( 0,0,0,0), mar =c(5,5,3,1)+0.1, mgp=c(3,1,0))
# plotLimsAndLabels(gene_matrix = ms, gene = "SCP-iSC genes", cex.main = 20, cex.lab= 10, cex.axis = 8)
plot(0,0, type= "n", 
     xlim = c(min(as.numeric(colnames(ms)))-5, max(as.numeric(colnames(ms)))+5),
     ylim = c(0,max((ms),na.rm = TRUE)), 
     xlab = "Stages", ylab = "Expression",
     cex.lab= 3,
     cex.axis = 3,
     xaxt ="n",
     bty="n")

mtext(text = "Stages",
      side = 1,line = 30,cex = 20)
mtext(text = "Expression", cex = 20,
      side = 2, 
      line = 25)
axis(side = 1, lwd =3, at= c(-10,-5,0,7,14,30,35), labels = c("","PSC","NC","day7","day14","day30",""), cex.axis= 3)
axis(side = 2, lwd = 3, labels = FALSE)

# k=3
for (k in 1:length(gene.color.green)){
  
  g <- names(gene.color.green)[k]
  
  A <- matrix.list[["NGFR"]]
  
  p <- colMedians(A,na.rm = T)
  colMedians()
  plotPointsAndMeans(gene_matrix = A, points_col = gene.color,  line_col = gene.color.green[g], lty = 1, lwd =5, plot_points = FALSE, plot_legend = F)
  text(x=35, y = mean(A[,"30"],na.rm = TRUE), labels =g, cex = 2, pos = 2, col = "black")
  arrows( x0 = as.numeric(names(p)), p-desv[[g]],
          x1 = as.numeric(names(p)), p+desv[[g]],angle=90, code= 3, lwd =2)
  # so turn off clipping:
  # par(xpd=TRUE)
  #legend("topright",xpd = T, legend =names(gene.color.orange),fill = gene.color.orange, cex = 6, ncol = 3)
  # legend_size = 6, lwd = 15, pch = 19, points_size = 15,legend_locus = "topright"
}
dev.off()






# # Barplot gene expression
# png(filename = file.path(file.path(markers.dir, paste0("Barplotsclineage.png"))),
#     width = 3000, height = 2500)
# par(mar = c(20, 20, 20, 20)+ 0.1)
fct <- norm.data
fff <- fct[grepl("NF1", rownames(fct)),]
fff <- fff[c(1,4,7,10,13,16,19,2,5,8,11,14,17,20,3,6,9,12,15,18,21),]
rownames(fct)
#Representation correcta

fct <- norm.data
ff <- fct[c(1,5:7,11:13,17:19,23:25,29:31,2,8,14,20,26,32,35,3,9,15,21,27,33,36,4,10,16,22,28,34,37),]
fff <- ff[c(8:10,14:16,27,30,29,34,37,36,20,23,22),]#14d2D-3Dhom-3DHet
sample.group <- fff$sample.group
data.group <- fff$data.group


fff <- data.matrix(fff)
rownames(fff) <- c("WT_14d","WT_14d","WT_14d","WT_Het_3D","WT_Het_3D","WT_Het_3D",
                  "NF1_B_14d","NF1_B_Hom_3D", "NF1_B_Het_3D",
                   "NF1_C_14d", "NF1_C_Hom_3D", "NF1_C_Het_3D",
                   "NF1_A_14d", "NF1_A_Hom_3D", "NF1_A_Het_3D")

#Colors
barplot.colors <- rownames(fff)
names(barplot.colors)[grepl("WT", x = barplot.colors)] <-"grey"
names(barplot.colors)[grepl("WT_Het", x = barplot.colors)] <-"grey17"

# "#FF8000"
names(barplot.colors)[grepl("NF1_A", x = barplot.colors)] <- "#A35200"
names(barplot.colors)[grepl("NF1_A_Hom", x = barplot.colors)] <- karyoploteR::lighter("#A35200")
names(barplot.colors)[grepl("NF1_A_Het", x = barplot.colors)] <- "#3F2000"


names(barplot.colors)[grepl("NF1_B", x = barplot.colors)] <- "#0F99B2"
names(barplot.colors)[grepl("NF1_B_Hom", x = barplot.colors)] <- karyoploteR::lighter("#0F99B2")
names(barplot.colors)[grepl("NF1_B_Het", x = barplot.colors)] <- "#053740"


names(barplot.colors)[grepl("NF1_C", x = barplot.colors)] <- "#610051"
names(barplot.colors)[grepl("NF1_C_Hom", x = barplot.colors)] <- karyoploteR::lighter("#610051")
names(barplot.colors)[grepl("NF1_C_Het", x = barplot.colors)] <-  "#25001F"


ff <- fff[,-c(ncol(fff),ncol(fff)-1)]
i =1
rownames(ff)[grepl("WT",rownames(ff))] <- "WT"
rownames(ff)[grepl("NF1_A",rownames(ff))] <- "NF1_A"
rownames(ff)[grepl("NF1_B",rownames(ff))] <- "NF1_B"
rownames(ff)[grepl("NF1_C",rownames(ff))] <- "NF1_C"

# Data to barplot
# 
# dm <- matrix(nrow = 4,ncol = 9)
# colnames(dm) <- c(rep("2D",3),rep("Hom",3),rep("Het",3))
# rownames(dm) <- c("WT", "NF1_B", "NF1_C", "NF1_A")
# fg <- ff[,g]
# i =2
# for(i in seq_len(nrow(dm))){
#   rn <- rownames(dm)[i]
#   f.data <- fg[grepl(rn,names(fg))]
#   dm[i,c(1:3)] <- f.data[grepl("14d", names(f.data))]
#   
#   if(!grepl("Hom", names(f.data))){
#     dm[i,c(4:6)] <- f.data[grepl("Het", names(f.data))]
#   }else{
#     dm[i,c(4:6)] <- f.data[grepl("Hom", names(f.data))]
#     
#   }
#   
#   if(rn =="WT"){
#     dm[i,c(7:9)] <- 0
#     
#   }else{
#     dm[i,c(7:9)] <- f.data[grepl("Het", names(f.data))]
#     
#     dm[i,c(1,2,4,5,7,8)] <-0
#   }
#   
# }

for(i in seq_len(ncol(ff))){
  g <- colnames(ff)[i]
  
  
  
  # png(filename = file.path(gene.expression.dir, paste0(g, "_","NF1_WT_14dAllSamples_Expression_individual_Barplot.png")), width = 1200, height = 800)
  svg(filename = file.path(gene.expression.dir, paste0(g, "_","NF1_WT_14dAllSamples_Expression_individual_Barplot.svg")), width = 7, height = 5)
  par(mar=c(5,2,2,0))
  barCenters <-barplot(ff[,g], beside=TRUE,
                       cex.names=2,
                       las=2,
                       ylim=c(0,1), 
                       yaxt = "n",
                       cex.axis = 3,
                       cex.lab = 3,
                       cex.main =1,
                       font=2,
                       col = names(barplot.colors),
                       main = g,
                       lwd = 3,
                       names.arg ="",
                       space = c(0,0,0,0,0,0,1.5,
                                 0,0,1.5,0,0,1.5,0,0))
  axis(side = 2, at = c(0,0.5,1), cex.axis = 3,lwd = 8, labels =TRUE,font = 2,las =1)
  
  text(cex=3, x=c(3,9,13.5,18), labels =  c("WT", "NF1_B", "NF1_C", "NF1_A"),
       y=-0.2, xpd=TRUE, srt=45, font =2)

  dev.off()
}
group




#All Sampres ordered
ff <- fct
ff <-ff[c(1,3,4,2,5:7,9,10,8,11:13,15,16,14,17:19,21,22,20,23:25,27,28,26,29:31,33,34,32,36,37,35),]#FiPS+NF1_2D
# ff$cell.tyype <- sample.data$Cell.Type

#Mean by stage
ff$cell.tyype[c(29:34)] <- c("Het","Het","Het","Het","Het","Het")

ff$cell.tyype[c(35:37)] <- c("Hom","Hom","Hom")
ff$cell.tyype <- paste0(ff$cell.tyype, "_",ff$sample.group)
ff <- ff[,-c(85:86)]
ff.l <- split(ff,ff$cell.tyype)
fd <- lapply(ff.l, function(x) c(colSds(x[,-ncol(x)])))
fd <- as.data.frame(do.call("rbind", fd))
ff.l <- lapply(ff.l, function(x) c(colMeans(x[,-ncol(x)])))
ff.l <- as.data.frame(do.call("rbind", ff.l))
ss <- unlist(m.list)
class(ss)

#14d samples
ff <- rbind(fct[grepl("14d",rownames(fct)),],fct[grepl("3D",rownames(fct)),])
bar.colors.2 <- c(bar.colors.2[grepl("14d",names(bar.colors.2))],bar.colors.2[grepl("3D",names(bar.colors.2))])
ff <- ff.l
#Get colors
bar.colors.2<-rownames(ff)
names(bar.colors.2)<- rownames(ff)
rownames(fct)

bar.colors.2 <- sample.data$colors.p1
names(bar.colors.2) <- sample.data$Diff.Day
bar.colors.2<-c(unique(bar.colors.2)[c(1:5)],unique(bar.colors.2))
bar.colors.2
ff.l <-ff.l[order(rownames(ff.l),decreasing =T),]
ff.l <- ff.l[c(1:3,5,4,6:8,10,9,12,13,11),]
ff <- ff.l
for(i in seq_len(length(colors))){
  colors
  b.col <- names(colors)[i]
  g <- colors[i]
  
  bar.colors.2[grepl(g,names(bar.colors.2))] <- b.col
  
}
rownames(ff)

#14d
ff <- ff[c(4,9,11:13),]
bar.colors.2 <- bar.colors.2[c(4,9,10:13)]
rownameS(ff)
ff <- fff[,-c(ncol(fff),ncol(fff)-1)]
for(i in seq_len(ncol(ff))){
  g <- colnames(ff)[i]
  # png(filename = file.path(gene.expression.dir, paste0(g, "_","NF1_Expression_individual_Barplot.png")), width = 1200, height = 800)
  # par(mar=c(15,10,10,10))
  # 
  # barCenters <-barplot(fff[,g], beside=TRUE,
  #                      cex.names=3,
  #                      las=2,
  #                      ylim=c(0,1), 
  #                      cex.axis = 3,
  #                      cex.lab = 3,cex.main =5,font=2,
  #                      col = bar.colors,main = g,
  #                      names.arg =rownames(fff))
  # dev.off()
  # 
  
  # png(filename = file.path(gene.expression.dir, paste0(g, "_","NF1_WT_14dAllSamples_Expression_individual_Barplot.png")), width = 1200, height = 800)
  # par(mar=c(20,10,10,10))
    barCenters <-barplot(ff[,g], beside=TRUE,
                       cex.names=2,
                       las=2,
                       ylim=c(0,1), 
                       cex.axis = 3,
                       cex.lab = 3,cex.main =5,font=2,
                       col = names(barplot.colors),main = g,
                       names.arg =rownames(ff))
  # dev.off()
}

# arrows( barCenters, boxplot.data.g-boxplot.data.sd,
#         barCenters, boxplot.data.g+boxplot.data.sd,angle=90, code= 3, lwd =2)
# mtext(side=2, line=8, "log2(Expression count)", font=2, cex=3)
# legend("topleft", fill = color.data, legend = names(color.data), horiz = T, cex = 5)
# dev.off()
# 

bar.colors <- c(bar.colors, rownames(fct)



colData(filtered.dds)$Diff.Day

matrix.list <- list()
ms <- data.frame()
for(i in 1:length(genes)){
  gene <- genes[i]
  
  matrix.list[[gene]]<- generationGeneMatrix(dds = filtered.dds, 
                                             gene = gene, 
                                             sname_variable = "graph.Names",
                                             contrast_group = "Diff.Day",
                                             normalize =T)
  
}  
maxs <-unlist(lapply(matrix.list,function(x) max(x,na.rm = T)))


for(i in names(matrix.list)){
  matrix.list[[i]] <- matrix.list[[i]]/maxs[i]
  ms <- rbind(ms, matrix.list[[i]]) 
}
desv <- lapply(matrix.list, function(x){
  apply(x,2, function(y) sd(y,na.rm = TRUE))
})

colnames(ms) <- c("-5", "0","7", "14", "30")
matrix.list <- lapply(matrix.list,function(x){
  colnames(x)<- c("-5", "0","7", "14", "30")
  x }) 






