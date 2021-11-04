####################################################
###########      Enrichment plots      #############
####################################################

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
# FiPS NF1 2D parameters
params <- yaml.load_file("./Parameters/FiPS_NF1_2D_pipeline_parameters.yaml")

#Model of Analysis
model <- params$model
samples.group <-params$samples.group
stages <- params$stages

# Getting the file names of the analysis
file.names <- sample.data$File.Name
sample.names <- sample.data$graph.Names
goterms <- params$go.terms

# Salmon alignement and quantification parameters
file1.suffix <- params$file1.suffix
file2.suffix <- params$file2.suffix
fastqdir <- params$fastqdir
transcript.index <- params$transcript.index
output.suffix <- params$output.suffix
output.quants <- params$output.quants

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


#### loading the sample data information 
sample.data <- read.table(file = params$sample.data.file , header = T, sep = "\t", stringsAsFactors = FALSE, comment.char = "")
sample.data <- sample.data[sample.data$Genotype %in% c(params$genotype),]
sample.data <- sample.data[!sample.data$graph.Names %in% c( "", "FB_E"),]
sample.data$graph.Names <- gsub("FB","Fb",sample.data$graph.Names)


########## Directories
model.dir <- params$model.dir

markers.dir <- file.path("Results", model.dir, "markers")
if(!file.exists(markers.dir)) dir.create(markers.dir)

gsea.dir <- file.path("Results", model.dir, "GSEAGO")
if(!file.exists(gsea.dir)) dir.create(gsea.dir)
for(term in goterms){
  if(!file.exists(file.path(gsea.dir, term))) dir.create(file.path(gsea.dir, term))
}

annot.dir <- file.path("Results",model.dir,  "Annotation")
if(!file.exists(annot.dir)) dir.create(annot.dir)

gene.expression.dir <- file.path("Results",model.dir,  "GeneExpressionPlot")
if(!file.exists(gene.expression.dir)) dir.create(gene.expression.dir)


################## Loading Markers
# #Markers FiPS
# stages <- c("PSC", "NC", "day7", "day14", "day30")
# gene.markers <- c()
# for(i in seq_len(length(stages))){
#   gr <- stages[i]
#   markers<- read.table(file = file.path(markers.dir, paste0("up_", gr,".txt")), header = FALSE,stringsAsFactors = FALSE)[,1]
#   names(markers) <- rep(gr, length(markers))
#   gene.markers <- c(gene.markers,markers)
# 
# }
# gene.markers <- gene.markers[!names(gene.markers)%in% c("PSC")]

# #Markers NF1 -/-
# stages <- c("NC","day7", "day14", "day30")
# cells <- c("3MM", "5MM", "D12")
# gene.markers <- c()
# gene.markers.list <- list()
# 
# for(j in seq_len(length(cells))){
#   cell <- cells[j]
#   for(i in seq_len(length(stages))){
#     gr <- stages[i]
#     markers<- read.table(file = file.path(markers.dir, paste0("up_", cell ,"_",gr,".txt")), header = FALSE,stringsAsFactors = FALSE)[,1]
#     names(markers) <- rep(gr, length(markers))
#     gene.markers <- c(gene.markers,markers)
#     
#   }
#   gene.markers.list[[cell]] <- gene.markers
# }
# 
# gene.markers <- gene.markers[!names(gene.markers)%in% c("PSC")]

################## Loading Enrichment
#Enriched terms 2D NF1-/-
cells <- c("3MM", "5MM","D12")
enrichment.2DNF1 <- list()

for (i in seq_len(length(cells))){
  cell <- cells[i]
  enrichment.2DNF1[[cell]] <- read.table(file = file.path(gsea.dir,"BP", paste0("Enrichment_", cell, "_BPmarkers_up_day30.csv")), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
}

#Enrichment FiPS
stages <- c("PSC", "NC", "day7", "day14", "day30")
enrichment.FiPS <- list()
i=5
for (i in seq_len(length(stages))){
  stage <- stages[i]
  enrichment.FiPS[[stage]] <- read.table(file = file.path(gsea.dir, "BP", paste0("Enrichment_BP_markers_up_", stage, "_FiPS",".csv")), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
}

#####  Point representation of enriched terms  #####


### Day 30 NF1(-/-) cell lines enrichment
# sg <- "3MM"
# sg <- "5MM"
# sg < "D12"
# df.bp <- enrichment.2DNF1[[sg]]

### Day 30 Markers FiPS Enrichment
df.bp <- enrichment.FiPS$day30
df.bp <- data.frame(df.bp) [c(2,3,6,10,11,14,15),]#30Days Markers



df.bpg <- head(data.frame( Description= df.bp$Description,"log10pvalue"=-log10(df.bp$pvalue), ngenes = df.bp$Count ,Position = c(1:nrow(df.bp))),7)

wr.lap <- wrap.labels(df.bpg$Description, 60)
wr.sintetizado <- wr.lap

# #3MM
# wr.sintetizado[c(2,3,4,7)] <-"" 
# 
# p.size <- c(1:7)
# names(p.size)<- wr.sintetizado
# p.size[names(p.size)==""]<- 3
# p.size[names(p.size)!=""]<- 6
# p.color <-wr.sintetizado
# names(p.color)<- wr.sintetizado
# p.color[names(p.color)==""]<- "grey"
# p.color[names(p.color)!=""]<-"black"

# #D12
# wr.sintetizado[c(1,2,3,5)] <- ""
# p.size <- c(1:7)
# names(p.size)<- wr.sintetizado
# p.size[names(p.size)==""]<- 3
# p.size[names(p.size)!=""]<- 6
# p.color <-wr.sintetizado
# names(p.color)<- wr.sintetizado
# p.color[names(p.color)==""]<- "grey"
# p.color[names(p.color)!=""]<-"black"

# #5MM
# wr.sintetizado[c(1,2,3,7)] <- ""
# p.size <- c(1:7)
# names(p.size)<- wr.sintetizado
# p.size[names(p.size)==""]<- 3
# p.size[names(p.size)!=""]<- 6
# p.color <-wr.sintetizado
# names(p.color)<- wr.sintetizado
# p.color[names(p.color)==""]<- "grey"
# p.color[names(p.color)!=""]<-"black"

# FiPS 30dM
wr.sintetizado[c(3,2,5,7)] <- ""
p.size <- c(1:7)
names(p.size)<- wr.sintetizado
p.size[names(p.size)==""]<- 3
p.size[names(p.size)!=""]<- 6
p.color <-wr.sintetizado
names(p.color)<- wr.sintetizado
p.color[names(p.color)==""]<- "grey"
p.color[names(p.color)!=""]<-"black"


################## Enrichment plot representation ###########################
# svg(filename = file.path(gsea.dir,"BP",  paste0("BP_Enrichment_Selected_FiPS_30day_Markers_", model, ".svg")), width = 16, height = 9)
# svg(filename = file.path(gsea.dir,"BP",  paste0("BP_Enrichment_D12_day30_", model, ".svg")), width = 13, height = 9)
svg(filename = file.path(gsea.dir,"BP",  paste0("BP_Enrichment_",sg,"_day30_", model, ".svg")), width = 20, height = 9)

par(mar=c(10,10,5,2), mgp=c(5,2,0), oma=c( 0,0,0,0))
# plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', cex.lab = 2, ylab = '-log10(pvalue)', xlab = 'Biological Process Enrichment', xlim= c(0,nrow(df.bpg)+2), ylim = c(min(df.bpg$log10pvalue)-2,max(df.bpg$log10pvalue)+2))
# ticks.labels <- seq(as.integer(min(df.bpg$log10pvalue)-2),as.integer(max(df.bpg$log10pvalue)+mean(df.bpg$log10pvalue)+50), by =3)

ticks.labels <- round(seq(min(df.bpg$log10pvalue),max(df.bpg$log10pvalue)+1.5,0.5),1)

ticks.values <- (max(ticks.labels)-min(ticks.labels))/2
round.tick.labels <- seq(min(ticks.labels),max(ticks.labels),ticks.values)


plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', xlab="",ylab="", cex.axis =3, font =2,
     ylim= c(nrow(df.bpg)*10+2,0), xlim = c(min(ticks.labels),max(ticks.labels)))

# y side
axis(side = 2, lwd = 3, line = TRUE, lwd.ticks = 0, labels = FALSE, at = c(nrow(df.bpg)*10+20,-10))
mtext(side=2, at =(nrow(df.bpg)*10/2), line =3.5, 'Enrichment', cex=5)

# x side
axis(side = 1, lwd = 3, line = TRUE, lwd.ticks =5, labels = round.tick.labels, cex.axis = 4, at =round.tick.labels)
mtext(side=1, at =mean(c(min(df.bpg$log10pvalue)-0.5,max(round.tick.labels))) ,line = 6, '-log10(pvalue)', cex=5)

# points(y = df.bpg$Position*10, x=df.bpg$log10pvalue,
#        pch=16, cex =as.numeric(names(point.size)), col = "black")

points(y = df.bpg$Position*10, x=df.bpg$log10pvalue,
       pch=16, cex =p.size, col = p.color)


text(y =seq_len(length(df.bpg$Description))*10, pos=4, 
     x =df.bpg$log10pvalue+0.2 , labels = wr.sintetizado, cex =4)



# legend(y =c(-20,-20), 
#        x =c(min(df.bpg$log10pvalue),as.integer(max(df.bpg$log10pvalue)+mean(df.bpg$log10pvalue)+10)) ,
#        inset=c(0,-0.2), xpd=TRUE,  bty = "n", title = "Number of genes        ", 
#        legend = p.size, pch = 16, x.intersp = 1.2,  text.width = c(0,seq(1,as.integer(max(df.bpg$log10pvalue)+mean(df.bpg$log10pvalue)+10),length(p.size))),
#        pt.cex=as.numeric(names(p.size)), cex =3.5,xjust=0, yjust=0,horiz = TRUE)
# legend("top",xpd = TRUE,legend = p.size, pt.cex=as.numeric(names(p.size)),
#        pch = 16, x.intersp = 1.2,inset=c(0,-0.08),cex =3.5, xjust=0, yjust=0,  horiz = T, bty = "n",title = "N of genes")

# legend("topright",xpd = TRUE,legend = p.size, pt.cex=as.numeric(names(p.size)),
#        pch = 16, y.intersp = 1.5, x.intersp = 1.2,inset=c(0,0),cex =3.5, xjust=0, yjust=0,  bty = "n",title = "N of genes")

# legend(y =mean(df.bpg$Position),
#        x = as.integer(max(df.bpg$log10pvalue)+mean(df.bpg$log10pvalue)+11),
#        inset=c(0,-0.2), xpd=TRUE,  bty = "n", title = "Number of genes        ", 
#        legend = p.size, pch = 16, x.intersp = 1.2,  text.width = c(0,seq(1,as.integer(max(df.bpg$log10pvalue)+mean(df.bpg$log10pvalue)+10),length(p.size))),
#        pt.cex=as.numeric(names(p.size)), cex =3.5, xjust=0, yjust=0)


dev.off()





# 
# ##################### Code not used for the paper #########################
# ############ New Representation of Enriched BP points representation
# 
# # png(filename = file.path(gsea.dir,"BP",  paste0("BP_Enrichment_day30_MarkersgenesInSC_3_5_6_", model, ".png")), width = 1200, height = 800)
# png(filename = file.path(gsea.dir,"BP",  paste0("BP_Enrichment_",sg, "_day30_", model, "_NEW.png")), width = 1600, height = 900)
# # png(filename = file.path(gsea.dir,"BP",  paste0("BP_Enrichment_Selected_FiPS_30day_Markers_", model, ".png")), width = 1600, height = 900)
# 
# # svg(filename = file.path(gsea.dir,"BP",  paste0("BP_Enrichment_Selected_FiPS_30day_Markers_", model, ".svg")), width = 16, height = 9)
# svg(filename = file.path(gsea.dir,"BP",  paste0("BP_Enrichment_",sg,"_day30_", model, ".svg")), width = 13, height = 9)
# 
# par(mar=c(10,10,5,0), mgp=c(6,2,0), oma=c( 0,0,0,0))
# # plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', cex.lab = 2, ylab = '-log10(pvalue)', xlab = 'Biological Process Enrichment', xlim= c(0,nrow(df.bpg)+2), ylim = c(min(df.bpg$log10pvalue)-2,max(df.bpg$log10pvalue)+2))
# # ticks.labels <- seq(as.integer(min(df.bpg$log10pvalue)-2),as.integer(max(df.bpg$log10pvalue)+mean(df.bpg$log10pvalue)+50), by =3)
# ticks.labels <- round(seq(min(df.bpg$log10pvalue),max(df.bpg$log10pvalue)+1.5,0.5),2)
# # n.ticks.labels <-round(quantile(unique(round(df.bpg$log10pvalue,digits = 2))),2)
# 
# point.pos <- df.bpg$log10pvalue
# 
# 
# plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', xlab="",ylab="", cex.axis =3, font =2,
#      xlim= c(nrow(df.bpg)*10+2,0), ylim = c(min(ticks.labels),max(ticks.labels)))
# 
# # x side
# axis(side = 1, lwd = 5, line = TRUE, lwd.ticks = 0, labels = FALSE, at = c(nrow(df.bpg)*10+20,-10), font =2)
# mtext(side=1, at =(nrow(df.bpg)*10/2), line =3.5, 'Enrichment',   cex=4)
# 
# # y side
# axis(side = 2, lwd = 5, line = TRUE, lwd.ticks =5, labels = TRUE, cex.axis = 4, font =2, at =c(ticks.labels))
# mtext(side=2, at =mean(ticks.labels) ,line = 6, '-log10(pvalue)', cex=4)
# 
# points(x = df.bpg$Position*10, y=df.bpg$log10pvalue,
#        pch=16, cex =p.size, col = p.color)
# # points(x=df.bpg$Position, pch=21, cex =df.bpg$ngenes, col = c("grey", "grey", "black", "black", "black", "grey", " black" ,"black", "black"))
# 
# 
# text(x =df.bpg$Position*10 , pos=4, 
#      y =df.bpg$log10pvalue , labels = wr.sintetizado, cex =3)
# 
# 
# dev.off()
# 
# 
# 
# 
# # Size of the points
# p.size <- sort(unique(df.bpg$ngenes))
# names(p.size) <- seq(1,length(unique(df.bpg$ngenes)),1)*2
# point.size <- df.bpg$ngenes
# 
# for (i in seq_len(length(df.bpg$ngenes))){
#   names(point.size)[i] <- names(p.size)[p.size%in% point.size[i]]
# }



