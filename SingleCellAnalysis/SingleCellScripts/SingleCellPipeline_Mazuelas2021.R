#########################################################
##    Single cell Pipeline Mazuelas et al. 2021        ##
#########################################################

######################## Packages Needed #########################
library(DropletUtils)
library(scater)
library(scran)
library(cowplot)
library(Seurat)
library(gplots)
library(EnsDb.Hsapiens.v86)
library(pheatmap)
library(SingleR)
library(TSCAN)
library(karyoploteR)
edb <- EnsDb.Hsapiens.v86
############################  Parameters  ##########################
genes.of.interest <- c("NF1", "EDNRB", "GAS7", "KLK8", "SOX8", "EN1", "EN2", "GFAP", "POU3F1", "GLB1" )
group <- "ClassicRoadmap_SC_DEG"

##### Markers Single cell ####### 
# mkrs.dir <- file.path("./Results/2D/ExpressionCluster",group,"Cluster_Exp_profile","Single_cell_Expression")
# if(!file.exists(mkrs.dir))dir.create(mkrs.dir)
# mkrs <- read.table(file.path("./Results/2D/ExpressionCluster", paste0(group, "_singleCell.csv")),sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# genes  <- unique(mkrs[,1])
# genes <- gene.markers[gene.markers%in% mkrs$GENES]
# mkrs[mkrs$GENES %in% g,]
# gene.markers <- c()
# stages <-  c("PSC", "NC","day7", "day14", "day30")
# for(i in seq_len(length(stages)-1)){
#   gr <- stages[i]
#   for(j in (i+1):(length(stages))){
#     gr2 <- stages[j]
#     markers<- read.table(file = file.path("/imppc/labs/eslab/mmagallon/Projects/RNA-Seq-timecourse.2/Analysis/RNA-Seq_2D_3D_ipsDifferentiation/Results/2D/DEGgenes", paste0("Genes_", gr,"vs",gr2,".csv")),header = FALSE,stringsAsFactors = FALSE)[,1]
#     names(markers) <- rep(gr, length(markers))
#     gene.markers <- c(gene.markers,markers)
#     gene.markers <- unique(gene.markers)
#   }
#   
# }
stages <- c("PSC", "NC", "day7", "day14", "day30")
# stages <- c( "NC", "day14", "day30")

# Loading the FiPS stage specific markers
gene.markers <- c() 
for(i in seq_len(length(stages))){
  gr <- stages[i]
  markers<- read.table(file = file.path('./Results/2D/markers', paste0("up_", gr,".txt")), header = FALSE,stringsAsFactors = FALSE)[,1]
  names(markers) <- rep(gr, length(markers))
  gene.markers <- c(gene.markers,markers)
  
}
color.plate <-  bluered(75)
sample.names <-  c("PNF19", "PNF20", "PNF23")
# sample.name <- "PNF19"
# sample.name <- "PNF20"
sample.name <- "PNF23"

#############################  FUNCTIONS  #########################

source(file= "./rna_seq_Functions.R")
source(file= "./SingleCellScripts/functionsSingleCellRnaSeq.R")

########################  Directories ########################
results.dir <-"./resultsSingleCell"
results.dir <- file.path(results.dir,sample.name)
if(!file.exists(results.dir)) dir.create(results.dir)
# out.dir <- "./results"
cell.cycle.dir <- file.path(results.dir,"CellCycleStatus")
if(!file.exists(cell.cycle.dir)) dir.create(cell.cycle.dir)
markers.roadmap.dir <- file.path(results.dir,"MarkersRoadmap")
if(!file.exists(markers.roadmap.dir)) dir.create(markers.roadmap.dir)

markers.scRNAseq <- file.path(results.dir,"MarkerSingleCell")
if(!file.exists(markers.scRNAseq)) dir.create(markers.scRNAseq)

trajectory.dir <- file.path(results.dir,"Trajectory")
if(!file.exists(trajectory.dir)) dir.create(trajectory.dir)

##############################  LOAD DATA ###########################################################

#Loading the procesed data
dd <- list()
for(sn in sample.names) {
  load(paste0("/imppc/labs/eslab/bgel/Projects/scRNAseq/TesisHelena/", sn, ".sce.hvg.Rdata"))
  dd[[sn]] <- sce.hvg
}
sce.hvg <- dd[[sample.name]]
colnames(sce.hvg) <- paste0("cell_",1:ncol(sce.hvg))

# #Load the FILTERED (only cells, no all droplets) UMI count matrix
# #sce <- read10xCounts("/imppc/labs/eslab/bgel/Projects/scRNAseq/AnalyseSamples/PNF20/PNF20/outs/filtered_feature_bc_matrix/")
# #sce <- read10xCounts("/imppc/labs/eslab/bgel/Projects/scRNAseq/AnalyseSamples/PNF23/PNF23/outs/filtered_feature_bc_matrix/")
# 
# data.file <- paste0("/imppc/labs/eslab/bgel/Projects/scRNAseq/AnalyseSamples/", sample.name, "/", sample.name, "/outs/filtered_feature_bc_matrix/")
# sce <- read10xCounts(data.file)
# # #Filtering those genes with 0 counts in all cells
# # genes.to.remove <- which(rowSums(assays(sce)$counts) == 0)
# # sce <- sce[-genes.to.remove,]
# 
# ###################################  PROCESS DATA ####################################################
# 
# #Quality Control
# gene.chr <- select(EnsDb.Hsapiens.v86,  keys=rownames(sce), keytype="GENEID", column="SEQNAME")
# head(gene.chr)
# is.mito <- which(gene.chr[,"SEQNAME"]=="MT")
# head(is.mito)
# qcstats <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))
# filtered <- quickPerCellQC(qcstats, percent_subsets="subsets_Mito_percent")
# colSums(as.matrix(filtered))
# 
# 
# sce <- sce[, !filtered$discard]
# # colnames(sce) <- paste0("cell_",1:ncol(sce))
# 
# 
# #rename genes to Symbols (so it's easier to select specific markers)
# symbols <- mapIds(edb, keys = rownames(sce), keytype = "GENEID", column = "SYMBOL")
# #Some genes did not get a symbol correctly, use the ENG for them
# no.symbol <- which(is.na(symbols))
# symbols[no.symbol] <- rownames(sce)[no.symbol]
# rownames(sce) <- symbols
# 
# # is.mito <- grepl(rownames(sce), pattern = "^MT-")
# # table(is.mito)
# # rownames(sce)[is.mito]

# #Normalzation
# size.factors <- librarySizeFactors(sce)
# hist(log10(size.factors))
# 
# sce <- logNormCounts(sce)#since no size.factors parameter is provided, it will internally run  librarySizeFactors and apply those
# 
# 
# #FEature Selection: select the 10% most variable genes
# gene.var <- modelGeneVar(sce)
# high.variance.genes <-  getTopHVGs(gene.var, prop=0.1)
# sce.hvg <- sce[unique(c(high.variance.genes, genes.of.interest))]
# # sce.hvg <- sce[rownames(sce) %in%unique(c(high.variance.genes, gene.markers))]
# # sce.hvg <- dd$PNF19
# str(high.variance.genes)
# 
# 
# # Dimensionality reduction.
# 
# set.seed(123456)
# sce.hvg <- runPCA(sce.hvg, ncomponents=25)
# sce.hvg <- runUMAP(sce.hvg, dimred = 'PCA') #, external_neighbors=TRUE)
# if(sample.name=="PNF23") {
#   set.seed(123456)
#   perplexity <- 80
# } else if(sample.name=="PNF20") {
#   set.seed(42) #123458)
#   perplexity <- 80
# } else if(sample.name=="PNF19") {
#   set.seed(42) #123458)
#   perplexity <- 110
# }
# 
# sce.hvg <- runTSNE(sce.hvg, perplexity=perplexity)
# #Clustering
# #Cluster with a graph based approach using k-nearest neighbours
# if(sample.name=="PNF23")  {
#   k <- 40
# } else if(sample.name=="PNF20") {
#   k=70
# } else if(sample.name=="PNF19") {
#   k=65
# }
# g <- buildSNNGraph(sce.hvg, k=k, use.dimred = 'PCA', type="jaccard") #applying jaccard PCA type
# # clust <- igraph::cluster_walktrap(g)$membership
# # clust.louvain <- igraph::cluster_louvain(g)$membership
# clust.infomap <- igraph::cluster_infomap(g)$membership
# 
# # colData(sce.hvg)[,"default"] <- factor(clust)
# # colData(sce.hvg)[,"louvain"] <- factor(clust.louvain)
# # colData(sce.hvg)[,"infomap"] <- factor(clust.infomap)
# 
# #### Assigning labels to sce
# colData(sce.hvg)[,"labels"] <- factor(clust.infomap)
# colLabels(sce.hvg) <-  factor(clust.infomap)

# png(filename = file.path(results.dir, paste0("clustering_TSNE_",sample.name,".4.png")))
plotReducedDim(sce.hvg, dimred="TSNE", colour_by = "labels",text_by="labels")
# dev.off()


#####################  Cell type association  ########################
#ref <- BlueprintEncodeData()
ref <- HumanPrimaryCellAtlasData()
pred <- SingleR(test=sce.hvg, ref=ref, labels=ref$label.fine) #label.main
# table(pred$labels)

# plotScoreHeatmap(pred)

#How do cell types relate to clusters?
tab <- table(Assigned=pred$pruned.labels, Cluster=colData(sce.hvg)[,"labels"])

# Adding a pseudo-count of 10 to avoid strong color jumps with just 1 cell.
# p4 <- pheatmap(log2(tab+10),fontsize = 10, color=colorRampPalette(c("white", "blue"))(101),cellwidth = 5,cellheight = 2)
# save_plot(file.path(results.dir, paste0(sample.name, ".cluster_identity_map.4.png")), plot = p4, base_height = 10, base_asp = 1)
png(file.path(results.dir, paste0(sample.name, ".cluster_identity_map.supp.png")),width =600, height = 950)
pheatmap(log2(tab+10),fontsize = 10,fontsize_col = 15, color=colorRampPalette(c("white", "blue"))(101))
dev.off()    


####### Cell type TSNE ######
#colors
sc.color <- "#8DD35F"
fb.color <- "#FF6600"
end.color <- "#D8AC24"
bc.color <- "#8CE4FC"

sample.name = "PNF20"
# sample.name = "PNF19"
# sample.name = "PNF23"
sce.hvg <- dd[[sample.name]]
colnames(sce.hvg) <- paste0("cell_",1:ncol(sce.hvg))

#PNF20
if(sample.name =="PNF20"){
  sc.cells <- which(sce.hvg$labels %in%c(2,16,14,17,5,3,10,7))#PNF20
  length(sc.cells)
  names(sce.hvg$labels)[sc.cells] <-sc.color
  fb.cells <- which(sce.hvg$labels %in%c(1,9,13,6,8))
  names(sce.hvg$labels)[fb.cells] <-fb.color
  length(fb.cells)
  end.cells <-which(sce.hvg$labels %in%c(12))
  length(end.cells)
  names(sce.hvg$labels)[end.cells] <-end.color
  blood.cells <- which(sce.hvg$labels %in%c(11,15,4))
  length(blood.cells)
  names(sce.hvg$labels)[blood.cells] <-bc.color
  
}

if(sample.name == "PNF19"){
  sc.cells <- which(sce.hvg$labels %in%c(1,2,5,7,13,11,14,12))#PNF20
  length(sc.cells)
  names(sce.hvg$labels)[sc.cells] <-sc.color
  fb.cells <- which(sce.hvg$labels %in%c(9,15,4,3))
  names(sce.hvg$labels)[fb.cells] <-fb.color
  length(fb.cells)
  end.cells <-which(sce.hvg$labels %in%c(16,6,8))
  length(end.cells)
  names(sce.hvg$labels)[end.cells] <-end.color
  blood.cells <- which(sce.hvg$labels %in%c(17,10,18))
  length(blood.cells)
  names(sce.hvg$labels)[blood.cells] <-bc.color
  
}

#PNF23
if(sample.name =="PNF23"){
  sc.cells <- which(sce.hvg$labels %in%c(5,7,2,6,3,4,10,17))
  length(sc.cells)
  names(sce.hvg$labels)[sc.cells] <-sc.color
  fb.cells <- which(sce.hvg$labels %in%c(9,1,13,12,15))
  names(sce.hvg$labels)[fb.cells] <-fb.color
  length(fb.cells)
  end.cells <-which(sce.hvg$labels %in%c(8,14))
  length(end.cells)
  names(sce.hvg$labels)[end.cells] <-end.color
  blood.cells <- which(sce.hvg$labels %in%c(11,16))
  length(blood.cells)
  names(sce.hvg$labels)[blood.cells] <-bc.color
  
}


png(filename = file.path(results.dir, sample.name, paste0(sample.name, "_cellTypeClustering.2.png")), width = 600, height = 500)
par(mar = c(5,5,3,1),mgp=c(2.8, 1.2, 0))

plotSC(sce = sce.hvg,
       dimred = "TSNE",
       col =names(sce.hvg$labels),
       main = sample.name, 
       cex = 1,cex.axis = 2,lwd=4,cex.lab = 2 )
dev.off()


#### AddmoduleScore to select modules of cells according to expression ####
# AddModuleScore: A positive score would suggest that 
# this module of genes is expressed in a particular cell more highly 
# than would be expected,
# given the average expression of this module across the population

#### genes to represent modules
# Markers Roadmap
# gene.markers.nc <- gene.markers[names(gene.markers) == "NC"]
# gene.markers.nc <- gene.markers.nc[gene.markers.nc%in% rownames(sce.hvg)]
# df <- data.frame(Gene= gene.markers.nc, Stage=names(gene.markers.nc), Tumor = sample.name)

# gene.markers.nc <- c("SOX10","NGFR","SSUH2","PTX3","ARHGAP5","SOX6","SCD","CIT","GFRA1","RXRG")
# gene.markers.nc <- c("SOX10","NGFR") #classics

# gene.markers.d7 <- gene.markers[names(gene.markers) == "day7"]
# gene.markers.d7 <- gene.markers.d7[gene.markers.d7%in% rownames(sce.hvg)]
# df <- rbind(df,data.frame(Gene= gene.markers.d7, Stage=names(gene.markers.d7), Tumor = sample.name))
# gene.markers.d7 <- c("CDH19","CADM1", "EHBP1","SPP1", "PMEPA1", "ITGB8")
# gene.markers.d7 <- c("CDH19","ITGA4","MPZ")#classics

# gene.markers.d14 <- gene.markers[names(gene.markers) == "day14"]
# gene.markers.d14 <- gene.markers.d14[gene.markers.d14%in% rownames(sce.hvg)]
# df <- rbind(df,data.frame(Gene= gene.markers.d14, Stage=names(gene.markers.d14), Tumor = sample.name))

# gene.markers.d14 <- c("PRSS23")
# gene.markers.d14 <- c("TFAP2A","HNK1","EGR2")#classics

# gene.markers.d30<- gene.markers[names(gene.markers) == "day30"]
# gene.markers.d30 <- gene.markers.d30[gene.markers.d30%in% rownames(sce.hvg)]
# df <- rbind(df,data.frame(Gene= gene.markers.d30, Stage=names(gene.markers.d30), Tumor = sample.name))

# gene.markers.d30 <- c("S100B","PLP1","ERBB3","DAG1","FXYD1","GFRA3","LGI4","PDLIM4","MIA","FXYD3")
# gene.markers.d30 <- c("S100B","PLP1","DAG1", "PMP22","ERBB3")

# write.table(df,file.path(markers.roadmap.dir,paste0(sample.name,"_RoadmapMarkersScoring.csv")),sep = "\t", col.names = TRUE,row.names = FALSE)
# we select those roadmap markers expressed in all Tumors.
road.markers.tumors <- list()
for(i in seq_len(length(sample.names))){
  sn <- sample.names[i] 
  df <- read.table(file.path(results.dir, sn,"MarkersRoadmap",paste0(sn,"_RoadmapMarkersScoring.csv")), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  road.markers.tumors[[sn]] <- df[,1]
  names(road.markers.tumors[[sn]]) <- df[,2]
}
common.roadmap.genes <- Reduce(intersect,road.markers.tumors)

# ##################  sup table paper ###############################
# gene.markers <- gene.markers[names(gene.markers) !="PSC"]
# df.g <- data.frame(genes = gene.markers, stage = names(gene.markers))
# df.g
# i =1
# for(i in seq_len(length(road.markers.tumors))){
#   tm.nm <- names(road.markers.tumors)[i]
#   tm <- road.markers.tumors[[tm.nm]]
#   df.g[,tm.nm] <- ""
#   for(j in seq_len(length(df.g$genes))){
#     g <- df.g$genes[j]
#     pos <- g %in% tm
#     if(pos){
#       df.g[j,tm.nm] <- "expressed"
#       
#     }else{
#       df.g[j,tm.nm] <- ""
#     }
#   }
# }
# df.g$commonSingleCell <- ""
# i=13
# for(i in seq_len(nrow(df.g))){
#   # if(is.null(ncol(df.g[i,grepl("expressed",df.g[i,])])))next
#   if(ncol(data.matrix(df.g[i,grepl("expressed",df.g[i,])])) ==3){
#     df.g$commonSingleCell[i] <- "commonExp"
#   }
# }
# write.table(df.g, file.path(results.dir,"supp_file1.csv"),sep = "\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
# table(gene.markers%in% road.markers.tumors$PNF19)

gene.markers.nc  <- common.roadmap.genes[common.roadmap.genes%in% road.markers.tumors$PNF19[names(road.markers.tumors$PNF19) =="NC"]]
names(gene.markers.nc) <- rep("NC",length(gene.markers.nc))
gene.markers.d7 <- common.roadmap.genes[common.roadmap.genes%in% road.markers.tumors$PNF19[names(road.markers.tumors$PNF19) =="day7"]]
names(gene.markers.d7) <- rep("day7",length(gene.markers.d7))

gene.markers.d14 <- common.roadmap.genes[common.roadmap.genes%in% road.markers.tumors$PNF19[names(road.markers.tumors$PNF19) =="day14"]]
names(gene.markers.d14) <- rep("day14",length(gene.markers.d14))

gene.markers.d30 <- common.roadmap.genes[common.roadmap.genes%in% road.markers.tumors$PNF19[names(road.markers.tumors$PNF19) =="day30"]]
names(gene.markers.d30) <- rep("day30",length(gene.markers.d30))

df <- data.frame(genes = c(gene.markers.nc,gene.markers.d7,gene.markers.d14,gene.markers.d30), stages = names(c(gene.markers.nc,gene.markers.d7,gene.markers.d14,gene.markers.d30)))

cut.off <- 0 # cut off for selecting the module of cells.
# sc.cells <- colnames(sce.hvg)[which(sce.hvg$infomap %in% c(1,2,5,7,11,12,13,14))]
####################

gene.set.list <- list(nc= gene.markers.nc,
                      d7= c(gene.markers.d7),
                      # d14= gene.markers.d14,
                      d30= gene.markers.d30)

analysis.name <- "roadmapMarkers"
if(sample.name == "PNF19"){
  num.clusters <- c(1,2,5,7,11,12,13,14) #PNF19 SC_clust
}else if(sample.name == "PNF20"){
  num.clusters <- c(2,16,14,17,5,3,10,7)#PNF20
}else if(sample.name == "PNF23"){
  num.clusters <- c(2,5,7,6,3,4,10,17)#PNF23
}



specific.clusters <- num.clusters
coldata.clusters <- "infomap"
ctrl = 50
marrow <- AddMarkersModuleScore(sce=sce.hvg,
                                gene.set.list = gene.set.list,
                                analysis.name = "roadmapMarkers",
                                specific.clusters=num.clusters,
                                coldata.clusters = coldata.clusters,
                                ctrl =ctrl)

cells.to.plot <- getCellgroupsByGenesExp(sce = sce.hvg,
                                         gene.set.list  = gene.set.list,
                                         coldata.clusters ="infomap", 
                                         specific.clusters = num.clusters,
                                         analysis.name = "roadmapMarkers"
)

### Representation
stage.color <- c(NC ="#1832F5", day7 ="#1BB52F", day14 ="#8FEC0D", day30 ="#E9E333")
analysis.names <- names(cells.to.plot)[1:3] #no14days
names(analysis.names) <- c("NC","day7","day30") #no14days

#cell colors asignation
#giving grey colors for all cells
i = 1
for(i in seq_len(length(analysis.names))){
  #General colors
  cell.color<- colnames(sce.hvg)
  names(cell.color)<- colnames(sce.hvg)
  cell.color[sce.hvg$labels %in%num.clusters] <- karyoploteR::transparent("grey60") #In SC cluster
  cell.color[!sce.hvg$labels %in%num.clusters] <- "grey88" #Transparent not in SC
  
  
  stg <- names(analysis.names)[i]
  an <- analysis.names[i]
  stg.cells <- names(cells.to.plot[[an]])
  cell.color[names(cell.color) %in%stg.cells] <- transparent(stage.color[stg])
  png(filename = file.path(markers.roadmap.dir,paste0(sample.name, "_", stg,"_CellStage_SpecificColor_Scoring_Common_RoadmapMarkers_AllClusters_No14d.L2Norm.png")), width = 600, height = 500)
  par(mar = c(5,5,3,1),mgp=c(2.8, 1.2, 0))
  plotSC(sce = sce.hvg, 
         dimred = "TSNE",
         col = cell.color,
         # main = paste0(stg,"scoring"),
         cex = 1, cex.axis = 2,lwd=4,cex.lab = 2)
  dev.off()
}

#General colors
cell.color<- colnames(sce.hvg)
names(cell.color)<- colnames(sce.hvg)
cell.color[sce.hvg$labels %in%num.clusters] <- karyoploteR::transparent("grey60") #In SC cluster
cell.color[!sce.hvg$labels %in%num.clusters] <- "grey88" #Transparent not in SC

for(i in seq_len(length(analysis.names))){
  
  
  stg <- names(analysis.names)[i]
  an <- analysis.names[i]
  stg.cells <- names(cells.to.plot[[an]])
  cell.color[names(cell.color) %in%stg.cells] <- transparent(stage.color[stg])
  
  
}
png(filename = file.path(markers.roadmap.dir,paste0(sample.name, "_CellStages_Scoring_Common_RoadmapMarkers_uniquePlot_AllClusters_transparentColors_No14d.L2Norm.png")), width = 600, height = 500)
par(mar = c(5,5,3,1),mgp=c(2.8, 1.2, 0))
plotSC(sce = sce.hvg, 
       dimred = "TSNE",
       col = cell.color,
       # main = "Stage scoring",
       cex = 1, cex.axis = 2,lwd=4,cex.lab = 2)
dev.off()





# #########Individual marker plots ######################
# genes <- c("PDGFRA","IL33","CRLF1","PDGFA")
# genes <- c("ANGPT2", "JAG1")
# genes <- c("NOV","INHBA","FIGF")
# genes <- c("SFRP2")
# genes <- c("PECAM1")
# 
# genes <- c("SMOC2")
# 
# for(sn in sample.names) {
#   for(m in genes) {
#     if(!m %in% rownames(dd[[sn]]))next
#     tryCatch({
#       p <- plotMarkers(dd[[sn]], ncol=1, genes = m, title.text.size = 25) 
#       # save_plot(file.path(mkrs.dir,paste0(sn, "_singleGene_Clust",n.clust,"_",group,"_", m, ".png")), plot = p, base_height = 4, base_width = 4)
#       save_plot(file.path(results.dir,sn,"MarkerSingleCell",paste0(sn, "_singleGene_Clust", "_", m, ".png")), plot = p, base_height = 4, base_width = 4)
#       
#     }, error=function(e){message(e)})
#   }
# }  

################## CellCycle status using cyclone from Scran ########################
# hm.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
# 
# hm.pairs$G1
# head(hm.pairs$G1)
# hm.pairs$G1$first <- mapIds(edb, keys = hm.pairs$G1$first, keytype = "GENEID", column = "SYMBOL")
# hm.pairs$G1$second <- mapIds(edb, keys = hm.pairs$G1$second, keytype = "GENEID", column = "SYMBOL")
# hm.pairs$S$first <- mapIds(edb, keys = hm.pairs$S$first, keytype = "GENEID", column = "SYMBOL")
# hm.pairs$S$second <- mapIds(edb, keys = hm.pairs$S$second, keytype = "GENEID", column = "SYMBOL")
# hm.pairs$G2M$first <- mapIds(edb, keys = hm.pairs$G2M$first, keytype = "GENEID", column = "SYMBOL")
# hm.pairs$G2M$second <- mapIds(edb, keys = hm.pairs$G2M$second, keytype = "GENEID", column = "SYMBOL")
# gene.var <- modelGeneVar(sce)
# high.variance.genes <-  getTopHVGs(gene.var, prop=0.1)
# # sce.hvg <- sce[unique(c(high.variance.genes, genes.of.interest))]
# sce.hvgs <- sce[rownames(sce) %in%unique(c(high.variance.genes, gene.markers))]
# 
# assigned <- cyclone(sce.hvgs, pairs=hm.pairs)
# # save(assigned,file =file.path(cell.cycle.dir, "cycle_23.RData"))
# load(file.path(cell.cycle.dir, "cycle_19.RData"))
# load(file.path(cell.cycle.dir, "cycle_23.RData"))
# colLabels(sce.hvg) <- factor(assigned$phases)
# # Representation of cellcycle status obtained by cyclone
# png(filename = file.path(results.dir, paste0("clustering_TSNE_",sample.name,"_cellCycleStatus_cyclone.png")))
# plotReducedDim(sce.hvg, dimred="TSNE", colour_by = "label",text_by="label")
# dev.off()

# ########### Seurat Cell cycle status ################
# # Cycle Genes
# s <- cc.genes.updated.2019$s.genes
# g2 <- cc.genes.updated.2019$g2m.genes
# cell.cycle.gens <- c(s,g2)
# sce.cellclycle <- sce.hvg[rownames(sce.hvg)%in%unique(cell.cycle.gens)]
# cell.cycle.gens <- cell.cycle.gens[cell.cycle.gens%in%rownames(sce.hvg)]
# s <- s[s%in%cell.cycle.gens]
# g2 <- g2[g2 %in% cell.cycle.gens]
# 
# # Create our Seurat object and complete the initalization steps
# pp <- assays(sce.hvg)$counts
# # colnames(pp) <- paste0("cell_",1:ncol(pp))
# marrow <- CreateSeuratObject(counts = pp)
# marrow <- NormalizeData(marrow)
# marrow <- FindVariableFeatures(marrow, selection.method = "vst")
# marrow <- ScaleData(marrow, features = rownames(marrow))
# 
# marrow <- RunPCA(marrow, features = VariableFeatures(marrow), ndims.print = 6:10, nfeatures.print = 10)
# marrow <- RunTSNE(marrow, features = VariableFeatures(marrow), ndims.print = 6:10, nfeatures.print = 10)
# 
# # DimHeatmap(marrow, dims = c(1, 2))
# marrow <- CellCycleScoring(marrow,
#                            s.features = s,
#                            g2m.features = g2,
#                            set.ident = TRUE)
# 
# colLabels(sce.hvg) <- factor(marrow$Phase)
# 
# png(filename = file.path(cell.cycle.dir, paste0("cellcycle_status_",sample.name,"_Seurat.3.png")))
# plotReducedDim(sce.hvg, dimred="TSNE", colour_by = "label",text_by="label")
# dev.off()

# 
# ##################### Markers Selection  ###################
# sce.hvg.nomt <-sce.hvg[!grepl(rownames(sce.hvg),pattern = "^MT-"),]
# colLabels(sce.hvg)
# 
# if(sample.name == "PNF19"){
#   sc.cells <- c(1,2,5,7,11,12,13,14) #PNF19 SC_clust
# }else if(sample.name == "PNF20"){
#   sc.cells <- c(2,16,14,17,5,3,10,7)#PNF20
# }else if(sample.name == "PNF23"){
#   sc.cells <- c(2,5,7,6,3,4,10,17)#PNF23
# }
# 
# 
# sce.hvg.nomt <- sce.hvg.nomt[,sc.cells] #Markers SC
# colLabels(sce.hvg.nomt)  <- colData(sce.hvg.nomt)$infomap
# colLabels(sce.hvg.nomt) <- factor(as.character(colLabels(sce.hvg.nomt)), levels=sort(unique(as.character(colLabels(sce.hvg.nomt)))))
# 
# # markers <- findMarkers(sce.hvg.nomt, direction = "up", lfc =1)
# # markers <- findMarkers(sce.hvg.nomt, test= "wilcox", direction = "up")
# # markers <- findMarkers(sce.hvg.nomt, test="binom", direction="up")
# markers.all.test <- multiMarkerStats(t=findMarkers(sce.hvg.nomt, direction="up"),
#                                      wilcox=findMarkers(sce.hvg.nomt, test="wilcox", direction="up"),
#                                      binom=findMarkers(sce.hvg.nomt, test="binom", direction="up"))
# 
# chosen.clust <- "9" #the cluster we want to investigate
# # interesting <- markers[[chosen.clust]]
# interesting <-markers.all.test[[chosen.clust]] 
# 
# 
# 
# # interesting[rownames(interesting) == "S100B",]
# # interesting <- markers[[chosen.clust]][rownames(markers[[chosen.clust]]) %in%road.markers.d30,]
# best.set <- interesting[interesting$Top <= 10,]
# # logFCs <- getMarkerEffects(best.set)
# 
# # pairwise and binom
# prefix <- "t.logFC"
# prefix <- "binom.logFC"
# 
# logFCs <- getMarkerEffects(best.set,prefix = prefix)
# # png(filename = file.path(markers.scRNAseq,
# #                          paste0(sample.name,"_Markers_SCcluster_multitest_pheatmap_",prefix,"_clust.",chosen.clust,"_scran.png")),
# #     width=1000, height = 800)
# pheatmap(logFCs, 
#          fontsize = 10,
#          fontsize_row = 8,
#          breaks=seq(-5, 5, length.out=100),
#          show_rownames = TRUE,
#          main= paste0("Cluster ", chosen.clust, " ", prefix))
# # dev.off()
# 
# #Wilcoxon
# prefix <-"wilcox.AUC"
# logFCs <- getMarkerEffects(best.set,prefix = prefix)
# # png(filename = file.path(markers.scRNAseq,
# #                          paste0(sample.name,"_Markers_SCcluster_multitest_pheatmap_",prefix,"_clust.",chosen.clust,"_scran.png")),
# #     width=1000, height = 800)
# pheatmap(logFCs,
#          fontsize = 10,
#          fontsize_row = 8,
#          # color = color.plate,
#          breaks=seq(0,1, length.out=100),
#          show_rownames = TRUE,
#          main= paste0("Cluster ", chosen.clust, " ", prefix))
# # dev.off()
# 
# #Makres in scRNA-seq in the roadmap markers?
# cluster.group <- as.numeric(unique(colLabels(sce.hvg.nomt)))
# real.markers <- c()
# for(i in seq_len(length(cluster.group))){
#   chosen.clust <- as.character(cluster.group[i]) #the cluster we want to investigate
#   interesting <-markers.all.test[[chosen.clust]] 
#   # best.set <- interesting[interesting$Top <=10,]
#   # int <-gene.markers[gene.markers%in%rownames(interesting)]
#   int <-common.roadmap.genes[common.roadmap.genes%in%rownames(interesting)]
#   
#   # real.markers <- c(real.markers,gene.markers[gene.markers%in%rownames(best.set)])
#   real.markers <- c(real.markers,int)
#   
# }
# length(int)
# length(real.markers)
# real.markers <- unique(real.markers)
# real.markers.road <- names(gene.markers)
# names(real.markers.road)<- gene.markers
# real.markers <- real.markers.road[real.markers]
# # names(real.markers) <- names(gene.markers)
# unique(real.markers)
# real.markers <- real.markers[real.markers!= "PSC"]
# 
# df.new.markers <- data.frame(genes = names(real.markers),stages = real.markers)
# write.table(df.new.markers,file = file.path(markers.scRNAseq,paste0(sample.name,"_newMarkersRoadmap","_",sample.name,".csv")), sep = "\t", col.names = TRUE,row.names = FALSE)
# df.new.markers <-df
# # genes <- c("HDAC9","LMO4","PPP1R10","SCML1","SNAI1","SSUH2","STMN1","TOB1","GFRA1",
# #            "BIRC2", "EHBP1","COTL1","FRMD6","GPAT3","RCAN1","SERPINE2","USP53",
# #            "ABCA8","CCL2","COL6A1","DAG1","ENDOD1","FXYD1","FXYD3","GAS7","GFRA3","KLF9","LGI4","MIA","PLP1","PLPP1","S100B")
# # names(genes)<- c(rep("NC",9),rep("day7",8),rep("day30",15))
# # df.new.markers <- data.frame(genes = genes, stages =names(genes))
# # # df.new.markers <- data.frame(genes = c("SOX10","CDH19","S100B","PLP1","GAP43"))
# for(sm in sample.names){
#   
#   for(i in seq_len(nrow(df.new.markers))){
#     m <- df.new.markers$genes[i]
#     if(!m %in%rownames(dd[[sm]]))next
#     p <- plotMarkers(dd[[sm]], ncol=1, genes = m, title.text.size = 25, x.axis.size = 15, y.axis.size = 15, size.x.axis.line = 1, size.y.axis.line = 1, ticks.size = 1,axis.title.size = 15) 
#     save_plot(file.path(results.dir,sm,"MarkerSingleCell",paste0(sm, "_SCclusterMarkers_singleGene_CommonRoadmap_",df.new.markers$stages[i],"_", m, ".png")), plot = p, base_height = 4, base_width = 4)
#     
#     # save_plot(file.path(markers.scRNAseq,paste0(sample.name, "_singleGene_CommonRoadmap_",df.new.markers$stages[i],"_", m, ".png")), plot = p, base_height = 4, base_width = 4)
#     # save_plot(file.path(markers.scRNAseq,paste0(sample.name, "_singleGene_ClassicMarkers","_", m, ".png")), plot = p, base_height = 4, base_width = 4)
#     
#   }
#   
# }
# 
# dev.off()
# #Counts of genes by clusters
# top.genes <- head(rownames(interesting),10)
# 
# png(filename = file.path(markers.scRNAseq,
#                          paste0(sample.name,"_TopMarkers_multitest_Expression_clust.",chosen.clust,"_scran.png")),
#     width=1000, height = 800)
# 
# plotExpression(sce.hvg.nomt, x="label", features = top.genes, log2_values = 1, colour_by = "label")
# dev.off()

