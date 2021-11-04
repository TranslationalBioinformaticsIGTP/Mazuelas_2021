################################################################################
#                               RNA-Seq Pipeline                               #
################################################################################

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
source(file = "/imppc/labs/eslab/mmagallon/Projects/RNA-Seq-timecourse.2/Analysis/rna_seq_Functions.R")

############################  Parameters  ##########################
# # FiPS parameters
# params <- yaml.load_file("/imppc/labs/eslab/mmagallon/Projects/RNA-Seq-timecourse.2/Analysis/RNA-Seq_2D_3D_ipsDifferentiation/Parameters/FiPS_pipeline_parameters.yaml")
## NF1 parameters
# params <- yaml.load_file("/imppc/labs/eslab/mmagallon/Projects/RNA-Seq-timecourse.2/Analysis/RNA-Seq_2D_3D_ipsDifferentiation/Parameters/NF1_2D_pipeline_parameters.yaml")
# # NF1_A 5MM parameters
# sn <- "NF1_A"
# # # NF1_B 3MM parameters
# # sn <- "NF1_B"
# NF1_C D12 parameters
# sn <- "NF1_C"
# # FiPS NF1 2D parameters
# params <- yaml.load_file("/imppc/labs/eslab/mmagallon/Projects/RNA-Seq-timecourse.2/Analysis/RNA-Seq_2D_3D_ipsDifferentiation/Parameters/FiPS_NF1_2D_pipeline_parameters.yaml")
# # FiPS + SC parameters
# params <- yaml.load_file("/imppc/labs/eslab/mmagallon/Projects/RNA-Seq-timecourse.2/Analysis/RNA-Seq_2D_3D_ipsDifferentiation/Parameters/FiPS_SC_pipeline_parameters.yaml")
# # Sphers hom vs het
# params <- yaml.load_file("/imppc/labs/eslab/mmagallon/Projects/RNA-Seq-timecourse.2/Analysis/RNA-Seq_2D_3D_ipsDifferentiation/Parameters/Spherers_pipeline_paramaters.yaml")
# # FiPS + Spheres parameters
# params <- yaml.load_file("/imppc/labs/eslab/mmagallon/Projects/RNA-Seq-timecourse.2/Analysis/RNA-Seq_2D_3D_ipsDifferentiation/Parameters/FiPS_2Dvs3D_pipeline_parameters.yaml")
#All samples parameters
params <- yaml.load_file("/imppc/labs/eslab/mmagallon/Projects/RNA-Seq-timecourse.2/Analysis/RNA-Seq_2D_3D_ipsDifferentiation/Parameters/AllSamples_pipeline_parameters.yaml")
# #FB vs SC parameters
# params <- yaml.load_file("/imppc/labs/eslab/mmagallon/Projects/RNA-Seq-timecourse.2/Analysis/RNA-Seq_2D_3D_ipsDifferentiation/Parameters/SCvsFB_pnf_pipeline_parameters.yaml")


#Model of Analysis
model <- params$model
samples.group <-params$samples.group
stages <- params$stages

#### loading the sample data information 
sample.data <- read.table(file = params$sample.data.file , header = T, sep = "\t", stringsAsFactors = FALSE, comment.char = "")
sample.data <- sample.data[sample.data$Genotype %in% c(params$genotype),]
sample.data <- sample.data[!sample.data$graph.Names %in% c( "", "FB_E"),]
sample.data$graph.Names <- gsub("FB","Fb",sample.data$graph.Names)

# NF1 samples alone
# sample.data <- sample.data[grepl(sn,sample.data$graph.Names),]

# Getting the file names of the analysis
file.names <- sample.data$File.Name
sample.names <- sample.data$graph.Names
goterms <- params$go.terms

# Salmon alignement and quantification parameters
# salmonDir <- "/soft/bio/salmon-1.1.0/bin/salmon"
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

# Genes FiPS pipelina
gene.markers <- c() 
stages <- c("PSC", "NC", "day7", "day14", "day30")
for(i in seq_len(length(stages))){
  gr <- stages[i]
  mks<- read.table(file = file.path('/imppc/labs/eslab/mmagallon/Projects/RNA-Seq-timecourse.2/Analysis/RNA-Seq_2D_3D_ipsDifferentiation/Results/2D/markers', paste0("up_", gr,".txt")), header = FALSE,stringsAsFactors = FALSE)[,1]
  # mks <- markers[[gr]]
  names(mks) <- rep(gr, length(mks))
  gene.markers <- c(gene.markers,mks)
  
}
#Deleting PSC markers
gene.markers <- gene.markers[names(gene.markers) != "PSC"]
# # allDEG
# gene.markers <- c()
# stages <-  c("PSC", "NC","day7", "day14", "day30")
# for(i in seq_len(length(stages)-1)){
#   gr <- stages[i]
#   for(j in (i+1):(length(stages))){
#     gr2 <- stages[j]
#     markers<- read.table(file = file.path(selected.genes.dir, paste0("Genes_", gr,"vs",gr2,".csv")),header = FALSE,stringsAsFactors = FALSE)[,1]
#     names(markers) <- rep(gr, length(markers))
#     gene.markers <- c(gene.markers,markers)
#     gene.markers <- unique(gene.markers)
#   }
#   
# }


######### Directories
model.dir <- params$model.dir

pca.dir <-  file.path("Results", model.dir, "PCAplot")
if(!file.exists(pca.dir)) dir.create(pca.dir)

distance.dir <- file.path("Results", model.dir, "DistancePlots")
if(!file.exists(distance.dir)) dir.create(distance.dir)

heatmap.QC.dir <- file.path("Results", "heatmap.QC")
if(!file.exists(heatmap.QC.dir)) dir.create(heatmap.QC.dir)

heatmap.dir <- file.path("Results",model.dir, "DEG_heatmap")
if(!file.exists(heatmap.dir)) dir.create(heatmap.dir)

DEG.dir <- file.path("Results", model.dir, "DESeqResults")
if(!file.exists(DEG.dir)) dir.create(DEG.dir)

selected.genes.dir <- file.path("Results", model.dir, "DEGgenes")
if(!file.exists(selected.genes.dir)) dir.create(selected.genes.dir)

markers.dir <- file.path("Results", model.dir, "markers")
if(!file.exists(markers.dir)) dir.create(markers.dir)

gsea.dir <- file.path("Results", model.dir, "GSEAGO")
if(!file.exists(gsea.dir)) dir.create(gsea.dir)
for(term in goterms){
  if(!file.exists(file.path(gsea.dir, term))) dir.create(file.path(gsea.dir, term))
}

revigo.dir <- file.path("Results", model.dir, "REViGO")
if(!file.exists(revigo.dir)) dir.create(revigo.dir)

annot.dir <- file.path("Results",model.dir,  "Annotation")
if(!file.exists(annot.dir)) dir.create(annot.dir)

gene.expression.dir <- file.path("Results",model.dir,  "GeneExpressionPlot")
if(!file.exists(gene.expression.dir)) dir.create(gene.expression.dir)


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

# We import the data from tximport output to analyse it using DESeq2
dds.txi <- DESeqDataSetFromTximport(txi = txi.salmon, colData = sample.data, formula(paste0("~", samples.group)))
dim(dds.txi)

# Estimation of Size Factors
dds <- estimateSizeFactors(dds.txi)

# Filtering rawdata: removing data without expression
keep <- rowSums(counts(dds)>=filt.min.reads)>=filt.min.samples
dds <- dds[keep,]
#Normalization
dds.deg <- DESeq(object = dds, test = "Wald")
dds.rlog <- rlog(dds.deg)

# Taking a look after filtering and normalizing
# boxplot(filtered.pseudocount, main = "BoxPlot_AllSamples_Filtered")

# heatmap QC

# mat.dist <- counts(dds.deg,normalize =T)[gene.markers,]
# mat.dist <- assay(dds.rlog)

# mat.dist <- as.matrix(dist(t(mat.dist), method = "manhattan"))
# mat.dist <- mat.dist/max(mat.dist)

# 
# 
# png(file.path(heatmap.QC.dir, paste0( model,"SCvsFb_QCHeatMap_manhattan_dist2_.png")), width =1000, height = 1000)
# pheatmap(mat.dist, color = color.plate, fontsize = 5, legend = TRUE, treeheight_row = 50, treeheight_col = 50, lwd = 5, border_color = NA)
# # cim(mat_dist, row.cex = 10, col.cex = 10, center = 24, margins = c(100,100), keysize = c(0.5,0.5), color = color.plate)
# dev.off()

###### Dendogram
# 
# # mat.dist <- counts(dds.deg,normalize =T)[gene.markers,]
# mat.dist <- assay(dds.rlog)[gene.markers,]
# mat.dist <- mat.dist[,!grepl("PSC",colnames(mat.dist))]
# 
# mat.dist <- mat.dist[,c(1:3,7:9,13:15,19:21,25:38)]
# 
# data.clust <- dist(t(mat.dist), method = "euclidean")
# cluster <- hclust(data.clust)
# cluster <- as.dendrogram(cluster)
# # Define nodePar
# # nodePar <- list(lab.cex = 0.6, pch = c(NA, 15), 
# #                 cex = 0.7, col = c("blue"))
# library(dendextend)
# 
# png(file.path(heatmap.QC.dir,
#               paste0( model,"_FiPS_3D_cells_DendogramMarkers_Euclidean_dist.png")), 
#     width =1500, height = 1000)
# par(mar =c(30,10,10,0),oma =c(0,0,0,0),mgp=c(6,1,0))
# cluster %>% 
#   set("labels_cex", 3)%>%
#   set("leaves_pch", 15)%>%
#   set("leaves_cex", 5)%>%
#   # #Euclidean All Samples
#   # set("leaves_col", c("white","white","white","white","white",
#   #                     "#E9E333","#E9E333","#E9E333",
#   #                     "#8FEC0D","#8FEC0D","#8FEC0D",
#   #                     "#1BB52F","#1BB52F","#1BB52F",
#   #                     "white","white","white","white",
#   #                     "#1BB52F","#8FEC0D","#E9E333","#8FEC0D",
#   #                     "white","white","white",
#   #                     "#1BB52F","#E9E333",
#   #                     "#1832F5","#1832F5","#1832F5","#1832F5","#1832F5","#1832F5",
#   #                     "#8FEC0D","#E9E333","#1BB52F","white","white"))%>%
# #Euclidean fips+3d+sc+fb
# set("leaves_col", c("#1832F5","#1832F5","#1832F5",
#                     "white","white","white","white",
#                     "white","white","white","white","white",
#                     "white","white","white","white","white",
#                     "#E9E333","#E9E333","#E9E333",
#                     "#8FEC0D","#8FEC0D","#8FEC0D",
#                     "#1BB52F","#1BB52F","#1BB52F"))%>%
#   
# 
# plot(ylab ="Height",
#        cex.lab = 5,
#        cex.main=10,
#        edgePar = list(col = "black", lwd = 10),
#        yaxt ="n")
# axis(side = 2, lwd = 10,at =c(0,40,80),cex.axis = 5)
# dev.off()


#### PCA
# pca.pseudocounts <-rlog(dds.deg)
pca.pseudocounts <- dds.rlog
# cell.type <- sample.data$Diff.Day
# sample.data$Diff.Day <- gsub(x = sample.data$Diff.Day, pattern = "-", replacement = "SC")

# PCA with just marker genes
pca.pseudocounts <- pca.pseudocounts[rownames(pca.pseudocounts) %in% gene.markers[!names(gene.markers)%in% "PSC"],]

# cell.type <- sample.data$Spheresvs2D
colData(pca.pseudocounts)$Diff.Day <- factor(as.character(colData(pca.pseudocounts)$Diff.Day), levels=c("PSC","NC","day7","day14","day30"))
Diff.Day <- colData(pca.pseudocounts)$Diff.Day


colData(pca.pseudocounts)$Cell.Type[sample.data$Sphere.Type %in% "Heterotypic" ] <- "Sphere"
colData(pca.pseudocounts)$Cell.Type[sample.data$Sphere.Type %in% "Homotypic" ] <- "Sphere"

colData(pca.pseudocounts)$Cell.Type <- factor(as.character(colData(pca.pseudocounts)$Cell.Type), levels=unique(colData(pca.pseudocounts)$Cell.Type))
Samples <- colData(pca.pseudocounts)$Cell.Type
                             
# mycolors <- c(PSC = "orange",SC = "skyblue", NC="#F8766D", day7 = "#00BF7D", day14 = "#A3A500", day30 = "purple1")
mycolors <- unique(sample.data$colors.p1)
names(mycolors) <- unique(colData(pca.pseudocounts)$Diff.Day)

sn <- sample.data$graph.Names
# sn[grepl(pattern = "WT",x=sn)] <- ""

png(filename=file.path(pca.dir,paste0(model, "_FiPS_MM_iPSC_NOlabel_newshape_pca.p1.png")), width = 1500, height = 800)
plotPCA(pca.pseudocounts, intgroup = "Diff.Day" ) +  
  
  # To change the fucking color
  scale_colour_manual(name = "Stages", values = mycolors)+
  # geom_text(aes(label = as.character(sn)), hjust=-0.21 , vjust=0.05, size = 10 )+
  # ggrepel::geom_text_repel(mapping = aes(label = as.character(sn),lwd = 2), hjust = -1, vjust = 5,
  # segment.size = 0.1, nudge_x = -0.9, size =10, direction ="x")+
 
  geom_point(aes(shape= Samples, group = "Cell.Type"),fill = sample.data$colors.p1, size = 20, stroke=5)+
  
  
  #Scale_shape used to give the specific forms to the groups
  scale_shape_manual(name= "Genotype",values = c(23))+ # rombo
  # scale_shape_manual(values = c(15,17))+ # square, triangle
  # scale_shape_manual(values = c(15,18,16,17))+ # square, rombo, circle, triangle
  # scale_shape_manual(values = c(15,18,16))+ # square, circle, rombo
  # scale_shape_manual(values = c(15,16,17))+ # square, circle, triangle
  # scale_shape_manual(values = c(15,16,21))+ # square, circle, empty circle
  # scale_shape_manual(values = c(15,16))+ # square, circle
  # scale_shape_manual(values = c(15,18))+ # square, rombo, triangle
  # scale_shape_manual(values = c(18,9))+
  # ggrepel::geom_text_repel(mapping = aes(label = as.character(sn), lwd = 2),
  #                          size =10, direction ="y", vjust=0.8, hjust = 1.3, color = "black")+

  theme(legend.title = element_text(size=30, 
                                    face="bold"), 
        legend.text = element_text( size=30, face="bold"), 
        axis.text=element_text(size=30),
        axis.title=element_text(size=30,face="bold"),
        plot.title = element_text(size=30),
        panel.background = element_rect(fill="white"), 
        panel.grid.major = element_line(colour = "grey"), plot.margin =margin(20,10,10,10))+
  guides(colour = guide_legend(override.aes = list(shape = 15)))
  #+xlim(-70,20)#+ ylim(-30,50)

 
dev.off()
#################     Differential gene expression      ########################

  msg(paste("Starting DESeq of", samples.group, "contrast"))
  
#Filtering the data
filtered.dds <- getFilteredDDS(tximport = txi.salmon, samples_group = samples.group, samples_df = sample.data,filter_min_reads = filt.min.reads,filter_min_samples = filt.min.samples)
# stages <- unique(sample.data$graph.Names)
stages <- params$stages
#Getting the Deseq results

  deseq.results <- list()
  for(i in seq_len(length(stages)-1)) {
    # for(i in seq_len(length(con1))) {
    ref <- stages[i]
    # ref <- con1[i]
    
    deseq.results[[ref]] <- list()
    
    for(j in (i+1):length(stages)) {
    # for(j in 1:length(con2)) {
      other_group <- stages[j]
      # other_group <- con2[j]
      
      # FIPS samples and time %in% c(ref, otro)
      t <- c(ref, other_group)
      deseq_samp <- sample.data[sample.data[,samples.group] %in% t, c("Sample.Name", samples.group)]
      
      #runDESeq function
      
      message("RunDeSEQ ", other_group, " vs ", ref)
      
      deseq.results[[ref]][[other_group]] <-runDESeq( sample_names = deseq_samp$Sample.Name, sname_variable = "Sample.Name", ref_group = ref, samples_df = sample.data, filtered_dds = filtered.dds, type_lfcShrink = "apeglm", test = "Wald")
      
      msg("runDESeq done")
      
    }
    
  }
  
#Saving DESeqResults
# save(deseq.results, file = file.path(DEG.dir, paste0(model,"_DESeqResults_LFCshrinkApeglm.RData")))
load(file = file.path(DEG.dir, paste0(model,"_DESeqResults_LFCshrinkApeglm.RData")))

# Saving the DEG
selected_genes <- list()

for( i in 1:length(deseq.results)){
  lists <- deseq.results[[i]]
  class(lists)
  ref <- names(deseq.results[i])
  contrast <- data.frame()
  for (j in 1:length(lists)){
    contrast <- lists[[j]]
    other_group <- names(lists[j])
    
    # We select those genes which have a pvale <= 0.05 from all analised genes.
    selected_genes[[ref]][[other_group]]<- getTopGenesLFC(deseq_results = contrast, padj_value = 0.05)
    
    selected_genes[[ref]][[other_group]]$genes <- rownames(selected_genes[[ref]][[other_group]])
    
    # We save the results obtained
    # write.table(selected_genes[[ref]][[other_group]],
    #             file.path(selected.genes.dir, paste0("ResultTable_", ref,"vs", other_group,".csv")),
    #             sep = "\t", dec = ".", quote = FALSE, col.names = TRUE, row.names = FALSE)
    # 
    # 
    # write.table(selected_genes[[ref]][[other_group]]$genes,
    #             file.path(selected.genes.dir, paste0("Genes_", ref,"vs", other_group,".csv")),
    #             sep = "\t", col.names = FALSE, row.names = FALSE)

  }
}



###############          Specific Markers               #######################
updown_markers <- "up"
markers <- list()
annotation <- list()

for (i in seq_len(length(updown_markers))){
  
  cond <- updown_markers[i]
  j=1
  for(j in 1:length(stages)){
    gr<- stages[j]
    markers[[gr]]<-getGroupSpecificGenes(select_genes = selected_genes, query_group =gr, updown_markers = cond)
    write.table(markers[[gr]], file.path(markers.dir, paste0(cond,"_", gr,".txt")), col.names = FALSE, row.names = FALSE, quote = F)

    # Annotation of the putative stage specific markers
    # annotation <- getBioMartGOAnnotation(genes = markers[[gr]])
    # 
    # annot <- data.frame()

    # for(k in seq_len(length(annotation))){
    #   annot <- as.data.frame(annotation[k])
    #   bioterm <- names(annotation[k])

      # if(cond == "all"){
      # 
      #   #colum generation to get in files containing up and down genes which one is each one.
      #   annot$stages <- annot[,1]
      #   annot$stages[annot$stages %in% markers$up[[gr]]] <- "up"
      #   annot$stages[annot$stages %in% markers$down[[gr]]] <- "down"
      # }
        # write.table(annot, file.path(annot.dir, paste0(gr, "_", cond, "_", bioterm,"_markers_annotation.csv")),
        # col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


    # }
  }
}



#################     Enrichment specific markers    #########################
updown_markers <- "up"
# genes.markers <- selected_genes$Het_day14$day14_2D_PP$genes[selected_genes$Het_day14$day14_2D_PP$log2FoldChange <0]

# Markers Enrichment 
for (i in seq_len(length(updown_markers))){
  
  cond <- updown_markers[i]
  j=2
  for(j in 1:length(stages)){
    gr<- stages[j]
    genes.markers <- markers[[gr]]
    enrich <- getAllEnrichTermList(genes.markers)
    
    # Saving enriched terms
  
    for(term in goterms){
      # revigo.data <- enrich[[term]][, c("ID", "pvalue")]
      # revigo.contrast.dir <- file.path(revigo.dir, paste0("markers_", gr))
      # write.table(revigo.data, file.path(revigo.dir, paste0("RevigoTable_", term,"_markers_up_", model.dir,".csv")), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
      # write.table(revigo.data, file.path(revigo.dir, paste0("RevigoTable_", term,"Markers_",model,"_", model.dir,"_FiPS.csv")), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
      
      
      # write.table(data.frame(enrich[[term]]), file.path(gsea.dir, term, paste0("Enrichment_", term, "_markers_up_", model.dir, "_FiPS",".csv")),
      #             col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

      write.table(data.frame(enrich[[term]]), file.path(gsea.dir, term, paste0("Enrichment", term, "Markers_",model,"_", gr, ".csv")),
                  col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

    }
  } 
}



###############################################################################
#                          Result's Representation                            #
###############################################################################

##########################################################################################
##################################### Markers FiPS plots  ################################  

#### Heatmap Annotation

#Creating annotation.row data to plot
# annotation.row <- data.frame(Stages =gene.markers)
# rownames(annotation.row) <- annotation.row$Stages 
# annotation.row$Stages <- names(gene.markers)
# colnames(dds.rlog)

#Annotation Column
# #WT
# Diff.Days <- factor(c("PSC",
#                       "NC","NC","NC",
#                       "7d","7d","7d",
#                       "14d","14d","14d",
#                       "30d","30d","30d",
#                       "SC", "SC", "SC"), levels = c("PSC","NC","7d","14d","30d","SC"))
# colData(dds.rlog)$Diff.Day <-Diff.Days

#NF1 and WT
Diff.Days <- factor(c("PSC","PSC","PSC","PSC",
                      "NC","NC","NC","NC","NC","NC",
                      "7d","7d","7d","7d","7d","7d",
                      "14d","14d","14d","14d","14d","14d",
                      "30d","30d","30d","30d","30d","30d"), levels = c("PSC","NC","7d","14d","30d"))
colData(dds.rlog)$Diff.Day <-Diff.Days

#Sphers and WT
# Diff.Days <- factor(c("PSC",
#                       "NC","NC","NC",
#                       "7d","7d","7d",
#                       "14d","14d","14d",
#                       "30d","30d","30d","WT_Heterotypic","WT_Heterotypic","WT_Heterotypic",
#                       "NF1_Heterotypic","NF1_Heterotypic","NF1_Heterotypic",
#                       "NF1_Homotypic","NF1_Homotypic","NF1_Homotypic"),
#                     levels = c("PSC","NC","7d","14d","30d","WT_Heterotypic","NF1_Heterotypic","NF1_Homotypic"))
# colData(dds.rlog)$Diff.Day <-Diff.Days

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

# sample_grououp <- sample_grououp[-c(1,6,7,8)]
# sample_grououp <- c(sample_grououp,"white")
names(sample_grououp)
# sample_grououp[5:7] <- "white"
# names(sample_grououp)[5:7] <- ""
# sample_grououp <- sample_grououp[-c(6,7)]

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
#columncolors
# mycolors.col <- c("grey44","grey79")
# names(mycolors.col) <- c("WT", "NF1(-/-)")
mycolors.col <- c("grey44","grey70","grey89")
names(mycolors.col) <- c("WT", "NF1(+/-)","NF1(-/-)" )
mycolors <- list("Samples   " =sample_grououp, "Genotype   " = mycolors.col)


# No PSC
annotcol <- annotcol[-which(grepl("PSC",rownames(annotcol))),]
# annotcol[,"Samples   "][13:nrow(annotcol)] <- ""
annotcol[,"Samples   "][25:nrow(annotcol)] <- ""

# annotcol[-which(grepl("Fb",annotcol$`Genotype   `)),]
# annotcol <- annotcol[-which(grepl("X",rownames(annotcol))),]
# 
# rn <- rownames(annotcol)
# annotcol <- data.matrix(annotcol[-which(grepl("PSC",rownames(annotcol))),])
# rownames(annotcol) <- rn[-1]
# colnames(annotcol) <-"Samples   "
# annotcol<-data.frame(annotcol)
# colnames(annotcol) <-"Samples   "

#Representation of WT_2D and WT and NF1 3D spheres
# col.selected <- colnames(dds.rlog)
# col.selected <- c(1:16,18,19,17,21,22,20) #Fips+Spheres
col.selected <- c(1,3,4,2,5:7,9,10,8,11:13,15,16,14,17:19,21,22,20,23:25,27,28,26,29:31,33,34,32,36,37,35,40:42,38:39) #Fips+Spheres+SC+FB

# col.selected <- c(1,3,4,2,5:7,9,10,8,11:13,15,16,14,17:19,21,22,20,23:25,27,28,26)#FiPS+NF1_2D
# col.selected <- c(1,3,4,2,5:7,9,10,8,11:13,15,16,14,17:19,21,22,20,23:25,27,28,26,29:31,33,34,32,36,37,35,42:44)#AllSamples

# dds.rlog <- dds.rlog[,col.selected]
# dds.rlog <- rlog(filtered.dds)

#Selecting markers genes 
# count.data <- assay(dds.rlog)[rownames(annotation.row),]
# count.data <- assay(dds.rlog)[rownames(dds.rlog)%in%gene.markers,]
count.data <- assay(dds.rlog)[gene.markers,]
colnames(count.data)
# Data normalized to plot in the heatmap
data_subset_norm <- t(apply(count.data, 1, cal_z_score))
data_subset_norm <- data_subset_norm[,col.selected]
colnames(data_subset_norm)
#without iPSC
data_subset_norm <- data_subset_norm[,!grepl(x=colnames(data_subset_norm),pattern = "PSC")]
zscore <- data_subset_norm
heat.col <- colorRampPalette(c("blue","white","red"))(75)
myBreaks <- c(seq(min(zscore), 0, length.out=ceiling(75/2) + 1), 
              seq(max(zscore)/75, max(zscore), length.out=floor(75/2)))


#WT+sc
data_subset_norm <- data_subset_norm[,c(1:3,7:9,13:15,19:21,34:36)]
colnames(data_subset_norm)
# annotcol.f <- rbind(annotcol[grepl("WT",rownames(annotcol)),],annotcol[grepl("^SC",rownames(annotcol)),])[-c(13:15),]
annotcol.f <- annotcol[c(1:3,7:9,13:15,19:21,36:38),]
annotcol <- annotcol.f
annotcol[grepl("^SC",rownames(annotcol)),2] <- "PNF SC"
annotcol[grepl("^SC",rownames(annotcol)),1] <- "NF1(-/-)"

sample.g <- sample_grououp
sample.g[length(sample.g)]<-"#F0CE77"
# sample.g[length(sample.g)]<-"white"
names(sample.g)[length(sample.g)] <- "PNF SC"
colnames(annotcol) <- c("Genotype   ","Samples   ")
mycolors <- list("Samples   " =sample.g, "Genotype   " = mycolors.col[-2])

#WT+NF1
# data_subset_norm <- data_subset_norm[,c(1:24)]
# colnames(data_subset_norm)
# annot.nf <- annotcol[c(1:24),]
# sg <- sample_grououp
# sg <- sg[-length(sg)]
# gn <- mycolors.col
# gn <- gn[-2]
# mycolors <- list("Samples   " =sample_grououp, "Genotype   " = mycolors.col)
# annotcol <- annot.nf

#WT+spheres+sc+fb
annotcol[grepl("^SC",rownames(annotcol)),2] <- "PNF SC"
annotcol[grepl("^Fb",rownames(annotcol)),2] <- "PNF Fb"
sample.g <- sample_grououp
# sample.g[length(sample.g)]<-"white"
names(sample.g)[length(sample.g)] <- "PNF SC"
sample.g <- c(sample.g,"")
sample.g[length(sample.g)]<-"white"
names(sample.g)[length(sample.g)] <- "PNF Fb"
sample.g <- c(sample.g,"")
sample.g[length(sample.g)]<-"white"
names(sample.g)[length(sample.g)] <- ""

colnames(annotcol) <- c("Genotype   ","Samples   ")
mycolors <- list("Samples   " =sample.g, "Genotype   " = mycolors.col)

data_subset_norm <- data_subset_norm[,c(1:3,7:9,13:15,19:21,25:38)]
colnames(data_subset_norm)


#heatmap of stage specific markers

# png(filename = file.path(file.path(heatmap.dir, paste0("FiPS_AllSpheres_supervised_","Markers_NoLegend_", model, ".colBar.2.png"))), width = 1000, height = 800)
# png(filename = file.path(file.path(heatmap.dir, paste0("FiPS_supervised_","Markers_NoLegend_", model, ".png"))), width = 1000, height = 800)
png(filename = file.path(file.path(heatmap.dir, paste0("FiPS_Markers_AllSamples_supervised_", model, "_3verticallLegend_colBar.2.png"))), width = 1000, height = 800)
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
         gaps_col = c(12,15,18,21,24)
         # gaps_col= c(24,27,30,33)
         
         # main = paste0("Clustering all DEG in all samples", height," height ", nrow(annotation_row), " genes"))
         # main = c("Stage Specific Markers in FiPS")
)
dev.off()


#############################    Heatmaps   ################################


gene.list <- c(list(deg.genes = deg.genes, mesenchime.dev.genes= mesenchime.dev.genes), markers)
#CD105 = "ENG"
#CD73 = "NT5E" meri
#CD90 = "THY1"
#CD13 = "ANPEP" meri
#CD29 = "ITGB1"
#CD44 = "CD44" meri
#CD10 = "MME"
#CD45 = "PTPRC" NO esta
#CD34 = "CD34" meri
#CD14 = CD14
#CD11B = "ITGAM"
#CD79A = "CD79A"
#CD19 = "CD19"

#mesenchimal genes
genes <- c("ENG", "NT5E", "THY1","ANPEP", "ITGB1", "CD44","MME", "PTPRC","CD34", "CD14","ITGAM", "CD79A", "CD19")
genes <- c("SIX1","FOXF1", "WT1", "HMGA2", "BASP1", "FOXD1", "PITX2", "GATA4", "ISL1", "NKX2-5","SEMA3A", "TWIST1", "TBX20", "THY1", "CD44", "MME","ANPEP")

# NCtoSC linage genes
genes <- c("POU5F1","POU3F1","NGFR","SOX10", "CDH19", "ITGA4","PLP1", "GAP43", "S100B", "EGR2", "PMP22","MPZ", "JUN", "ADGRG6", "HDAC4")
genes <- gene.markers
deg.genes[deg.genes %in% genes]
# deg.genes <- selected_genes$Fb$

genes <- deg.genes
#rlog to estandarize the diferences.
dds.rlog <- rlog(filtered.dds)
dds.rlog <- rlog(dds.deg)
gene.list <- markers
gl <- "3D_DEGenes"
gl <- "SC_differentiation"
gl <- "mesenchimal_genes"

for(gl in names(gene.list)){
  genes <- gene.list[[gl]]
  
  counts.data <- assay(dds.rlog)[rownames(dds.rlog) %in% genes,]
  data_subset_norm <- t(apply(counts.data, 1, cal_z_score))
  
  if(gl == "deg.genes"){
    png(filename = file.path(file.path(heatmap.dir, paste0(model,"_", gl, "_DEG_heatmap.png"))),
        width = 5000, height = 3500)
    pheatmap(data_subset_norm, color = color.plate, fontsize = 50, legend = TRUE,
             fontsize_number = 8 ,treeheight_row = 250, treeheight_col = 250, lwd = 10, show_rownames = TRUE, cutree_rows = 2, fontsize_row = 50)
    dev.off()
  }
  if (gl == "markers"){
    png(filename = file.path(file.path(markers.dir, paste0(gl, "_markers_heatmap.png"))),
        width = 5000, height = 3500)
    pheatmap(data_subset_norm, color = color.plate, fontsize = 50, legend = TRUE,
             fontsize_number = 8 ,treeheight_row = 250, treeheight_col = 250, lwd = 10, show_rownames = TRUE, fontsize_row = 10)
    dev.off()
  }
  if(gl == "mesenchime.dev.genes"){
    png(filename = file.path(file.path(markers.dir, paste0(gl, "_mesenchime_heatmap.png"))),
        width = 3000, height = 2500)
    pheatmap(data_subset_norm, color = color.plate, fontsize = 60, legend = TRUE,
             fontsize_number = 1 ,treeheight_row = 250, treeheight_col = 250, lwd = 10, show_rownames = TRUE,legend_breaks = c(1.5, 0, -1.5))
    dev.off()
  }
}

######################## Gene expression Plot #############################
genes.lineage <- c("NGFR","SOX10", "EGR2","CDH19", "ITGA4","MPZ","GAP43","PLP1", "PMP22", "S100B")
names(genes.lineage) <- c("NC","NC","NC", "SCP","SCP","SCP","SCP","iSC","iSC", "SC")
genes <- c("NGFR","SOX10", "EGR2","CDH19", "ITGA4","MPZ","GAP43","PLP1", "PMP22", "S100B")
colData(dds.deg)$Diff.Day <- factor(as.character(colData(dds.deg)$Diff.Day), levels=c("PSC","NC","day7","day14","day30", "SC"))

########Classic and new stage markers#####################
gene.markers.nc <- c("SOX10","NGFR","SSUH2","PTX3","ARHGAP5","SOX6","SCD","CIT","GFRA1","RXRG")

# gene.markers.d7 <- gene.markers[names(gene.markers) == "day7"]
gene.markers.d7 <- c("CDH19","CADM1", "EHBP1","SPP1", "PMEPA1", "ITGB8")

# gene.markers.d14 <- gene.markers[names(gene.markers) == "day14"]
gene.markers.d14 <- c("PRSS23")

# gene.markers.d30<- gene.markers[names(gene.markers) == "day30"]
gene.markers.d30 <- c("S100B","PLP1","ERBB3","DAG1","FXYD1","GFRA3","LGI4","PDLIM4","MIA","FXYD3")
genes <- c(gene.markers.nc,gene.markers.d7,gene.markers.d14,gene.markers.d30)
genes <- c("EHBP1","GAS7","FXYD1","GFRA3")
filtered.dds <- getFilteredDDS(tximport = txi.salmon, samples_group = samples.group, samples_df = sample.data,filter_min_reads = filt.min.reads,filter_min_samples = filt.min.samples)

# Boxplot for gene expression
filt.counts <- counts(filtered.dds, normalize = TRUE)
# filtered.pseudocount <- rlog(counts(dds.deg))
colnames(filt.counts) <- c("PSC","NC","NC", "NC","day7","day7","day7","day14","day14","day14","day30","day30","day30")

filt.counts <- filt.counts[genes,]
# stages <- colnames(filt.counts)
stages <- colData(filtered.dds)$Diff.Day
fit.counts.t <- as.data.frame(t(filt.counts))

fit.counts.t <- cbind(fit.counts.t, stages)

numbers <- c( 0, 7, 14, 30)
names(numbers) <- c("NC", "day7","day14","day30")
i=1
fit.counts.t$stages <- factor(fit.counts.t$stages, levels = c("PSC","NC", "day7","day14","day30"))
for(i in seq_len(length(fit.counts.t$stages))){
  fit.counts.t$stages[i]  <- numbers[names(numbers) %in% fit.counts.t$stages[i]]
}
fit.counts.t$stages[1]

par(mar=c(3,3,3,3))
for(i in seq_len(length(genes.lineage))){
  png(filename = file.path(gene.expression.dir,paste0(genes.lineage[i],"_boxplot_FiPS.png")),width = 1000, height = 800)
  par(mar=c(3,3,3,3))
  
  boxplot(fit.counts.t[,genes.lineage[i]]~fit.counts.t$stages, frame=FALSE,cex.lab =2,cex.main=2, cex.axis=2, col = mycolors[-1], xlab = "Stages", ylab= "TPM", main = names(fit.counts.t)[i])
  boxplot(fit.counts.t[,genes.lineage[2]]~fit.counts.t$stages, frame=FALSE,cex.lab =2,cex.main=2, cex.axis=2, col = mycolors[-1], xlab = "Stages", ylab= "TPM", main = names(fit.counts.t)[i])
  
  # boxplot(x = fit.counts.t$stages, y = fit.counts.t[,genes.lineage[i]]~fit.counts.t$stages, frame=FALSE,cex.lab =2,cex.main=2, cex.axis=2, col = mycolors[-1], xlab = "Stages", ylab= "TPM", main = names(fit.counts.td)[i])
  
  # png(filename = file.path(gene.expression.dir,paste0(genes.lineage[i],"_FiPS.png")),width = 1000, height = 800)
  # library("ggpubr")
  # ggviolin(fit.counts.t, x = "stages", y = genes.lineage[1], fill = "stages", add = "boxplot", main = genes.lineage[1], cex = 4)
  dev.off()
}

# Barplot gene expression
png(filename = file.path(file.path(markers.dir, paste0("Barplotsclineage.png"))),
    width = 3000, height = 2500)
par(mar = c(20, 20, 20, 20)+ 0.1)
barCenters <-barplot(t(boxplot.data)[,1], beside=TRUE,  
        cex.names=0.5, las=2, ylim=c(0,20),  cex.axis = 0.5, cex.lab = 0.5)
arrows( barCenters, boxplot.data.g-boxplot.data.sd,
        barCenters, boxplot.data.g+boxplot.data.sd,angle=90, code= 3, lwd =2)
mtext(side=2, line=8, "log2(Expression count)", font=2, cex=3)
legend("topleft", fill = color.data, legend = names(color.data), horiz = T, cex = 5)
dev.off()


# Lines as gene expression

gene.color <- c("#388D36", "#388D36","#EF171A","#EF171A","#650A61","#EF171A","#650A61","#EE9D3A","#650A61","#EF171A")
genes <- c("NGFR","SOX10", "CDH19", "ITGA4","PLP1", "GAP43", "S100B", "EGR2", "PMP22","MPZ")
names(gene.color) <- genes
gene.color.red <- gene.color [which(gene.color =="#EF171A")]
gene.color.green <- gene.color[which(gene.color =="#388D36")]
gene.color.purple <- gene.color[which(gene.color =="#650A61")]
gene.color.orange <- gene.color[which(gene.color =="#EE9D3A")]

sample.data

matrix.list <- list()
ms <- data.frame()
for(i in 1:length(genes)){
  gene <- genes[i]
  
  matrix.list[[gene]]<- generationGeneMatrix(dds = dds.deg, 
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


boxplot(fit.counts.t[,genes.lineage[i]]~fit.counts.t$stages, frame=FALSE,cex.lab =2,cex.main=2, cex.axis=2, col = mycolors[-1], xlab = "Stages", ylab= "TPM", main = names(fit.counts.t)[i])
gene.color.green <- fit.counts.t[names(gene.color.green)]
gene.color.green$stages <- fit.counts.t$s
mv <- max(gene.color.green[,1])
gene.color.green[,2] <- gene.color.green[,2]/mv
boxplot(gene.color.green[,1]~gene.color.green$stage, boxfill = NA, border = NA) #invisible boxes - only axes and plot area
boxplot(gene.color.green[,1]~gene.color.green$stages, xaxt = "n", add = TRUE, boxfill="#388D36", 
        boxwex=0.25) #shift these left by -0.15
boxplot(gene.color.green[,2]~gene.color.green$stages, xaxt = "n", add = TRUE, boxfill="#CEFFCC", 
        boxwex=0.25 ) #shift to the right by +0.15

lines(x =c(1,2,3,4,5), y = p/mv, col = "#388D36")
text(x=6, y = p[5], labels ="NGFR", cex = 2, pos = 2, col = "black")

karyoploteR::lighter("#388D36")


png(filename = file.path(gene.expression.dir, "GreenColor_Markers_RTPCR.2.png"), width = 1000, height = 800)
par(oma=c( 0,0,0,0), mar =c(5,5,3,1)+0.1, mgp=c(3,1,0))
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




######################## GO Terms Representation  #############################
cont <- "Fb"
ont <- "BP"
for(cont in names(markers)){
  for(ont in goterms){
    
    # if(!file.exists(file.path(gsea.dir, ont, paste0(cont,"_markers_", ont, "_table.csv")))) next
    if(!file.exists(file.path(gsea.dir, ont, paste0("Enrichment", ont, "markers_up_", cont, ".csv"))))next
    if(is.null(count.fields(file.path(gsea.dir, ont, paste0("Enrichment", ont, "markers_up_", cont, ".csv"))))) next
    # go.ont <- read.table(file.path(gsea.dir, ont, paste0(cont,"_markers_", ont, "_table.csv")), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    go.ont <- read.table(file.path(gsea.dir, ont, paste0("Enrichment", ont, "markers_up_", cont, ".csv")), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    
    # colnames(go.ont) <- c("ID","Description",	"GeneRatio",	"BgRatio",	"pvalue",	"p.adjust",	"qvalue",	"geneID", "Count")
    if(nrow(go.ont)== 0)next
    mat <- t(matrix(head(go.ont$Count,10), dimnames = list(head(go.ont$Description,10), "Count")))
    wr.lap <- wrap.labels(colnames(mat), 60)
    colnames(mat)<- wr.lap
    
    png(filename = file.path(gsea.dir, ont, paste0("Barplot_enrich_", ont,"_", cont, "_markers.png")),
        width = 4000, height = 2000)
    
    par(mar=c(5,20,4,2)+0.1, oma=c(10,100,10,10))
    barplot(height = mat,
            horiz = TRUE, 
            las=2, 
            xlim =c(0,max(mat[1,])+5),
            main = paste0(ont," ", cont),
            # main = paste(enrichterm, " ", group, ".SC ", clust),
            col= "lightblue", 
            cex.names = 5,
            cex.axis = 6, 
            cex.main= 6)
    abline(v=seq(0,max(as.numeric(mat[1,]))+5, 5))
    
    dev.off()
  }
}


############################  Gene expression plots ############################
dds.deg <- DESeq(filtered.dds)
dds.rlog <- rlog(dds.deg)
genes.lineage <- c("NGFR","SOX10", "EGR2","CDH19", "ITGA4","MPZ","GAP43","PLP1", "PMP22", "S100B")
for (g in genes.lineage){
  count.matrix <- generationGeneMatrix(dds.deg, gene = g, sname_variable = "Sample.Name", contrast_group = samples.group)
  count.matrix <- count.matrix [c(4,3,1,2)]
  names(count.matrix) <- c("0","7","14","30")
  png(filename = file.path(gene.expression.dir, paste0(g, "_", model, ".png")), width = 1000, height = 1000)
  
  plotLimsAndLabels(gene_matrix = count.matrix, gene = g,xlab = "Stages", xaxt = "n", bty= "n", cex.lab = 2.5, cex.axis = 2.5, cex.main = 2.5, ylab = "")
  axis(side = 1, lwd = 2, at= c(0,7,14,30), labels = c("NC","day7","day14","day30"), cex.axis= 2.5)
  
  # plotPointsAndMeans(gene_matrix = count.matrix[5:16,], plot_points = TRUE, plot_legend = FALSE, line_col = "black", lwd = 2)
  # plotPointsAndMeans(gene_matrix = count.matrix[17:20,], plot_points = TRUE, plot_legend = FALSE, line_col = "red", lwd = 2)
  # plotPointsAndMeans(gene_matrix = count.matrix[21:24,], plot_points = TRUE, plot_legend = FALSE, line_col = "green", lwd = 2)
  # plotPointsAndMeans(gene_matrix = count.matrix[1:4], plot_points = TRUE, plot_legend = FALSE, line_col = "purple", lwd = 2)
  # legend("topright",legend = c("FiPS","3MM","5MM","D12MM"),
  #        fill = c("black","red", "green", "purple") ,col =c("black","red", "green", "purple"))
  # 
  plotPointsAndMeans(gene_matrix = count.matrix[14:25,], plot_points = TRUE, plot_legend = FALSE, line_col = "black", lwd = 1)
  plotPointsAndMeans(gene_matrix = count.matrix[26:29,], plot_points = TRUE, plot_legend = FALSE, line_col = "red", lwd = 1)
  plotPointsAndMeans(gene_matrix = count.matrix[30:33,], plot_points = TRUE, plot_legend = FALSE, line_col = "green", lwd = 1)
  plotPointsAndMeans(gene_matrix = count.matrix[10:13,], plot_points = TRUE, plot_legend = FALSE, line_col = "purple", lwd = 1)
  plotPointsAndMeans(gene_matrix = count.matrix[c(1,2,7),], plot_points = TRUE, plot_legend = FALSE, points_col = "pink", points_size = 1)
  plotPointsAndMeans(gene_matrix = count.matrix[3,], plot_points = TRUE, plot_legend = FALSE, points_col = "blue", points_size = 1)
  plotPointsAndMeans(gene_matrix = count.matrix[4,], plot_points = TRUE, plot_legend = FALSE, points_col = "lightblue", points_size = 1)
  plotPointsAndMeans(gene_matrix = count.matrix[5,], plot_points = TRUE, plot_legend = FALSE, points_col = "grey", points_size = 1)
  plotPointsAndMeans(gene_matrix = count.matrix[6,], plot_points = TRUE, plot_legend = FALSE, points_col = "orange", points_size = 1)
  plotPointsAndMeans(gene_matrix = count.matrix[8,], plot_points = TRUE, plot_legend = FALSE, points_col = "brown", points_size = 1)
  plotPointsAndMeans(gene_matrix = count.matrix[9,], plot_points = TRUE, plot_legend = FALSE, points_col = "yellow", points_size = 1)
  legend("topright",legend = c("FiPS","3MM","5MM","D12MM","FiPS3D_5MMFb", "3MM3D","3MM3D_Het","5MM3D","5MM3D_Het","D123D","D123D_het"),
         fill = c("black","red", "green", "purple", "pink","blue", "lightblue","grey","orange","brown","yellow") ,col =c("black","red", "green", "purple", "pink","blue", "lightblue","grey","orange","brown","yelow") )
  dev.off()
  
}













#Barplot gene expression
genes.lineage <- c("NGFR","SOX10", "ITGA4","CDH19", "PMP22","MPZ","PLP1","EGR2","S100B","GAP43")
samples.group <- "Diff.Day"
g < "CDKN2A"
for (g in genes.lineage){
count.matrix <- generationGeneMatrix(dds.deg, gene = g, sname_variable = "Sample.Name", contrast_group = samples.group)
# count.matrix <- count.matrix [c(4,3,1,2)]
# names(count.matrix) <- c("0","7","14","30")
# 
# fips <- count.matrix[14:25,]
# fips <- colMeans(fips, na.rm = TRUE)
# MM3 <- count.matrix[26:29,]
# 
# MM5 <-count.matrix[30:33,]
# D12 <-count.matrix[10:13,]
# fips3d <- count.matrix[c(1,2,7),]
# fips3d <- colMeans(fips3d, na.rm = TRUE)
# MM3_3D <- count.matrix[3,]
# MM3_3D_het <- count.matrix[4,]
# MM5_3D <- count.matrix[5,]
# MM5_3D_het <- count.matrix[6,]
# D12_3D <- count.matrix[8,]
# D12_3D_het <- count.matrix[9,]
count(dds.deg,normalize=TRUE)
 
count.matrix <- count.matrix [c(5,4,2, 3, 1)]
names(count.matrix) <- c("0","7","14","30", "40")

fips <- count.matrix[19:30,]
fips <- colMeans(fips, na.rm = TRUE)
MM3 <- count.matrix[31:34,]
MM5 <-count.matrix[35:38,]
D12 <-count.matrix[15:18,]
fips3d <- count.matrix[c(6,7,12),]
fips3d <- colMeans(fips3d, na.rm = TRUE)
MM3_3D <- count.matrix[8,]
MM3_3D_het <- count.matrix[9,]
MM5_3D <- count.matrix[10,]
MM5_3D_het <- count.matrix[11,]
D12_3D <- count.matrix[13,]
D12_3D_het <- count.matrix[14,]
fibros <- count.matrix[c(1:4),]
SC <- count.matrix[c(5,39:40),]


nrow(count.matrix)



png(filename = file.path(gene.expression.dir, paste0(g, "_", model, "AllSamples_bars.png")), width = 1000, height = 1000)

# plotbars <- as.matrix(rbind(fips, MM3,MM5,D12, fips3d,MM3_3D, MM3_3D_het, MM5_3D, MM5_3D_het, D12_3D, D12_3D_het))
plotbars <- as.matrix(rbind(fips, MM3,MM5,D12, fips3d,MM3_3D, MM3_3D_het, MM5_3D, MM5_3D_het, D12_3D, D12_3D_het,SC,fibros))

class(plotbars)
# day14 <- plotbars[!is.na(plotbars[,3]),3]
# names(day14) <- c("FiPS","3MM","5MM","D12MM","FiPS3D_5MMFb", "3MM3D","3MM3D_Het","5MM3D","5MM3D_Het","D123D","D123D_het")
# barplot(day14, beside = TRUE, las=2, col = c("black","red", "green", "purple", "pink","blue", "lightblue","grey","orange","brown","yellow"), main =g)
# legend("topleft",legend = c("FiPS","3MM","5MM","D12MM","FiPS3D_5MMFb", "3MM3D","3MM3D_Het","5MM3D","5MM3D_Het","D123D","D123D_het"),
#        fill = c("black","red", "green", "purple", "pink","blue", "lightblue","grey","orange","brown","yellow"))

day14.cells <- plotbars[,c(3,5)]
rownames(day14.cells)[1] <- "FiPS"
rownames(day14.cells)[14] <- "FiPS3D_5MMFb"
day14.cells <- day14.cells[c(1,4,8,10,14:20,21:27),]




day14.cells <- t(day14.cells)
day14.cells[1,c(12:18)] <-day14.cells[2,c(12:18)]

par(mar=c(10,10,10,20))
barplot(day14.cells[1,], beside = TRUE, las=2, col = c("black","red", "green",  "purple", "darkgreen","blue", "cyan","orange","brown","yellow", "aquamarine","lightblue", "lightblue", "lightblue","grey","grey","grey","grey"), main =g)
legend("topright", inset=c(-0.25,0),xpd = TRUE,legend = colnames(day14.cells),ncol=1,
       fill = c("black","red", "green",  "purple", "darkgreen","blue", "cyan","orange","brown","yellow", "aquamarine","lightblue", "lightblue", "lightblue","grey","grey","grey","grey"))

dev.off()

}



####################################################################
#                                 FB and SC                        #
####################################################################
gene.markers <- markers$Fb
counts.deg  <-counts(dds.deg, normalize = TRUE)
counts.deg <- counts.deg[rownames(counts.deg) %in% gene.markers,]

counts.deg.fb <- counts.deg[,grepl(pattern ="_F", x=colnames(counts.deg))]
sd.deg.fb <- apply(counts.deg.fb,1,sd)
sd.deg.fb <- sd.deg.fb[gene.markers]
counts.deg.fb <- rowMeans(counts.deg.fb)
counts.deg.fb <- counts.deg.fb[gene.markers]

counts.deg.sc <- counts.deg[,grepl(pattern ="SC", x=colnames(counts.deg))]
sd.deg.sc <- apply(counts.deg.sc,1,sd)
sd.deg.sc <- sd.deg.sc[gene.markers]
counts.deg.sc <- rowMeans(counts.deg.sc)
counts.deg.sc <- counts.deg.sc[gene.markers]


pm.terms <- c("GO:0005886", "GO:0009897")
term <- "CC"
annotated.pm <- read.table(file = file.path("Results/FbvsSC/Annotation/Fb_up_CC_markers_annotation.csv"), header = TRUE, sep = "\t",  stringsAsFactors = FALSE)
annotated.pm <- annotated.pm[annotated.pm$CC.go_id %in% pm.terms,]

res.markers <- deseq.results$Fb$SC[rownames(deseq.results$Fb$SC) %in%gene.markers,]
# res.markers <- dds.deg[rownames(dds.deg)%in% gene.markers,]
res.markers.lfc <- res.markers$log2FoldChange
names(res.markers.lfc) <- rownames(res.markers)
res.markers.lfc <- res.markers.lfc[gene.markers]



g <- annotated.pm$CC.hgnc_symbol
names(g) <- annotated.pm$CC.hgnc_symbol
g <- g[gene.markers]
name_1006 <- annotated.pm$CC.name_1006
names(name_1006) <- annotated.pm$CC.hgnc_symbol
name_1006 <- name_1006[gene.markers]


table.PlasmaMembrane <- data.frame(genes = gene.markers, term = name_1006, 
                                   lfc = abs(res.markers.lfc), counts.fb=counts.deg.fb,
                                   sd.fb= sd.deg.fb, counts.sc=counts.deg.sc, sd.sc = sd.deg.sc)
table.PlasmaMembrane <- table.PlasmaMembrane[!is.na(table.PlasmaMembrane$term),]
table.PlasmaMembrane <- table.PlasmaMembrane[order(table.PlasmaMembrane$lfc,decreasing = T),]

write.table(x = table.PlasmaMembrane, file = file.path(annot.dir, "Annotation_PlasmaMembrane_Fb_Markers.csv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)




