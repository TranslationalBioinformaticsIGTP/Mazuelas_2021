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
source(file = "./rna_seq_Functions.R")

############################  Parameters  ##########################
# # FiPS parameters
# params <- yaml.load_file("./Parameters/FiPS_pipeline_parameters.yaml")
## NF1 parameters
# params <- yaml.load_file("./Parameters/NF1_2D_pipeline_parameters.yaml")
# # NF1_A 5MM parameters
# sn <- "NF1_A"
# # # NF1_B 3MM parameters
# # sn <- "NF1_B"
# NF1_C D12 parameters
# sn <- "NF1_C"
# # FiPS NF1 2D parameters
# params <- yaml.load_file("./Parameters/FiPS_NF1_2D_pipeline_parameters.yaml")
# # FiPS + SC parameters
# params <- yaml.load_file("./Parameters/FiPS_SC_pipeline_parameters.yaml")
# # Sphers hom vs het
# params <- yaml.load_file("./Parameters/Spherers_pipeline_paramaters.yaml")
# # FiPS + Spheres parameters
# params <- yaml.load_file("./Parameters/FiPS_2Dvs3D_pipeline_parameters.yaml")
#All samples parameters
params <- yaml.load_file("./Parameters/AllSamples_pipeline_parameters.yaml")
# #FB vs SC parameters
# params <- yaml.load_file("./Parameters/SCvsFB_pnf_pipeline_parameters.yaml")


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
  mks<- read.table(file = file.path('./Results/2D/markers', paste0("up_", gr,".txt")), header = FALSE,stringsAsFactors = FALSE)[,1]
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



####################### Code not used in the paper ###################
# ####################################################################
# #                                 FB and SC                        #
# ####################################################################
# gene.markers <- markers$Fb
# counts.deg  <-counts(dds.deg, normalize = TRUE)
# counts.deg <- counts.deg[rownames(counts.deg) %in% gene.markers,]
# 
# counts.deg.fb <- counts.deg[,grepl(pattern ="_F", x=colnames(counts.deg))]
# sd.deg.fb <- apply(counts.deg.fb,1,sd)
# sd.deg.fb <- sd.deg.fb[gene.markers]
# counts.deg.fb <- rowMeans(counts.deg.fb)
# counts.deg.fb <- counts.deg.fb[gene.markers]
# 
# counts.deg.sc <- counts.deg[,grepl(pattern ="SC", x=colnames(counts.deg))]
# sd.deg.sc <- apply(counts.deg.sc,1,sd)
# sd.deg.sc <- sd.deg.sc[gene.markers]
# counts.deg.sc <- rowMeans(counts.deg.sc)
# counts.deg.sc <- counts.deg.sc[gene.markers]
# 
# 
# pm.terms <- c("GO:0005886", "GO:0009897")
# term <- "CC"
# annotated.pm <- read.table(file = file.path("Results/FbvsSC/Annotation/Fb_up_CC_markers_annotation.csv"), header = TRUE, sep = "\t",  stringsAsFactors = FALSE)
# annotated.pm <- annotated.pm[annotated.pm$CC.go_id %in% pm.terms,]
# 
# res.markers <- deseq.results$Fb$SC[rownames(deseq.results$Fb$SC) %in%gene.markers,]
# # res.markers <- dds.deg[rownames(dds.deg)%in% gene.markers,]
# res.markers.lfc <- res.markers$log2FoldChange
# names(res.markers.lfc) <- rownames(res.markers)
# res.markers.lfc <- res.markers.lfc[gene.markers]
# 
# 
# 
# g <- annotated.pm$CC.hgnc_symbol
# names(g) <- annotated.pm$CC.hgnc_symbol
# g <- g[gene.markers]
# name_1006 <- annotated.pm$CC.name_1006
# names(name_1006) <- annotated.pm$CC.hgnc_symbol
# name_1006 <- name_1006[gene.markers]
# 
# 
# table.PlasmaMembrane <- data.frame(genes = gene.markers, term = name_1006, 
#                                    lfc = abs(res.markers.lfc), counts.fb=counts.deg.fb,
#                                    sd.fb= sd.deg.fb, counts.sc=counts.deg.sc, sd.sc = sd.deg.sc)
# table.PlasmaMembrane <- table.PlasmaMembrane[!is.na(table.PlasmaMembrane$term),]
# table.PlasmaMembrane <- table.PlasmaMembrane[order(table.PlasmaMembrane$lfc,decreasing = T),]
# 
# write.table(x = table.PlasmaMembrane, file = file.path(annot.dir, "Annotation_PlasmaMembrane_Fb_Markers.csv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)




