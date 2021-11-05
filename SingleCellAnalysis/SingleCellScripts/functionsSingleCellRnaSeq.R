#############################  FUNCTIONS  #########################

#Scoring cell selection functions
AddMarkersModuleScore <- function(sce,
                                  gene.set.list,
                                  analysis.name,
                                  specific.clusters,
                                  coldata.clusters,
                                  ctrl = 50,...){
  
  if(class(sce) != "SingleCellExperiment"){
    stop("sce parameter must be a SingleCellExperiment")
  }
  
  if(!is.list(gene.set.list)){
    stop("gene.set.list must be a list")
  }
  
  if(!is.character(analysis.name)){
    stop("analysis.name must be a character")
  }
  
  if(!is.null(specific.clusters)){
    if(!is.numeric(specific.clusters)){
      stop("num.clusters must be a numerical vector")
    }
  }
  
  #coldata.cluster test
  if(!is.character(coldata.clusters)){
    stop("coldata.clusters must be a character name of the sce colData colums")
  }
  
  if(!coldata.clusters %in% names(colData(sce))){
    stop("coldata.clusters must be a character name of the sce colData colums")
  }
  
  
  #If specific.clusters is not null that means the user wants to select specific clusters
  #First we stract all clusters
  num.all.clusters <- sort(unique(as.numeric(colData(sce)[,coldata.clusters])))
  
  if(is.null(specific.clusters)){
    message("seleccting all clusters")
    specific.clusters <- num.all.clusters
    
  }else{
    message ("selecting specific clusters")
    specific.clusters <- specific.clusters
  }
  # sce <- sce.hvg
  
  # Then we select the cells of interst before scoring them
  sce.sel <- colnames(sce)[colData(sce)[[coldata.clusters]] %in% specific.clusters]
  sce.sel <- sce [,sce.sel]
  colData(sce.sel)[[coldata.clusters]] <- factor(as.numeric(colData(sce.sel)[[coldata.clusters]]), levels=sort(unique(as.numeric(colData(sce.sel)[[coldata.clusters]]))))
  colLabels(sce.sel)[[coldata.clusters]] <-colData(sce.sel)[[coldata.clusters]]
  
  
  #We process sce.hvg data to convert it to saurat objet
  marrow <- as.Seurat(sce.sel)
  
  #Prior filtering of the feature list before scoring.
  to.del <- NULL
  for(i in seq_len(length(gene.set.list))){
    x <- gene.set.list[[i]]
    z <- rownames(marrow)[rownames(marrow) %in% x]
    if(length(z)==0){
      msg(names(gene.set.list)[i], " has no genes in the seurat object to get the cell scoring")
      to.del <- i
    }else{
      gene.set.list[[i]] <-z
    }
    
  }
  
  if(!is.null(to.del)){
    gene.set.list <- gene.set.list[-to.del]
    
  }
  
  
  marrow <-AddModuleScore(object = marrow,
                          features = gene.set.list,
                          name = analysis.name,
                          ctrl = ctrl,...)
  
  
  return(marrow)
}

#Example


#L2.norm function
l2.norm <- function(x){
  norm(x, type = "2")
}




# GetCell groups By Gene Exp
getCellgroupsByGenesExp <- function(sce,
                                    gene.set.list,
                                    analysis.name,
                                    specific.clusters = NULL,
                                    coldata.clusters ="",
                                    num.clusters = NULL,
                                    ctrl = 50,...){
  
  
  #First we obtain the module score given a group of genes
  module.score.markers <- AddMarkersModuleScore(sce = sce,
                                                gene.set.list = gene.set.list,
                                                analysis.name = analysis.name,
                                                specific.clusters= specific.clusters,
                                                coldata.clusters= coldata.clusters,
                                                ctrl = 50,...)
  
  #names of seurat objec variables 
  analysis.names <- paste0(analysis.name, seq_len(length(gene.set.list))) #names of seurat objec variables
  
  #Normalizing the scoring before specifying the scoring group
  # L2-norming of ech group of module scores before making any comparison.
  normalize.score <- list()
  i=1
  for(i in seq_len(length(analysis.names))){
    an <- analysis.names[i]
    
    #Normalization
    normalizer <- l2.norm(module.score.markers[[an]])
    normalize.score[[an]] <-module.score.markers[[an]]/normalizer
    
  }
  
  # for(i in seq_len(length(gene.set.list))){
  #   nm <- names(gene.set.list)[i]
  #   gene.set <- gene.set.list[[nm]]
  # 
  #   #obtaining pseudo bulk rna-seq data by gene set
  #   p.bulk <- rowSums(assay(sce.sel))[gene.set]
  #   #the mean of pseudo-bulk counts
  #   mean.p.bulk <- mean(p.bulk)
  # 
  #   #normalizer calculation
  #   normalizer  <- p.bulk/mean.p.bulk
  #   an <- analysis.names[i]
  #   #Normalization
  #   normalize.score[[an]] <-module.score.markers[[an]]/normalizer
  # }
  
  #Selecting those cells with a scoring >0
  for(i in seq_len(length(analysis.names))){
    an <- analysis.names[i]
    analysis.data <-normalize.score[[an]]
    selected.cells <- rownames(analysis.data)[analysis.data > 0]# we want to select those groups with positive expresion (overexpresed genes)
    
    # We select those cells present in the specific clusters
    normalize.score[[an]]<- normalize.score[[an]][selected.cells,]
    names(normalize.score[[an]]) <- selected.cells
  }
  
  #common cells in gene.set.list terms
  
  for (i in seq_len(length(analysis.names))){
    a1 <- analysis.names[i]
    
    # scoring.cells.1 <- rownames(normalize.score[[a1]])
    
    anms <- analysis.names[!analysis.names %in% a1]
    
    common.a1.del<- c()
    
    for(j in seq_len(length(anms))){
      a2 <- anms [j]
      common.cells <- Reduce(intersect, list(names(normalize.score[[a1]]),names(normalize.score[[a2]])))
      common.a1 <- normalize.score[[a1]][common.cells] 
      common.a2 <- normalize.score[[a2]][common.cells] 
      common.a1.del<- c(common.a1.del,names(common.a1)[common.a1 < common.a2])
      
    }
    
    normalize.score[[a1]] <-normalize.score[[a1]][-which(names(normalize.score[[a1]])%in%common.a1.del )]
    
  }
  return(normalize.score)
}

# sce <- sce.hvg

plotSC <- function(sce, dimred="TSNE", col="grey", pch=16, cex=1, main="", main.cex=4, cex.axis =1,  xlab="TSNE1", ylab="TSNE2", lwd.axis = 2, cex.lab =1) {
  plot(SingleCellExperiment::reducedDim(sce, "TSNE"), 
       xlab=xlab, 
       ylab=ylab, 
       main=main,
       cex=cex,
       pch=pch,
       col=col,
       bty="L", 
       yaxt = "n", 
       xaxt = "n",
       cex.lab = cex.lab,axes=FALSE)
  
  #TODO: Muntar els axes a part, ontrolar nums sÃ�/no, tamanys lletres
  ticks.x <- seq(from =round(min(reducedDim(sce.hvg, "TSNE")[,1])-10), to= round(max(reducedDim(sce.hvg, "TSNE")[,1])+10),10)
  ticks.y <- seq(from =round(min(reducedDim(sce.hvg, "TSNE")[,2])-10), to= round(max(reducedDim(sce.hvg, "TSNE")[,2])+10),10)
  
  axis(1, 
       at = ticks.x,
       cex.axis =cex.axis,
       lwd = lwd.axis,
       gap.axis = 0)
  axis(2,
       at = ticks.y,
       cex.axis =cex.axis, 
       lwd = lwd.axis,
       gap.axis = 0)
  
  #TODO: Afegir main
}

# plotSC(sce, cex=0.4)
# 
# g1 <- counts(sce)["IL6",]
# cr <- colorRamp(c("white", "red"))
# 
# plotSC(sce, col=rgb(cr(g1/max(g1))/255), cex=0.4)
# plotSC(sce, col=cr(g1), cex=0.4)
# 
# gg <-df.new.markers$genes[df.new.markers$stages == "NC"]
# g2 <- counts(sce)["CDH19",]
# cr <- colorRamp(c("white", "blue"))
# cr <- colorRamp(c("grey", "red"),)
# 
# g2 <- d30cells
# cr <-cr[!is.na(cr(sce$label)[,1]),]
# rgb(cr(sce$label)/255)
# 
# plotSC(sce, col=rgb(cr(g2/max(g2))/255), cex=0.4)
# plotSC(sce, col=rgb(cr(g2)/255), cex=0.4)
# 
# rgb(cr(sce$label)/255)
# 
# 
# col1 <- cr(g1/max(g1))/255
# col2 <- cr(g2/max(g2))/255
# col1 <- cr(sce$label)/255
# plotSC(sce, col=rgb(col1 <- cr(sce$label)/255))
# 
# mm <-(1- (1 - col1) + (1- col2))/2
# merge.col <- 1-((1-col1)+(1-col2))
# 
# mm <- (col1+col2)/2
# plotSC(sce, col=rgb(mm))
# sce$label
# 


plotMarkers <- function(sce, genes, dimred="TSNE", ncol=3, add.legend=FALSE, title.text.size=15, x.axis.size = 1, y.axis.size= 1, axis.title.size =15, size.x.axis.line = 1, size.y.axis.line = 1, ticks.size = 1) {
  plotlist <- list()
  for(i in genes) {
    plotlist[[i]] <- plotReducedDim(sce,dimred = dimred, colour_by = i,by_exprs_values = "logcounts") +
      scale_fill_gradientn(colours = colorRampPalette(c("grey90","orange3","firebrick","firebrick","red","red" ))(10)) +
      ggtitle(label = i)+ theme(plot.title = element_text(size=title.text.size), 
                                axis.text.x = element_text(face = "bold",size = x.axis.size),#x.axis text size
                                axis.text.y = element_text(face = "bold",size = y.axis.size),#y.axis text size
                                axis.line.y.left = element_line(size= size.y.axis.line, linetype ="solid") , 
                                axis.line.x.bottom = element_line(size= size.y.axis.line, linetype ="solid"),
                                axis.ticks = element_line(size = ticks.size),
                                axis.title = element_text(face = "bold", size= axis.title.size)) 
    if(!add.legend) {
      plotlist[[i]] <- plotlist[[i]] + theme(legend.position = "none") 
    }
  }
  
  return(plot_grid(ncol=ncol, plotlist = plotlist))
}

