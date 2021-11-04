# funtion to obtain the normalization of the expresion for expression plots
getNormalizedFiltExp <- function(filtered.dds,
                                 genes,
                                 sample.data,
                                 norm.data.column,
                                 sample.groups){
  
  # Counts normalization
  filt.counts <- counts(object = filtered.dds, normalize = TRUE)
  
  #selecting genes to normalize 
  filt.counts <- filt.counts[rownames(filt.counts) %in% genes,]

  # Preparing the data for normalization
  fct <- t(filt.counts)

  # Normalization of the values
  for(i in seq_len(ncol(fct))){
    gene <- colnames(fct)[i]
    mv <- max(fct[,gene])
    fct[,gene] <- fct[,gene]/mv
  }
  
  fct <- data.frame(fct)
  fct<- cbind(fct, data.group = sample.data[,norm.data.column])

  fct$sample.group <- ""
  
   for(i in seq_len(length(sample.groups))){
    fct$sample.group[grepl(x=rownames(fct),pattern = sample.groups[i])] <- sample.groups[i]
   }
  
  return(fct)
}

# # Test getNormalizedFiltExp funtion
# norm.data.column <- "Diff.Day"
# filtered.dds <- filtered.dds
# genes <- genes
# sample.data <- sample.data
# # sample.groups <- c("WT","NF1_A","NF1_B", "NF1_C")
# sample.groups <- c("WT")
# # sample.groups <- c("WT","NF1_A","NF1_B", "NF1_C","WT_Het","NF1_A_Het","NF1_B_Het", "NF1_C_Het","NF1_A_Hom","NF1_B_Hom", "NF1_C_Hom")
# 
# norm.data <- getNormalizedFiltExp(filtered.dds = filtered.dds,
#                                   genes = genes,
#                                   sample.data = sample.data,
#                                   norm.data.column = "Diff.Day",
#                                   sample.groups = sample.groups )

i=1
#Funtion to get Mean and Desviation of the samples 
getMeanNomrFiltExp <-function(norm.data){
  
  # Getting norm.data by mean.sample.group
  norm.data.list <- split( norm.data , f = norm.data$sample.group)

  m.list <- list()
  
  #Computing mean of norm.data.list by gene
  for (i in seq_len(length(norm.data.list))){
    gr <- names(norm.data.list)[i]
    f <-norm.data.list[[gr]]
    f <- f[,-c(ncol(f), ncol(f)-1)]
      
    m <- matrix(ncol = ncol(f), nrow=length(unique(norm.data.list[[gr]]$data.group)))
    j=6
    for(j in seq_len(ncol(f))){
      gene <- colnames(f)[j]
      
      for(s in seq_len(length(unique(norm.data.list[[gr]]$data.group)))){
        stg <- unique(norm.data.list[[gr]]$data.group)[s]
        smp.stg <- norm.data.list[[gr]][norm.data.list[[gr]]$data.group ==stg,]
        if(nrow(smp.stg) ==1){
          m[s,j] <- as.numeric(as.character(f[norm.data.list[[gr]]$data.group %in% stg, colnames(f) == gene]))
        }else{
          m[s,j]<- mean(as.numeric(as.character(f[norm.data.list[[gr]]$data.group %in% stg, colnames(f) == gene])))
        }
      }
      
    }
    rownames(m) <- unique(norm.data.list[[gr]]$data.group)
    colnames(m) <- colnames(f)
    m.list[[gr]] <- data.frame( m, data.group = unique(norm.data.list[[gr]]$data.group))
    
  }
  
   return(m.list)
}

# # Test getMeanNomrFiltExp
# m.list <- getMeanNomrFiltExp(norm.data = norm.data,genes = genes)

getSdNomrFiltExp <-function(norm.data){
  # Getting norm.data by mean.sample.group
  norm.data.list <- split( norm.data , f = norm.data$sample.group)
  desv.list <- list()
  
  #Computing sd of norm.data.list by gene
  
  for (i in seq_len(length(norm.data.list))){
    gr <- names(norm.data.list)[i]
    f <-norm.data.list[[gr]]
    f <- f[,-c(ncol(f), ncol(f)-1)]
    
    desv <- matrix(ncol = ncol(f), nrow=length(unique(norm.data.list[[gr]]$data.group)))
    
    for(j in seq_len(ncol(f))){
      gene <- colnames(f)[j]
      
      for(s in seq_len(length(unique(norm.data.list[[gr]]$data.group)))){
        stg <- unique(norm.data.list[[gr]]$data.group)[s]
        smp.stg <- norm.data.list[[gr]][norm.data.list[[gr]]$data.group ==stg,]
        if(nrow(smp.stg) ==1){
          desv[s,j] <- 0
        }else{
          desv[s,j] <- sd(as.numeric(as.character(f[norm.data.list[[gr]]$data.group %in% stg, colnames(f) == gene])))
        }
      }
      
    }
    rownames(desv) <- unique(norm.data.list[[gr]]$data.group)
    colnames(desv) <- colnames(f)
    desv.list[[gr]] <- data.frame( desv, data.group = unique(norm.data.list[[gr]]$data.group))

  }
  
  return(desv.list)
}

# #test getSdNomrFiltExp
# desv.list <- getSdNomrFiltExp(norm.data = norm.data, genes)

# 
# 
# filt.counts <- counts(filtered.dds, normalize = TRUE)
# 
# 
# #selecting the PCR marker genes 
# filt.counts <- filt.counts[rownames(filt.counts) %in% genes,]
# # colnames(filt.counts) <- c("PSC","NC","NC", "NC","day7","day7","day7","day14","day14","day14","day30","day30","day30")
# 
# # Preparing the data for boxplot representation
# fct <- t(filt.counts)
# # rownames(fct) <- genes
# # Normalization of the values
# for(i in seq_len(ncol(fct))){
#   gene <- colnames(fct)[i]
#   mv <- max(fct[,gene])
#   fct[,gene] <- fct[,gene]/mv
# }
# fct <- data.frame(fct)
# fct<- cbind(fct, stage = sample.data$Diff.Day)
# # fct<- cbind(fct, stage = sample.data$Diff.Day)
# 
# fct$group <- ""
# group <- c("WT","NF1_A","NF1_B", "NF1_C")
# # group <- c("WT")
# # group <- c("WT","NF1_A","NF1_B", "NF1_C","WT_Het","NF1_A_Het","NF1_B_Het", "NF1_C_Het","NF1_A_Hom","NF1_B_Hom", "NF1_C_Hom")
# for(i in seq_len(length(group))){
#   fct$group[grepl(x=rownames(fct),pattern = group[i])] <- group[i]
# }
# 
# #Computing mean and deviation of the group
# 
# fct.list <- split( fct , f = fct$group )
# m.list <- list()
# desv.list <- list()
# for (g in seq_len(length(group))){
#   gr <- group[g]
#   f <-fct.list[[gr]]
#   f <- f[,-c(ncol(f), ncol(f)-1)]
#   m <- matrix(ncol = length(genes), nrow=length(unique(fct.list[[gr]]$stage)))
#   desv <- matrix(ncol = length(genes), nrow=length(unique(fct.list[[gr]]$stage)))
#   
#   
#   for(i in seq_len(length(genes))){
#     gene <- genes[i]
#     for(s in seq_len(length(unique(fct.list[[gr]]$stage)))){
#       stg <- unique(fct.list[[gr]]$stage)[s]
#       if(length(f[f$stage %in% stg, colnames(f) == gene]) ==1){
#         m[s,i] <- as.numeric(as.character(f[fct.list[[gr]]$stage %in% stg, colnames(f) == gene]))
#         desv[s,i] <- 0
#       }else{
#         m[s,i]<- mean(as.numeric(as.character(f[fct.list[[gr]]$stage %in% stg, colnames(f) == gene])))
#         desv[s,i] <- sd(as.numeric(as.character(f[fct.list[[gr]]$stage %in% stg, colnames(f) == gene])))
#       }
#     }
#     
#   }
#   rownames(m) <- unique(fct.list[[gr]]$stage)
#   colnames(m) <- genes
#   m.list[[gr]] <- data.frame( m, stage = unique(fct.list[[gr]]$stage))
#   rownames(desv) <- unique(fct.list[[gr]]$stage)
#   colnames(desv) <- genes
#   desv.list[[gr]] <- data.frame( desv, stage = unique(fct.list[[gr]]$stage))
# }
