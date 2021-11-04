############################ Functions RNA-Seq #########################
source(file ="/imppc/labs/eslab/mmagallon/Projects/RNA-Seq-ips-differenciation/utils.R")
source(file = "/imppc/labs/eslab/mmagallon/Projects/RNA-Seq-ips-differenciation/gsea_functions.R")
source(file ="/imppc/labs/eslab/mmagallon/Projects/RNA-Seq-ips-differenciation/annotationBioMart.R")
source (file = "/imppc/labs/eslab/mmagallon/Projects/RNA-Seq-ips-differenciation/plotfunctions.R")
source(file = "/imppc/labs/eslab/mmagallon/Projects/RNA-Seq-ips-differenciation/deg_result_functions.R")
source(file ="/imppc/labs/eslab/mmagallon/Projects/RNA-Seq-ips-differenciation/annotationBioMart.R")

### salmonAlignment
salmonAlignment <- function(sample.name, salmonDir,
                            file1.suffix, file2.suffix,
                            fastqdir,
                            transcript.index,
                            output.suffix,
                            output.quants, 
                            libtype = "IU",
                            threads = 4,
                            verbose = TRUE){
  now.msg("Salmon starting...")
  file.name <- file.path(output.quants, paste0(sample.name, output.suffix))
  
  if (!file.exists(file.name)){
    message("missing file", file.name)
    full.command <- paste0(salmonDir, " quant -i",
                           transcript.index," -l ",
                           libtype, " -1 ", 
                           fastqdir,
                           sample.name,
                           file1.suffix, " -2 ", 
                           fastqdir, sample.name, 
                           file2.suffix, " --validateMappings -p ",
                           as.character(threads),
                           " -o ", 
                           output.quants,
                           sample.name, output.suffix )
    system(full.command, wait = TRUE)
  }
}

####importQuantData
importQuantsData <- function(quant.files, orgdb, orgdb.keytype, org.columns, tximpot.type= "salmon", verbose = TRUE, ...){
  now.msg("Tximport starting...", verbose = verbose)
  #translating transctipts to gene symbol
  k <- keys(orgdb, keytype = org.keytype)
  tx2gene <- AnnotationDbi::select(x = orgdb,
                                   keys = k,
                                   columns = org.columns,
                                   keytype = org.keytype)
  
  
  #importing qunats data
  txi.salmon <- tximport(files = salmonquants.fl,
                         type = "salmon",
                         tx2gene = tx2gene, ...)
  now.msg("  Tximport done", verbose = verbose)
  return(txi.salmon)
}

###selectDataFromTximport  

selectDataFromTximport <- function(tximport, sample_names){
  tximport$abundance <- data.matrix(data.frame(tximport$abundance)[colnames(tximport$abundance) %in% sample_names])
  colnames(tximport$abundance)<- sample_names
  tximport$counts <-data.matrix(data.frame(tximport$counts)[colnames(tximport$counts) %in% sample_names])
  colnames(tximport$counts)<- sample_names
  tximport$length <- data.matrix(data.frame(tximport$length)[colnames(tximport$length) %in% sample_names])
  colnames(tximport$length)<- sample_names
  return(tximport)
}

#we create a function to calculate z-score
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

#Functions to wrap labels of barplot 
## Core wrapping function
wrap.it <- function(x, len)
{ 
  sapply(x, function(y) paste(strwrap(y, len), 
                              collapse = "\n"), 
         USE.NAMES = FALSE)
}


## Call this function with a list or vector
wrap.labels <- function(x, len)
{
  if (is.list(x))
  {
    lapply(x, wrap.it, len)
  } else {
    wrap.it(x, len)
  }
}
