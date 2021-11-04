#######################################################################
#                        Salmon Alignment                             #
#######################################################################

####### Packages Needed
library(DESeq2)
library(tximport)
library(org.Hs.eg.db)

####### Parameters
analysis.dir <- "./RNA-Seq_2D_3D_ipsDifferentiation"
sample.data <- read.table(file = "./Sample.Info.AllSamples.3.csv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
file.names <- sample.data$File.Name
# Salmon alignement and quantification parameters
salmonDir <- "/soft/bio/salmon-1.1.0/bin/salmon"
file1.suffix <- "_1.fastq.gz"
file2.suffix <- "_2.fastq.gz"
fastqdir <- "/imppc/labs/eslab/mmagallon/Projects/RNA-Seq-timecourse.2/Data/"
transcript.index <- "./Results/Salmon/salmon_indexes_UCSC_hg38"
output.suffix <- "_quant"
output.quants <- "./Results/Salmon/"
threads <- 8

source(file ="./utils.R")

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


###################          Salmon alignment           #####################

# ## Executing Salmon alignement by Selective allignement and quantification
# sn <- file.names[1]
for(sn in file.names){
  salmonAlignment(sample.name = sn, salmonDir,
                  file1.suffix = file1.suffix,
                  file2.suffix = file2.suffix,
                  fastqdir = fastqdir,
                  transcript.index = transcript.index,
                  output.suffix = output.suffix,
                  output.quants = output.quants,
                  threads = threads
  )
}
now.msg("   Salmon done")
