#define the repository
nz_repo = c("https://cran.stat.auckland.ac.nz/")

if (!require("tidyverse")) {
  install.packages("tidyverse", repos = nz_repo)
}
library("tidyverse")

if (!require("openxlsx")) {
  install.packages("openxlsx", repos = nz_repo)
}
library("openxlsx")

#dependency of DEP
if (!requireNamespace("mzR", quietly=TRUE)) {
  if (!requireNamespace('BiocManager', quietly = TRUE))
    install.packages('BiocManager', repos = nz_repo)
  BiocManager::install("mzR")
}
library("mzR")

#dependency of DEP
if (!require("iterators")) {
  install.packages("iterators", repos = nz_repo)
}
library("iterators")

#dependency of DEP
if (!require("png")) {
  install.packages("png", repos = nz_repo)
}
library("png")

#dependency of DEP
if (!require("rjson")) {
  install.packages("rjson", repos = nz_repo)
}
library("rjson")

#dependency of DEP
if (!require("hexbin")) {
  install.packages("hexbin", repos = nz_repo)
}
library(hexbin)

#dependency of DEP
if (!require("magick")) {
  install.packages("magick", repos = nz_repo)
}
library(magick)

#DEP
if (!requireNamespace("DEP", quietly=TRUE)) {
  if (!requireNamespace('BiocManager', quietly = TRUE))
    install.packages('BiocManager', repos = nz_repo)
  BiocManager::install("DEP")
}
library("DEP")

if (!requireNamespace("mixOmics", quietly=TRUE)) {
  if (!requireNamespace('BiocManager', quietly = TRUE))
    install.packages('BiocManager', repos = nz_repo)
  BiocManager::install("mixOmics")
}
library("mixOmics")

if (!require("cowplot")) {
  install.packages("cowplot", repos = nz_repo)
}
library("cowplot")

if (!require("ComplexUpset")) {
  install.packages("ComplexUpset", repos = nz_repo)
}
library("ComplexUpset")

if (!require("kableExtra")) {
  install.packages("kableExtra", repos = nz_repo)
}
library("kableExtra")

if (!require("SummarizedExperiment")) {
  if (!requireNamespace('BiocManager', quietly = TRUE))
    install.packages('BiocManager', repos = nz_repo)
  BiocManager::install("SummarizedExperiment")
}
library("SummarizedExperiment")

if (!require("GGally")) {
  install.packages("GGally", repos = nz_repo)
}
library(GGally)

# if (!require("MSnbase")) {
#   if (!requireNamespace('BiocManager', quietly = TRUE))
#     install.packages('BiocManager', repos = nz_repo)
#   BiocManager::install("MSnbase")
# }
# library("MSnbase")
# 
# if (!require("MSnset.utils")) {
#   if (!requireNamespace('BiocManager', quietly = TRUE))
#     install.packages('BiocManager', repos = nz_repo)
#   BiocManager::install("MSnset.utils")
# }
# library("MSnset.utils")


###############################################################################
#diann package
###############################################################################
if (!require("diann")) {
  if (!require("devtools")) {
    install.packages("devtools")  
  }
  library(devtools)
  install_github("https://github.com/vdemichev/diann-rpackage")
}
library(diann)


if (!require("limma")) {
  if (!requireNamespace('BiocManager', quietly = TRUE))
    install.packages('BiocManager', repos = nz_repo)
  BiocManager::install("limma")
}
library(limma)


#ComplexHeatmap
if (!require("ComplexHeatmap")) {
  if (!requireNamespace('BiocManager', quietly = TRUE))
    install.packages('BiocManager', repos = nz_repo)
  BiocManager::install("ComplexHeatmap")
}
library("ComplexHeatmap")

#ClusterProfiler
if (!require("clusterProfiler")) {
  if (!requireNamespace('BiocManager', quietly = TRUE))
    install.packages('BiocManager', repos = nz_repo)
  BiocManager::install("clusterProfiler")
  BiocManager::install("org.Hs.eg.db")
}

library("clusterProfiler")
library("org.Hs.eg.db")

set.seed(223422)

###############################################################################
#Some common functions
###############################################################################

subset_summarizedexperiment <- function(se, metadata) {
  #subsets a summarized experiment
  #accepts a SE, and a metadata table
  #metadata table must contain the 'label' column
  #reconstructs the SE
  
  #queries the incomming se
  id_name <- get_df_wide(se) %>%
    dplyr::select(ID, name)
  
  temp_assay <- get_df_wide(se) %>%
    dplyr::select(-c(ID, name)) %>%
    #subset the columns
    dplyr::select(matches(metadata$label))
  
  #transform the temp_assay
  temp_assay <- 2^temp_assay %>%
    #convert to df
    data.frame() %>%
    #put back the ids
    mutate(ID) = id_name$ID %>%
    #put back the names
    mutate(name) = id_name$name %>%
    #move the columns
    dplyr::relocate(ID) %>%
    dplyr::relocate(name, .after(ID))
  
  #get the sample column positions
  sample_column_positions <- seq(3, length(colnames(temp_assay)))
  
  #convert back to a SE
  SE <- make_se(temp_assay, sample_column_positions, metadata)
  return(SE)
}





