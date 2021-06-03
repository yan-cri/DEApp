packages <- c("shinydashboard", "DT","shiny", "ggplot2", "gplots")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

BCpackages <- c("edgeR", "DESeq2", "limma")
if (length(setdiff(BCpackages, rownames(installed.packages()))) > 0) {
  # source("http://bioconductor.org/biocLite.R")
  # biocLite(setdiff(BCpackages, rownames(installed.packages())))
  BiocManager::install(setdiff(BCpackages, rownames(installed.packages())))
}

sapply(c(packages, BCpackages), require, character.only=T)

print(sapply(c(packages, BCpackages), require, character.only=T))
