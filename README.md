# DEAPP: a web implementation of differential expression analysis using shiny

##Introduction
This DE analysis App is developed to 1). conduct differential expression (DE) analysis with edgeR, limma-voom, and DESeq2, and 2). cross-validate the DE analysis results with these 3 different DE analysis methods.

##Input files
The input of this App is 2 plain text files, one includes the count results of all tags with respect to each sample, another one includes the experimental group factor information for each sample. 

###Test data
Initial DE analysis is based on the provided RNA-seq (http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0031229) pre-analyais count results. Detailed information about these 2 input files can be seen in the "Introduction" tab of this App.

##App execution
To run this App, please follow the steps as below:

1. Open R or RStudio with R version >= 3.1.2
2. Set your working directory to where the App programs (ui.R and server.R) are located. 
   For example, I downloaded the App package to my desktop, then I need to set the working directory in R by setwd("~/Desktop/RNAseq-DE").
3. Install all depended CRAN R packages (shiny, shinydashboard, DT, ggplot2) and R Bioconductor packages (edgeR, limma, DESeq2) by sourcing the R installation program with source("install/prep.R").
If all depended packages are successfully installed, logical value “TRUE” should be returned and printed on the R Console pane for each package after sourcing the installation program. 
   If not, please manually install all above requested depended packages in R, restart R/RStudio and source the package installation code again with source("install/prep.R") to make sure that logical value “TRUE” will be returned and printed on the R Console. If some packages still cannot be installed, try to update your R and corresponding Bioconductor version, and repeat this step again until you can see logical value “TRUE” are returned and printed on the R Console for each requested R package.
4. run this DE analysis App by shiny::runApp()
5. A web page will be open in you Browser to display all DE analysis results with initial provided test data, and the "Data Input" tab could allow you to uploaded you own count results together with the experimental factor information for fast efficient DE analysis with 3 different analysis methods.

##Feedback
If you have further questions regarding this App, please contact Yan Li at yli22@bsd.uchicago.edu from the bioinformatics core at the Center for Research Informatics (CRI), University of Chicago.
