# DEApp: an interactive web application of differential expression analysis

##What is DEApp
This DE analysis interactive web application (App) is developed in R with Shiny, aiming to 1). conduct differential expression (DE) analysis with edgeR, limma-voom, and DESeq2, and 2). cross-validate the DE analysis results with these 3 different DE analysis methods based on your own provided input files.

##How To Use It
To run this App after downloading this App from github to your local PC, please follow the steps as below:

1. If 'Download ZIP', upzip the downloaded file to the Desktop, if 'Clone in Desktop' via git, the whole App folder ("DEAPP-shiny") should be downloaded to the Desktop. 

2. Open R or RStudio (if installed) with R version >= 3.2.

3. Set your working directory to where this App are unzipped or downloaded. 

   For example, I downloaded this App via 'Download ZIP' and unzipped it to the Desktop, then I need to set the working directory in R by `setwd("~/Desktop/DEAPP-shiny-master")`.
   
4. Install all depended CRAN R packages (shiny, shinydashboard, DT, ggplot2) and R Bioconductor packages (edgeR, limma, DESeq2) by sourcing the R installation program with `source("install/prep.R")`. 

    If all depended packages are successfully installed, logical value “TRUE” should be returned and printed on the R Console pane for each package after sourcing the installation program. Otherwise, please check whether you have the R version >=3.2 installed.
    
5. run this App by `shiny::runApp()`

6. A web page will be open in you Browser to display all DE analysis results with initial provided test data, and the "Data Input" tab could allow you to uploaded you own count results together with the experimental factor information for fast efficient DE analysis with 3 different analysis methods.

##Input Files
The input of this App is 2 files in '.txt' or '.csv' format.

1. 'Raw Count Data' : includes the count results of all tags with respect to each sample.  
2. 'Meta-data Table' : includes the experimental group factor information for each sample.

**Demo input files**: 3 sets of data were used to test this App, they are under 'data' folder, named as:

1. TestData-featureCount.txt (input data 1) + TestData-featureCount-meta.txt (input data 2)
2. pnas-count_singleFactor.txt (input data 1) + pnas-count_singleFacotr-meta.txt (input data 2)
3. ReadCounts-Chen-edgeRSpringer-multiFactor.csv (input data 1) + ReadCounts-Chen-edgeRSpringer-multiFactor-meta.csv (input data 2)

Among them, the first and the second data's input files are saved in tab delimited text format, and they both are from single-factor experiment; whereas the third data's input files are saved as csv comma delimited format set, and it is from multi-factor experiment. The 'Raw Count Data' file of the first test data set is pre-processed results with 'featureCount' based on the RNA-seq experiment published [here](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0031229), 
whereas 'Raw Count Data' input files of the second and the third data set are downloaded from [https://sites.google.com/site/davismcc/useful-documents](https://sites.google.com/site/davismcc/useful-documents) and [http://bioinf.wehi.edu.au/edgeRSpringer/](http://bioinf.wehi.edu.au/edgeRSpringer/) respectively.

##Feedback
If you have further questions or suggestions regarding this App, please contact Yan Li at yli22@bsd.uchicago.edu from the bioinformatics core at the Center for Research Informatics (CRI), biological science division (BSD), University of Chicago.
