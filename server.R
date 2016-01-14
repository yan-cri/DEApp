## This App is used for RNA-seq DE analysis with different methods together with method comparison
## Developed by Yan Li, last update on June, 2015
packages <- c("shinydashboard", "DT","shiny", "ggplot2")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

BCpackages <- c("edgeR", "DESeq2", "limma")
if (length(setdiff(BCpackages, rownames(installed.packages()))) > 0) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(setdiff(BCpackages, rownames(installed.packages())))
}

sapply(c(packages, BCpackages), require, character.only=T)

print(sapply(c(packages, BCpackages), require, character.only=T))

trim <- function (x) gsub("^\\s+|\\s+$", "", x)


shinyServer(function(input, output, session) {
  dataObs <- reactiveValues(
    orgCount = read.delim(paste(getwd(),"data/pnas-count_singleFactor.txt",sep="/"), header=T, row.names=1),
    orgMeta = read.delim(paste(getwd(),"/data/pnas-count_singleFactor-meta.txt",sep=""), header=T) 
    )
  
  observeEvent(input$MultiSubmit, {
    if (is.null(input$countFileMulti) & is.null(input$metaTabMulti)) {
      dataObs$orgCount <- read.delim(paste(getwd(),"data/ReadCounts-Chen-edgeRSpringer-multiFactor.csv",sep="/"), header=T, sep=input$coutFileSepMulti, row.names=1)
      dataObs$orgMeta <- read.delim(paste(getwd(),"/data/ReadCounts-Chen-edgeRSpringer-multiFactor-meta.csv",sep=""), header=T, sep=input$metaSepMulti)
    } else if (!is.null(input$countFileMulti) & is.null(input$metaTabMulti)) {
      dataObs$orgCount <- read.delim(paste(getwd(),"data/ReadCounts-Chen-edgeRSpringer-multiFactor.csv",sep="/"), header=T, sep=input$coutFileSepMulti, row.names=1)
      dataObs$orgMeta <- NULL
    } else if (is.null(input$countFileMulti) & !is.null(input$metaTabMulti)) {
      dataObs$orgCount <- NULL
      dataObs$orgMeta <- read.delim(paste(getwd(),"/data/ReadCounts-Chen-edgeRSpringer-multiFactor-meta.csv",sep=""), header=T, sep=input$metaSepMulti)
    } else if (!is.null(input$countFileMulti) & !is.null(input$metaTabMulti)) {
      dataObs$orgCount <- try(read.delim(input$countFileMulti$datapath, header=T, sep=input$coutFileSepMulti, row.names=1 ), T)
      dataObs$orgMeta <- try(read.delim(input$metaTabMulti$datapath, header=T, sep=input$metaSepMulti), T)
    }
    #print(input$dataInputMulti == "dataInputMulti")
    #print("===*****====")
    #print(head(dataObs$orgMeta))
    #print("===*****====")
  })
  
  observeEvent(input$dataSubmit, {
    if (is.null(input$countFile) & is.null(input$metaTab) ) {
      dataObs$orgCount <- read.delim(paste(getwd(),"data/pnas-count_singleFactor.txt",sep="/"), header=T, sep=input$coutFileSep, row.names=1)
      dataObs$orgMeta <- read.delim(paste(getwd(),"/data/pnas-count_singleFactor-meta.txt",sep=""), header=T, sep=input$metaSep)
    } else if (!is.null(input$countFile) & is.null(input$metaTab) ) {
      dataObs$orgCount <- read.delim(paste(getwd(),"data/pnas-count_singleFactor.txt",sep="/"), header=T, sep=input$coutFileSep, row.names=1)
      dataObs$orgMeta <- NULL
    } else if (is.null(input$countFile) & !is.null(input$metaTab) ) {
      dataObs$orgCount <- NULL
      dataObs$orgMeta <- read.delim(paste(getwd(),"/data/pnas-count_singleFactor-meta.txt",sep=""), header=T, sep=input$metaSep)
    } else if (!is.null(input$countFile) & !is.null(input$metaTab) ){
      dataObs$orgCount <- try(read.delim(input$countFile$datapath, header=T, sep=input$coutFileSep, row.names=1 ), T)
      dataObs$orgMeta <- try(read.delim(input$metaTab$datapath, header=T, sep=input$metaSep), T) 
    }   
    #print(input$dataInputSingle == "dataInputSingle")
    #print("====*****")
    #print(head(dataObs$orgMeta))
    #print("====*****")
  })
  
  metaUpdateMulti <- eventReactive(input$MultiSubmit, {
    if (is.null(input$countFileMulti) & is.null(input$metaTabMulti)) {
      dataObs$orgCount <- read.delim(paste(getwd(),"data/ReadCounts-Chen-edgeRSpringer-multiFactor.csv",sep="/"), header=T, sep=input$coutFileSepMulti, row.names=1)
      dataObs$orgMeta <- read.delim(paste(getwd(),"/data/ReadCounts-Chen-edgeRSpringer-multiFactor-meta.csv",sep=""), header=T, sep=input$metaSepMulti)
    } else if (!is.null(input$countFileMulti) & is.null(input$metaTabMulti)) {
      dataObs$orgCount <- read.delim(paste(getwd(),"data/ReadCounts-Chen-edgeRSpringer-multiFactor.csv",sep="/"), header=T, sep=input$coutFileSepMulti, row.names=1)
      dataObs$orgMeta <- NULL
    } else if (is.null(input$countFileMulti) & !is.null(input$metaTabMulti)) {
      dataObs$orgCount <- NULL
      dataObs$orgMeta <- read.delim(paste(getwd(),"/data/ReadCounts-Chen-edgeRSpringer-multiFactor-meta.csv",sep=""), header=T, sep=input$metaSepMulti)
    } else if (!is.null(input$countFileMulti) & !is.null(input$metaTabMulti)) {
      dataObs$orgCount <- try(read.delim(input$countFileMulti$datapath, header=T, sep=input$coutFileSepMulti, row.names=1 ), T)
      dataObs$orgMeta <- try(read.delim(input$metaTabMulti$datapath, header=T, sep=input$metaSepMulti), T)
    }
    list(meta=as.data.frame(dataObs$orgMeta), count=as.data.frame(dataObs$orgCount) )
  })
  
  metaUpdate <- eventReactive(input$dataSubmit, {
    if (is.null(input$countFile) & is.null(input$metaTab) ) {
      dataObs$orgCount <- read.delim(paste(getwd(),"data/pnas-count_singleFactor.txt",sep="/"), header=T, sep=input$coutFileSep, row.names=1)
      dataObs$orgMeta <- read.delim(paste(getwd(),"/data/pnas-count_singleFactor-meta.txt",sep=""), header=T, sep=input$metaSep)
    } else if (!is.null(input$countFile) & is.null(input$metaTab) ) {
      dataObs$orgCount <- read.delim(paste(getwd(),"data/pnas-count_singleFactor.txt",sep="/"), header=T, sep=input$coutFileSep, row.names=1)
      dataObs$orgMeta <- NULL
    } else if (is.null(input$countFile) & !is.null(input$metaTab) ) {
      dataObs$orgCount <- NULL
      dataObs$orgMeta <- read.delim(paste(getwd(),"/data/pnas-count_singleFactor-meta.txt",sep=""), header=T, sep=input$metaSep)
    } else if (!is.null(input$countFile) & !is.null(input$metaTab) ){
      dataObs$orgCount <- try(read.delim(input$countFile$datapath, header=T, sep=input$coutFileSep, row.names=1 ), T)
      dataObs$orgMeta <- try(read.delim(input$metaTab$datapath, header=T, sep=input$metaSep), T) 
    }   
    list(meta=as.data.frame(dataObs$orgMeta), count=as.data.frame(dataObs$orgCount) )
  })
  
  progress <- reactiveValues(time=shiny::Progress$new())
  
  observeEvent(input$dataSubmit, {     
   progress$time$set(message = "Data input (single-factor)", value = 0)
   progress$time$set(value = 0.5, detail = "processing 50%")
  })
  
  observeEvent(input$MultiSubmit, { 
   progress$time$set(message = "Data input (multi-factor)", value = 0)
   progress$time$set(value = 0.5, detail = "processing 50%")
  })
  
  observeEvent(input$rmlow, { 
    progress$time$set(message = "Filtering", value = 0)
    progress$time$set(value = 0.5, detail = "processing 50%")
  })
  
  observeEvent(input$decompAnalysis, {
    progress$time$set(message = "Comparison analysis", value = 0)
    progress$time$set(value = 0.2, detail = "processing 20%")
  })
  
  observeEvent(input$edgerdeAnalysis, { 
    progress$time$set(message = "edgeR analysis", value = 0)
    progress$time$set(value = 0.2, detail = "processing 20%")
  })
  
  observeEvent(input$voomdeAnalysis, { 
    progress$time$set(message = "Limma-voom analysis", value = 0)
    progress$time$set(value = 0.2, detail = "processing 20%")
  })
  
  observeEvent(input$deseq2deAnalysis, {
    progress$time$set(message = "DESeq2 analysis", value = 0)
    progress$time$set(value = 0.2, detail = "processing 20%")
  })
  
  ##Reactive expression object for original row count
  datareactive <- reactive ({              
    org.counts <- dataObs$orgCount
    metadata <- dataObs$orgMeta
    metadata <- metadata[match(colnames(org.counts), metadata$Sample),]
    #print(head(org.counts))
    #print(head(metadata))
    #print("*********")
    if (input$decompAnalysis) {
      progress$time$set(message = "Comparison analysis", value = 0)
      progress$time$set(value = 0.2, detail = "processing 20%")
    }          
    if (dim(metadata)[2]>2) {
      groupinfo <- metadata[,2]
      for (i in 3:length(metadata[1,])) {
        groupinfo <- paste(groupinfo, metadata[,i], sep=".")
      }
      Group <- factor(groupinfo)
    } else {
      Group <- factor(metadata[,2])
    }
    
    org.count <- DGEList(counts=org.counts, group=Group)
    org.count$samples$lib.size <- colSums(org.count$counts)
    org.count <- calcNormFactors(org.count, lib.size=T, method="TMM")
    
    org.count
  })
  
  
  #Reactive expression object for the counts remove low expression tags
  rmlowReactive <- reactive({    
    if(input$edgerdeAnalysis) { 
      progress$time$set(message = "edgeR analysis", value = 0)
      progress$time$set(value = 0.2, detail = "processing 20%")
    }
    
    dge.count <- calcNormFactors(datareactive())
    cpm.count <- cpm(dge.count$counts)
    if (as.numeric(input$gThreshold) > length(colnames(datareactive()$counts)) ) stop("Cutoff sample number exceeds the total number of samples")
    keep = rowSums(cpm.count > as.numeric(input$cpmVal) ) >= as.numeric(input$gThreshold)
    dge.count.rmlow <- dge.count[keep,]
    dge.count.rmlow$samples$lib.size <- colSums(dge.count.rmlow$counts)
    dge.count.rmlow <- calcNormFactors(dge.count.rmlow, lib.size=T, method="TMM")
    dge.count.rmlow
  })
  
  #Reactive expression object for edgeR glm dispersion estimation
  edgerDispersionEst <- reactive({  
    
    dge <- rmlowReactive()
    if(!is(dge,"DGEList")) stop("Dispersion estimation input must be a DGEList.")
    Group <- rmlowReactive()$samples$group
    f <- factor(Group, levels=levels(Group))
    design <- model.matrix(~0+f)
    rownames(design) <- rownames(dge$samples)
    colnames(design) <- levels(Group)
    dge <- estimateDisp(dge, design)
    if(input$decompAnalysis) {
      progress$time$set(message = "Comparison analysis", value = 0.2)
      progress$time$set(value = 0.4, detail = "processing 40%")
    }   
    if(input$edgerdeAnalysis) { 
      progress$time$set(message = "edgeR analysis", value = 0.2)
      progress$time$set(value = 0.4, detail = "processing 40%")
    }
    dge
  })
  
  #Reactive expression object for edgeR glm fit
  edgerglmFit <- reactive({  
    Group <- edgerDispersionEst()$samples$group
    f <- factor(Group, levels=levels(Group))
    design <- model.matrix(~0+f)
    rownames(design) <- rownames(edgerDispersionEst()$samples)
    colnames(design) <- levels(Group)
    glm.fit <- glmFit(edgerDispersionEst(), design)
    glm.fit
  })
  
  #Reactive expression object for edgeR glmLRT result
  edgerDEres <- reactive({    
    Group <- edgerDispersionEst()$samples$group
    f <- factor(Group, levels=levels(Group))
    design <- model.matrix(~0+f)
    rownames(design) <- rownames(edgerDispersionEst()$samples)
    colnames(design) <- levels(Group)
    if (as.character(trim(input$edgercompGroup2)) =="" | as.character(trim(input$edgercompGroup1)) == "") {
      stop("Please select 2 groups levels for DE analysis!")
    } else if (as.character(trim(input$edgercompGroup2)) == as.character(trim(input$edgercompGroup1)) ) {
      stop("Group level 1 and level 2 are the same, please select 2 different group levels for DE analysis!")
    } else if( (! as.character(trim(input$edgercompGroup2)) %in% levels(Group)) | (! as.character(trim(input$edgercompGroup1)) %in% levels(Group)) ) {
      stop("Group level 1 or level 2 are not in the available group levels!")
    }     
    comp <- makeContrasts(contrasts=paste(as.character(trim(input$edgercompGroup2)), as.character(trim(input$edgercompGroup1)), sep="-"), levels=design )
    test <- glmLRT(edgerglmFit(), contrast=comp)
    test
  })
  
  #Reactive expression object for edgeR desideTestsDEG results
  edgerDEfilter <- reactive({
    observeEvent(input$edgerdeAnalysis, { 
      progress$time$set(value = 1, detail = "processing 100%")
    })
    
    Group <- edgerDispersionEst()$samples$group
    if (as.character(trim(input$edgercompGroup2)) =="" | as.character(trim(input$edgercompGroup1)) == "") {
      stop("Please select 2 group levels for DE analysis!")
    } else {
      if (as.character(trim(input$edgercompGroup2)) == as.character(trim(input$edgercompGroup1)) ) {
        stop("Group level 1 and level 2 are the same, please select 2 different group levels for DE analysis!")
      } else if( (! as.character(trim(input$edgercompGroup2)) %in% levels(Group)) | (! as.character(trim(input$edgercompGroup1)) %in% levels(Group)) ) 
      {stop("Group level 1 or level 2 are not in the available group levels!")
      }
    } 
    
    fcval <- as.numeric(input$edgerfc)
    fdr <- as.numeric(input$edgerfdr)
    if (input$edgerP == 'normp') {
      tp <- topTags(edgerDEres(), n=Inf, adjust.method="BH", sort.by="none")
      filter <- decideTestsDGE(edgerDEres(), adjust.method="none", p.value=fdr, lfc=log2(fcval) )
    } else if (input$edgerP == 'fdrp') {
      tp <- topTags(edgerDEres(), n=Inf, adjust.method="BH", sort.by="none")
      filter <- decideTestsDGE(edgerDEres(), adjust.method="BH", p.value=fdr, lfc=log2(fcval) )
    }
    tp$table$filter <- as.factor(filter)
    tp
    
  }) 
  
  ##Reactive expression object for limma-voom voom function
  voomRes <- reactive({      
    if (input$voomdeAnalysis) { 
      progress$time$set(message = "Limma-voom analysis", value = 0)
      progress$time$set(value = 0.2, detail = "processing 20%")
    }
    Group <- datareactive()$samples$group
    f <- factor(Group, levels=levels(Group))
    design <- model.matrix(~0+f)  
    rownames(design) <- rownames(rmlowReactive()$samples)
    colnames(design) <- levels(Group)
    par(mar=c(0,0,0,0))
    v <- voom(rmlowReactive(), design=design, plot=T, normalize="quantile")
    fit <- lmFit(v, design)
    if (input$decompAnalysis) {
      progress$time$set(message = "Comparison analysis", value = 0.4)
      progress$time$set(value = 0.6, detail = "processing 60%")
    }  
    if (input$voomdeAnalysis) { 
      progress$time$set(message = "Limma-voom analysis", value = 0.2)
      progress$time$set(value = 0.6, detail = "processing 60%")
    }
    fit
  })
  
  ##Reactive expression object for limma-voom DE analysis after contrast.fit
  voomDEres <- reactive({
    Group <- datareactive()$samples$group
    f <- factor(Group, levels=levels(Group))
    design <- model.matrix(~0+f)
    rownames(design) <- rownames(rmlowReactive()$samples)
    colnames(design) <- levels(Group)
    if (as.character(trim(input$voomcompGroup2)) =="" | as.character(trim(input$voomcompGroup1)) == "") {
      stop("Please select 2 groups levels for DE analysis!")
    } else if (as.character(trim(input$voomcompGroup1)) == as.character(trim(input$voomcompGroup2)) ) {
      stop("Group level 1 and level 2 are the same, please select 2 different group levels for DE analysis!")
    } else if( (! as.character(trim(input$voomcompGroup2)) %in% levels(Group)) | (! as.character(trim(input$voomcompGroup1)) %in% levels(Group)) ) {
      stop("Group level 1 or level 2 are not in the available group levels!")
    } 
    
    comp <- makeContrasts(contrasts=paste(as.character(trim(input$voomcompGroup2)), as.character(trim(input$voomcompGroup1)), sep="-"), levels=design )
    
    contrast.fit <- contrasts.fit(voomRes(), contrasts=comp)
    contrast.fit <- eBayes(contrast.fit)  
    contrast.fit
  })
  
  ##Reactive expression object for limma-voom decideTests results
  voomDEfilter <- reactive({
    Group <- datareactive()$samples$group
    if (as.character(trim(input$voomcompGroup2)) =="" | as.character(trim(input$voomcompGroup1)) == "") {
      stop("Please select 2 groups levels for DE analysis!")
    } else if (as.character(trim(input$voomcompGroup1)) == as.character(trim(input$voomcompGroup2)) ) {
      stop("Group level 1 and level 2 are the same, please select 2 different group levels for DE analysis!")
    } else if( (! as.character(trim(input$voomcompGroup2)) %in% levels(Group)) | (! as.character(trim(input$voomcompGroup1)) %in% levels(Group)) ) {
      stop("Group level 1 or level 2 are not in the available group levels!")
    } 
    
    voomfcval <- as.numeric(input$voomfc)
    voomfdr <- as.numeric(input$voomfdr)
    if (input$voomP == 'normp') {
      tpvoom <- topTable(voomDEres(), number=Inf, adjust="BH", sort.by="none") 
      filtervoom <- decideTests(voomDEres(), adjust.method="none", p.value=voomfdr, lfc=log2(voomfcval))
    } else if (input$voomP == 'fdrp') {
      tpvoom <- topTable(voomDEres(), number=Inf, adjust="BH", sort.by="none") 
      filtervoom <- decideTests(voomDEres(), adjust.method="BH", p.value=voomfdr, lfc=log2(voomfcval))
    }
    observeEvent(input$voomdeAnalysis, { 
      progress$time$set(value = 1, detail = "processing 100%")
    })
    tpvoom$filter <- as.factor(filtervoom)
    tpvoom
  })
  
  ##Reactive expression object for DEseq2 DESeq analysis results
  deseq2Res <- reactive({    
    if(input$deseq2deAnalysis) { 
      progress$time$set(message = "DESeq2 analysis", value = 0)
      progress$time$set(value = 0.2, detail = "processing20%")
    }
    Group <- datareactive()$samples$group
    colData <- data.frame(Group)
    dds <- DESeqDataSetFromMatrix(rmlowReactive()$count, colData=colData, design=formula(~Group) )
    dds <- DESeq(dds, test="Wald") 
    if(input$decompAnalysis) {
      progress$time$set(message = "Comparison analysis", value = 0.4)
      progress$time$set(value = 0.7, detail = "processing 70%")}
    if(input$deseq2deAnalysis) { 
      progress$time$set(message = "DESeq2 analysis", value = 0.4)
      progress$time$set(value = 0.7, detail = "processing 70%")
    }
    dds
  })
  
  ##Reactive expression object for DESeq2 results()
  deseq2DEres <- reactive({
    Group <- datareactive()$samples$group
    if (as.character(trim(input$deseq2compGroup2)) =="" | as.character(trim(input$deseq2compGroup1)) == "") {
      stop("Please select 2 groups levels for DE analysis!")
    } else if (as.character(trim(input$deseq2compGroup1)) == as.character(trim(input$deseq2compGroup2)) ) {
      stop("Group level 1 and level 2 are the same, please select 2 different group levels for DE analysis!")
    } else if( (! as.character(trim(input$deseq2compGroup2)) %in% levels(Group)) | (! as.character(trim(input$deseq2compGroup1)) %in% levels(Group)) ) {
      stop("Group level 1 or level 2 are not in the available group levels!")
    } 
    res <- results(deseq2Res(), contrast=c("Group", as.character(trim(input$deseq2compGroup2)), as.character(trim(input$deseq2compGroup1))), format="DataFrame")
    as.data.frame(res)
  })
  
  ##Reactive expression object for DESeq2 decideTests results
  deseq2DEfilter <- reactive({
    deseq2fcval <- as.numeric(input$deseq2fc)
    deseq2fdr <- as.numeric(input$deseq2fdr)
    deseq2res <- deseq2DEres()
    deseq2res$filter <- 0
    if (input$deseq2P == 'normp') {
      deseq2res$filter[deseq2res$pvalue < deseq2fdr & deseq2res$log2FoldChange > log2(deseq2fcval)] <- 1
      deseq2res$filter[deseq2res$pvalue < deseq2fdr & deseq2res$log2FoldChange < -log2(deseq2fcval)] <- -1
    } else if (input$deseq2P== 'fdrp') {
      deseq2res$filter[deseq2res$padj < deseq2fdr & deseq2res$log2FoldChange > log2(deseq2fcval)] <- 1
      deseq2res$filter[deseq2res$padj < deseq2fdr & deseq2res$log2FoldChange < -log2(deseq2fcval)] <- -1
    }
    if (input$deseq2deAnalysis) { 
      progress$time$set(value = 1, detail = "processing 100%")
    }
    deseq2res$filter <- as.factor(deseq2res$filter)
    deseq2res
  })
  
  ##Reactive expression object of edgeR for methods comparison
  edgerDecomp <- reactive({
    Group <- edgerDispersionEst()$samples$group
    f <- factor(Group, levels=levels(Group))
    design <- model.matrix(~0+f)
    rownames(design) <- rownames(edgerDispersionEst()$samples)
    colnames(design) <- levels(Group)
    
    if (as.character(trim(input$decompGroup2)) =="" | as.character(trim(input$decompGroup1)) == "") {
      stop("Please select 2 groups levels for DE analysis!")
    } else if (as.character(trim(input$decompGroup1)) == as.character(trim(input$decompGroup2)) ) {
      stop("Group level 1 and level 2 are the same, please select 2 different group levels for DE analysis!")
    } else if( (! as.character(trim(input$decompGroup2)) %in% levels(Group)) | (! as.character(trim(input$decompGroup1)) %in% levels(Group)) ) {
      stop("Group level 1 or level 2 are not in the available group levels!")
    } 
    
    comp <- makeContrasts(contrasts=paste(as.character(trim(input$decompGroup2)), as.character(trim(input$decompGroup1)), sep="-"), levels=design )
    compRes <- glmLRT(edgerglmFit(), contrast=comp)
    fcval <- as.numeric(input$decompfc)
    fdr <- as.numeric(input$decompfdr)
    if (input$decompP == 'normp') {
      tp <- topTags(compRes, n=Inf, adjust.method="BH", sort.by="none")
      filter <- decideTestsDGE(compRes, adjust.method="none", p.value=fdr, lfc=log2(fcval) )
    } else if (input$decompP == 'fdrp') {
      tp <- topTags(compRes, n=Inf, adjust.method="BH", sort.by="none")
      filter <- decideTestsDGE(compRes, adjust.method="BH", p.value=fdr, lfc=log2(fcval) )
    }
    tp$table$filter <- filter
    tp
  })
  
  ##Reactive expression object of limma-voom for methods comparison
  voomDecomp <- reactive({
    Group <- datareactive()$samples$group
    f <- factor(Group, levels=levels(Group))
    design <- model.matrix(~0+f)
    rownames(design) <- rownames(rmlowReactive()$samples)
    colnames(design) <- levels(Group)
    
    if (as.character(trim(input$decompGroup2)) =="" | as.character(trim(input$decompGroup1)) == "") {
      stop("Please select 2 groups levels for DE analysis!")
    } else if (as.character(trim(input$decompGroup1)) == as.character(trim(input$decompGroup2)) ) {
      stop("Group level 1 and level 2 are the same, please select 2 different group levels for DE analysis!")
    } else if( (! as.character(trim(input$decompGroup2)) %in% levels(Group)) | (! as.character(trim(input$decompGroup1)) %in% levels(Group)) ) {
      stop("Group level 1 or level 2 are not in the available group levels!")
    } 
    
    comp <- makeContrasts(contrasts=paste(as.character(trim(input$decompGroup2)), as.character(trim(input$decompGroup1)), sep="-"), levels=design )
    contrast.fit <- contrasts.fit(voomRes(), contrasts=comp)
    contrast.fit <- eBayes(contrast.fit)  
    
    voomfcval <- as.numeric(input$decompfc)
    voomfdr <- as.numeric(input$decompfdr)
    if (input$decompP == 'normp') {
      tpvoom <- topTable(contrast.fit, number=Inf, adjust="BH", sort.by="none") 
      filtervoom <- decideTests(contrast.fit, adjust.method="none", p.value=voomfdr, lfc=log2(voomfcval))
    } else if (input$decompP == 'fdrp') {
      tpvoom <- topTable(contrast.fit, number=Inf, adjust="BH", sort.by="none") 
      filtervoom <- decideTests(contrast.fit, adjust.method="BH", p.value=voomfdr, lfc=log2(voomfcval))
    }
    tpvoom$filter <- filtervoom
    tpvoom
  })
  
  ##Reactive expression object of DESeq2 for methods comparison
  deseq2Decomp <- reactive({
    
    Group <- datareactive()$samples$group
    if (as.character(trim(input$decompGroup2)) =="" | as.character(trim(input$decompGroup1)) == "") {
      stop("Please select 2 groups levels for DE analysis!")
    } else if (as.character(trim(input$decompGroup1)) == as.character(trim(input$decompGroup2)) ) {
      stop("Group level 1 and level 2 are the same, please select 2 different group levels for DE analysis!")
    } else if( (! as.character(trim(input$decompGroup2)) %in% levels(Group)) | (! as.character(trim(input$decompGroup1)) %in% levels(Group)) ) {
      stop("Group level 1 or level 2 are not in the available group levels!")
    }        
    
    deseq2res  <- results(deseq2Res(), contrast=c("Group", as.character(trim(input$decompGroup2)), as.character(trim(input$decompGroup1)) ), format="DataFrame")
    deseq2fcval <- as.numeric(input$decompfc)
    deseq2fdr <- as.numeric(input$decompfdr)
    deseq2res$filter <- 0
    if (input$decompP == 'normp') {
      deseq2res$filter[deseq2res$pvalue < deseq2fdr & deseq2res$log2FoldChange > log2(deseq2fcval)] <- 1
      deseq2res$filter[deseq2res$pvalue < deseq2fdr & deseq2res$log2FoldChange < -log2(deseq2fcval)] <- -1
    } else if (input$decompP== 'fdrp') {
      deseq2res$filter[deseq2res$padj < deseq2fdr & deseq2res$log2FoldChange > log2(deseq2fcval)] <- 1
      deseq2res$filter[deseq2res$padj < deseq2fdr & deseq2res$log2FoldChange < -log2(deseq2fcval)] <- -1
    }
    deseq2res
  })
  
  ##Reactive expression object of DE comparison results
  decompRes <- reactive({
    Group <- datareactive()$samples$group
    if (as.character(trim(input$decompGroup2)) =="" | as.character(trim(input$decompGroup1)) == "") {
      stop("Please select 2 groups levels for DE analysis!")
    } else if (as.character(trim(input$decompGroup1)) == as.character(trim(input$decompGroup2)) ) {
      stop("Group level 1 and level 2 are the same, please select 2 different group levels for DE analysis!")
    } else if( (! as.character(trim(input$decompGroup2)) %in% levels(Group)) | (! as.character(trim(input$decompGroup1)) %in% levels(Group)) ) {
      stop("Group level 1 or level 2 are not in the available group levels!")
    }   
    
    edgerDEgenename <- rownames(rbind(subset(edgerDecomp()$table, filter==1),subset(edgerDecomp()$table, filter==-1)))
    voomDEgenename <- rownames(rbind(subset(voomDecomp(), filter==1),subset(voomDecomp(), filter==-1)))
    deseq2DEgenename <- rownames(rbind(subset(deseq2Decomp(), filter==1),subset(deseq2Decomp(), filter==-1)))
    
    #res.matrix <- matrix(0, nrow=7, ncol=1)
    res.matrix <- data.frame(res=rep(0,7))
    rownames(res.matrix) <- c("edgeR", "limma-voom", "DESeq2", "edgeR & limma-voom", "edgeR & DESeq2", "limma-voom & DESeq2", "All")
    
    res.matrix[1,1] <- length(edgerDEgenename)
    res.matrix[2,1] <- length(voomDEgenename)
    res.matrix[3,1] <- length(deseq2DEgenename)
    res.matrix[4,1] <- length(intersect(edgerDEgenename,voomDEgenename))
    res.matrix[5,1] <- length(intersect(edgerDEgenename,deseq2DEgenename))
    res.matrix[6,1] <- length(intersect(voomDEgenename,deseq2DEgenename))
    res.matrix[7,1] <- length(intersect(intersect(edgerDEgenename,voomDEgenename),deseq2DEgenename))
    
    res.matrix
  })
  
  ##End of all reactive expression object
  ################################################################################################
  ################################################################################################
  ##tab Panel Data input
  output$countTabSamp <- renderTable({     
    org.counts <- read.delim(paste(getwd(),"www/dataInput1-exp.txt",sep="/"), header=T, row.names=1)
    org.counts
  }, align="l|cccccc")
  
  output$metaTabSamp <- renderTable({
    metadata <- read.delim(paste(getwd(),"/www/dataInput2-exp0.txt",sep=""), header=T)
    metadata
  }, include.rownames=F, align="llc")
  
  output$multimetaTabSamp22 <- renderTable({
    metadata <- read.delim(paste(getwd(),"/www/Multi-dataInput2-exp2.txt",sep=""), header=T)
    metadata
  }, include.rownames=F, align="llcccc")
  
  
  
  ################################################################################################
  ################################################################################################
  
  output$overallDataSummary <- renderTable({ 
    if(input$dataSubmit)
      isolate({
        #print("SSSS==number=======")
        #print(head(count))
        #print(dim(count))
        #print("=========")
        meta <- metaUpdate()$meta
        count <- metaUpdate()$count
        if (is.null(meta) & !is.null(count)) {stop("Please provide the corresponding input 2: Meta-data Table!")}
        else if (!is.null(meta) & is.null(count)) {stop("Please provide the input 1: Raw Count Data!")}
        else if (!is.null(meta) & !is.null(count)) {
          if (length(grep('Error',count[1]))==1) { 
            stop(paste(count[1]) )
          } else if (length(grep('Error',meta[1]))==1) {
            stop(paste(meta[1]) )
          } else {
            orgCount <- count
            orgMeta <- meta
            print(dim(orgCount)[2])
            print(dim(orgMeta)[1])
            print(sum(colnames(orgCount) == orgMeta[,1]))
            if (dim(orgCount)[2] != dim(orgMeta)[1] ) {stop("Input files do not correspond with each other. ")}
            #else if ( sum(colnames(orgCount) == orgMeta[,1]) !=dim(orgMeta)[1] ) {stop("Input files do not correspond with each other. ")}
            else if (dim(orgMeta)[2] < 2 ) {stop("Input 2 file format is wrong.") }
          }                 
        }
        
        no.samples <- length(colnames(count))
        no.gene <- dim(count)[1]
        observeEvent(input$dataSubmit, { 
          progress$time$set(value = 1, detail = "processing 100%")
        })
        res.summary <- rbind(no.samples, no.gene)
        rownames(res.summary) <- c("Samples", "Tags")
        colnames(res.summary) <- "Number"
        res.summary
      })
    
  },digits=0, align="l|c")
  
  output$sampleGroup <- renderTable({ 
    if(input$dataSubmit)
      isolate({      
        meta <- metaUpdate()$meta
        count <- metaUpdate()$count
        if (is.null(meta) & !is.null(count)) {stop("Please provide the corresponding input 2: Meta-data Table!")}
        else if (!is.null(meta) & is.null(count)) {stop("Please provide the input 1: Raw Count Data!")}
        else if (!is.null(meta) & !is.null(count)) {
          if (length(grep('Error',count[1]))==1) { 
            stop(paste(count[1]) )
          } else if (length(grep('Error',meta[1]))==1) {
            stop(paste(meta[1]) )
          } else {
            orgCount <- count
            orgMeta <- meta
            if (dim(orgCount)[2] != dim(orgMeta)[1] ) {stop("Input files do not correspond with each other. ")}
            #else if ( sum(colnames(orgCount) == orgMeta[,1]) !=dim(orgMeta)[1] ) {stop("Input files do not correspond with each other. ")}
            else if (dim(orgMeta)[2] < 2 ) {stop("Input 2 file format is wrong.") }
          }                  
        }
        
        groupinfo <- as.matrix(summary((datareactive())$samples$group))
        
        colnames(groupinfo) <- "No. in each group"
        groupinfo
      })
  },digits=0, align="l|c")
  
  output$sampleInfo <- renderText({ 
    if(input$dataSubmit)
      isolate({      
        #print("SSSSSS======")
        #print(metaUpdate()$meta)
        #print(head(metaUpdate()$count))
        meta <- metaUpdate()$meta
        count <- metaUpdate()$count
        if (is.null(meta) & !is.null(count)) {stop("Please provide the corresponding input 2: Meta-data Table!")}
        else if (!is.null(meta) & is.null(count)) {stop("Please provide the input 1: Raw Count Data!")}
        else if (!is.null(meta) & !is.null(count)) {
          if (length(grep('Error',count[1]))==1) { 
            stop(paste(count[1]) )
          } else if (length(grep('Error',meta[1]))==1) {
            stop(paste(meta[1]) )
          } else {
            orgCount <- count
            orgMeta <- meta
            if (dim(orgCount)[2] != dim(orgMeta)[1] ) {stop("Input files do not correspond with each other. ")}
            #else if ( sum(colnames(orgCount) == orgMeta[,1]) !=dim(orgMeta)[1] ) {stop("Input files do not correspond with each other. ")}
            else if (dim(orgMeta)[2] < 2 ) {stop("Input 2 file format is wrong.") }
          }                 
        }
        
        paste(as.character(meta[,1]), collapse=", " )
      })
  })
  
  output$sampleTitle <- renderText({ 
    if(input$dataSubmit)
      isolate({      
        #print("SSSSSS=Total=====")
        #print(metaUpdate()$meta)
        #print(head(metaUpdate()$count))
        
        meta <- metaUpdate()$meta
        count <- metaUpdate()$count
        if (is.null(meta) & !is.null(count)) {stop("Please provide the corresponding input 2: Meta-data Table!")}
        else if (!is.null(meta) & is.null(count)) {stop("Please provide the input 1: Raw Count Data!")}
        else if (!is.null(meta) & !is.null(count)) {
          if (length(grep('Error',count[1]))==1) { 
            stop(paste(count[1]) )
          } else if (length(grep('Error',meta[1]))==1) {
            stop(paste(meta[1]) )
          } else {
            orgCount <- count
            orgMeta <- meta
            if (dim(orgCount)[2] != dim(orgMeta)[1] ) {stop("Input files do not correspond with each other. ")}
            #else if ( sum(colnames(orgCount) == orgMeta[,1]) !=dim(orgMeta)[1] ) {stop("Input files do not correspond with each other. ")}
            else if (dim(orgMeta)[2] < 2 ) {stop("Input 2 file format is wrong.") }
          }        
          
        }
        
        no.samples <- length(colnames(count))
        paste("There are ", as.character(no.samples), " samples in the experiment:", sep="")
      })  
  })
  
  output$expDesign <- renderText({
    if(input$dataSubmit)
      isolate({
        meta <- metaUpdate()$meta
        count <- metaUpdate()$count
        if (is.null(meta) & !is.null(count)) {stop("Please provide the corresponding input 2: Meta-data Table!")}
        else if (!is.null(meta) & is.null(count)) {stop("Please provide the input 1: Raw Count Data!")}
        else if (!is.null(meta) & !is.null(count)) {
          if (length(grep('Error',count[1]))==1) { 
            stop(paste(count[1]) )
          } else if (length(grep('Error',meta[1]))==1) {
            stop(paste(meta[1]) )
          } else {
            orgCount <- count
            orgMeta <- meta
            if (dim(orgCount)[2] != dim(orgMeta)[1] ) {stop("Input files do not correspond with each other. ")}
            #else if ( sum(colnames(orgCount) == orgMeta[,1]) !=dim(orgMeta)[1] ) {stop("Input files do not correspond with each other. ")}
            else if (dim(orgMeta)[2] < 2 ) {stop("Input 2 file format is wrong.") }
          }                 
        }
        if (dim(meta)[2]>2) {
          paste("ERROR: Your input data is multi-factor experiment, please restart the APP and use 'Multi-factor Experiment' tab to input your data.")
        } else {
          paste("This is a single-factor experiment with factor - '", as.character(colnames(meta)[2]), "', the levels of this factor are:",sep="")
        }  
      })
  })
  
  
  output$GroupLevel <- renderText({ 
    if(input$dataSubmit)
      isolate({
        meta <- metaUpdate()$meta
        count <- metaUpdate()$count
        if (is.null(meta) & !is.null(count)) {stop("Please provide the corresponding input 2: Meta-data Table!")}
        else if (!is.null(meta) & is.null(count)) {stop("Please provide the input 1: Raw Count Data!")}
        else if (!is.null(meta) & !is.null(count)) {
          if (length(grep('Error',count[1]))==1) { 
            stop(paste(count[1]) )
          } else if (length(grep('Error',meta[1]))==1) {
            stop(paste(meta[1]) )
          } else {
            orgCount <- count
            orgMeta <- meta
            if (dim(orgCount)[2] != dim(orgMeta)[1] ) {stop("Input files do not correspond with each other. ")}
            #else if ( sum(colnames(orgCount) == orgMeta[,1]) !=dim(orgMeta)[1] ) {stop("Input files do not correspond with each other. ")}
            else if (dim(orgMeta)[2] < 2 ) {stop("Input 2 file format is wrong.") }
          }                 
        }
        
        if (dim(meta)[2]>2) {
          paste("ERROR!!!")
        } else {
          paste(as.character(levels(as.factor(meta[,2]))), collapse=", ") 
        }
      })
  })

  output$errorInputSingle <- renderText({
    if (input$dataSubmit )
      isolate({
        meta <- metaUpdate()$meta
        count <- metaUpdate()$count
        if (is.null(meta) & !is.null(count)) {stop("Please provide the corresponding input 2: Meta-data Table!")}
        else if (!is.null(meta) & is.null(count)) {stop("Please provide the input 1: Raw Count Data!")}
        else if (!is.null(meta) & !is.null(count)) {
          if (length(grep('Error',count[1]))==1) { 
            stop(paste(count[1]) )
          } else if (length(grep('Error',meta[1]))==1) {
            stop(paste(meta[1]) )
          } else {
            orgCount <- count
            orgMeta <- meta
            if (dim(orgCount)[2] != dim(orgMeta)[1] ) {stop("Input files do not correspond with each other. ")}
            #else if ( sum(colnames(orgCount) == orgMeta[,1]) !=dim(orgMeta)[1] ) {stop("Input files do not correspond with each other. ")}
            else if (dim(orgMeta)[2] < 2 ) {stop("Input 2 file format is wrong.") }
          }                 
        }
      })
    
  })
  #################################################
  
  #################################################
  ##Multi-factor Exp Res Summary
  output$overallDataSummaryMulti <- renderTable({ 
    if ( input$MultiSubmit )       
      isolate({
        #print("=========")
        #print(dataObs$orgMeta)
        #print("=========")
        meta <- metaUpdateMulti()$meta
        count <- metaUpdateMulti()$count
        if (is.null(meta) & !is.null(count)) {stop("Please provide the corresponding input 2: Meta-data Table!")}
        else if (!is.null(meta) & is.null(count)) {stop("Please provide the input 1: Raw Count Data!")}
        else if (!is.null(meta) & !is.null(count)) {
          if (length(grep('Error',count[1]))==1) { 
            stop(paste(count[1]) )
          } else if (length(grep('Error',meta[1]))==1) {
            stop(paste(meta[1]) )
          } else {
            orgCount <- count
            orgMeta <- meta
            if (dim(orgCount)[2] != dim(orgMeta)[1] ) {stop("Input files do not correspond with each other. ")}
            #else if ( sum(colnames(orgCount) == orgMeta[,1]) !=dim(orgMeta)[1] ) {stop("Input files do not correspond with each other. ")}
            else if (dim(orgMeta)[2] < 3 ) {stop("Input 2 file format is wrong.") }
          }          
        }
        no.samples <- length(colnames(count))
        no.gene <- dim(count)[1]
        observeEvent(input$MultiSubmit, { 
          progress$time$set(value = 1, detail = "processing 100%")
        })
        res.summary <- rbind(no.samples, no.gene)
        rownames(res.summary) <- c("Samples", "Tags")
        colnames(res.summary) <- "Number"
        res.summary
      })
    
  },digits=0, align="l|c")
  
  output$sampleGroupMulti <- renderTable({ 
    if ( input$MultiSubmit )     
      isolate({
        meta <- metaUpdateMulti()$meta
        count <- metaUpdateMulti()$count
        if (is.null(meta) & !is.null(count)) {stop("Please provide the corresponding input 2: Meta-data Table!")}
        else if (!is.null(meta) & is.null(count)) {stop("Please provide the input 1: Raw Count Data!")}
        else if (!is.null(meta) & !is.null(count)) {
          if (length(grep('Error',count[1]))==1) { 
            stop(paste(count[1]) )
          } else if (length(grep('Error',meta[1]))==1) {
            stop(paste(meta[1]) )
          } else {
            orgCount <- count
            orgMeta <- meta
            if (dim(orgCount)[2] != dim(orgMeta)[1] ) {stop("Input files do not correspond with each other. ")}
            else if ( sum(colnames(orgCount) == orgMeta[,1]) !=dim(orgMeta)[1] ) {stop("Input files do not correspond with each other. ")}
            else if (dim(orgMeta)[2] < 3 ) {stop("Input 2 file format is wrong.") }
          }          
        }
        
        groupinfo <- as.matrix(summary((datareactive())$samples$group))
        colnames(groupinfo) <- "No. in each group"
        groupinfo
      })
    
    
  },digits=0, align="l|c")
  
  
  
  output$sampleInfoMulti <- renderText({ 
    if ( input$MultiSubmit ) 
      isolate({
        #print("MMMMMMMM")
        #print(metaUpdateMulti()$meta)
        #print(head(metaUpdateMulti()$count))

        meta <- metaUpdateMulti()$meta
        count <- metaUpdateMulti()$count
        if (is.null(meta) & !is.null(count)) {stop("Please provide the corresponding input 2: Meta-data Table!")}
        else if (!is.null(meta) & is.null(count)) {stop("Please provide the input 1: Raw Count Data!")}
        else if (!is.null(meta) & !is.null(count)) {
          if (length(grep('Error',count[1]))==1) { 
            stop(paste(count[1]) )
          } else if (length(grep('Error',meta[1]))==1) {
            stop(paste(meta[1]) )
          } else {
            orgCount <- count
            orgMeta <- meta
            if (dim(orgCount)[2] != dim(orgMeta)[1] ) {stop("Input files do not correspond with each other. ")}
            else if ( sum(colnames(orgCount) == orgMeta[,1]) !=dim(orgMeta)[1] ) {stop("Input files do not correspond with each other. ")}
            else if (dim(orgMeta)[2] < 3 ) {stop("Input 2 file format is wrong.") }
          }        
          
        }
        
        paste(as.character(rownames((datareactive())$samples)), collapse=", " )
      })
    
  })
  
  output$sampleTitleMulti <- renderText({ 
    if ( input$MultiSubmit ) 
      
      isolate({
        meta <- metaUpdateMulti()$meta
        count <- metaUpdateMulti()$count
        if (is.null(meta) & !is.null(count)) {stop("Please provide the corresponding input 2: Meta-data Table!")}
        else if (!is.null(meta) & is.null(count)) {stop("Please provide the input 1: Raw Count Data!")}
        else if (!is.null(meta) & !is.null(count)) {
          if (length(grep('Error',count[1]))==1) { 
            stop(paste(count[1]) )
          } else if (length(grep('Error',meta[1]))==1) {
            stop(paste(meta[1]) )
          } else {
            orgCount <- count
            orgMeta <- meta
            if (dim(orgCount)[2] != dim(orgMeta)[1] ) {stop("Input files do not correspond with each other. ")}
            else if ( sum(colnames(orgCount) == orgMeta[,1]) !=dim(orgMeta)[1] ) {stop("Input files do not correspond with each other. ")}
            else if (dim(orgMeta)[2] < 3 ) {stop("Input 2 file format is wrong.") }
          }        
          
        }
        
        no.samples <- length(colnames(count))
        paste("A total of ", as.character(no.samples), " samples in the experiment, they are:", sep="")
      }) 
    
  })
  
  output$expDesignMulti <- renderText({
    if (input$MultiSubmit) 
      isolate({        
        meta <- metaUpdateMulti()$meta
        count <- metaUpdateMulti()$count
        if (is.null(meta) & !is.null(count)) {stop("Please provide the corresponding input 2: Meta-data Table!")}
        else if (!is.null(meta) & is.null(count)) {stop("Please provide the input 1: Raw Count Data!")}
        else if (!is.null(meta) & !is.null(count)) {
          if (length(grep('Error',count[1]))==1) { 
            stop(paste(count[1]) )
          } else if (length(grep('Error',meta[1]))==1) {
            stop(paste(meta[1]) )
          } else {
            orgCount <- count
            orgMeta <- meta
            if (dim(orgCount)[2] != dim(orgMeta)[1] ) {stop("Input files do not correspond with each other. ")}
            else if ( sum(colnames(orgCount) == orgMeta[,1]) !=dim(orgMeta)[1] ) {stop("Input files do not correspond with each other. ")}
            else if (dim(orgMeta)[2] < 3 ) {stop("Input 2 file format is wrong.") }
          }          
        }
        
        if (dim(meta)[2]>2) {
          paste("This is a multi-factor experiment with ", (dim(meta)[2]-1), " factors levels, they are", sep="" )
        } else {
          paste("ERROR: Your data input is a single-factor experiment, please restart the App and use 'Single-facotr Experiment' tab to input your data.")
        }       
      })
    
    
  })
  
  output$GroupLevelMulti <- renderText({ 
    if (input$MultiSubmit)      
      isolate({ 
        #org.counts <- dataObs$orgCount
        #print(head(org.counts))
        #metadata <- dataObs$orgMeta
        meta <- metaUpdateMulti()$meta
        count <- metaUpdateMulti()$count
        if (is.null(meta) & !is.null(count)) {stop("Please provide the corresponding input 2: Meta-data Table!")}
        else if (!is.null(meta) & is.null(count)) {stop("Please provide the input 1: Raw Count Data!")}
        else if (!is.null(meta) & !is.null(count)) {
          if (length(grep('Error',count[1]))==1) { 
            stop(paste(count[1]) )
          } else if (length(grep('Error',meta[1]))==1) {
            stop(paste(meta[1]) )
          } else {
            orgCount <- count
            orgMeta <- meta
            if (dim(orgCount)[2] != dim(orgMeta)[1] ) {stop("Input files do not correspond with each other. ")}
            else if ( sum(colnames(orgCount) == orgMeta[,1]) !=dim(orgMeta)[1] ) {stop("Input files do not correspond with each other. ")}
            else if (dim(orgMeta)[2] < 3 ) {stop("Input 2 file format is wrong.") }
          }          
        }
        
        if (dim(meta)[2]>2) {
          factor <- paste(as.character(colnames(meta)[-1]), collapse=", ")        
          group <- paste(as.character(levels((datareactive())$samples$group)), collapse=", ") 
          
          factor.group <- paste("\nFactor - ", as.character(colnames(meta)[2]), " includes factore levels of " , paste(as.character(levels(as.factor(meta[,2]))), collapse=", "), ".",sep="")          
          for (i in 3:length(meta[1,])) factor.group <- paste(factor.group, "\n", paste("Factor - ", as.character(colnames(meta)[i]), " include factor levels of " ,paste(as.character(levels(as.factor(meta[,i]))), collapse=", "), ".", sep=""), sep="")
          
          paste(as.character(factor), ".", as.character(factor.group), "\n\nThe combined factor levels are ", as.character(group), ".", sep="")
        } else {
          paste("ERROR!") 
        }
        
      })
    
    
  })
  
  output$errorInputMulti <- renderText({
    if ( input$MultiSubmit )
      isolate({
        meta <- metaUpdateMulti()$meta
        count <- metaUpdateMulti()$count
        if (is.null(meta) & !is.null(count)) {stop("Please provide the corresponding input 2: Meta-data Table!")}
        else if (!is.null(meta) & is.null(count)) {stop("Please provide the input 1: Raw Count Data!")}
        else if (!is.null(meta) & !is.null(count)) {
          if (length(grep('Error',count[1]))==1) { 
            stop(paste(count[1]) )
          } else if (length(grep('Error',meta[1]))==1) {
            stop(paste(meta[1]) )
          } else {
            orgCount <- count
            orgMeta <- meta
            if (dim(orgCount)[2] != dim(orgMeta)[1] ) {stop("Input files do not correspond with each other. ")}
            else if ( sum(colnames(orgCount) == orgMeta[,1]) !=dim(orgMeta)[1] ) {stop("Input files do not correspond with each other. ")}
            else if (dim(orgMeta)[2] < 3 ) {stop("Input 2 file format is wrong.") }
          }          
        }        
      })  
  })
  ##End data input tab panel
  ################################################
  
  ################################################
  ###tab panel Data summary
  output$orgLibsizeNormfactor <- renderTable({
    if(input$rmlow)
      isolate({        
        tab <- datareactive()$samples[,-1]
        colnames(tab) <- c("Library sizes", "Normalization factors")
        tab
      })
  }, align="l|cc", digits=c(0,0,2), display=c("s", "e", "f"))
  
  output$rmlowLibsizeNormfactor <- renderTable({ 
    if ( input$rmlow )
      isolate({ 
        tab <- (rmlowReactive())$samples[, -1]
        colnames(tab) <- c("Library sizes", "Normalization factors")
        tab
      })
  }, align="l|cc", digits=c(0,0,2), display=c("s", "e", "f"))
  
  output$orgSamplesize <- renderTable({
    if(input$rmlow)
      isolate({   
        no.samples <- length(colnames(datareactive()$counts))
        no.gene <- dim((datareactive())$counts)[1]
        res.summary <- rbind(no.samples, no.gene)
        rownames(res.summary) <- c("Samples", "Tags")
        colnames(res.summary) <- "Number"
        res.summary
      })
  },digits=0, align="l|c")
  
  output$rmlowSamplesize <- renderTable({ 
    if(input$rmlow) 
      isolate({ 
        if (as.numeric(input$gThreshold) > length(colnames(datareactive()$counts)) ) 
          stop("Cutoff sample number exceeds the total number of samples")
        
        no.samples <- length(colnames(rmlowReactive()$counts))
        no.gene <- dim(rmlowReactive()$counts)[1]
        res.summary <- rbind(no.samples, no.gene)
        rownames(res.summary) <- c("Samples", "Tags")
        colnames(res.summary) <- "Number"
        res.summary
      })
  },digits=0, align="l|c")
  
  output$sampleBoxplot <- renderPlot({ 
    if(input$rmlow) 
      isolate({
        if (as.numeric(input$gThreshold) > length(colnames(datareactive()$counts)) ) 
          stop("Cutoff sample number exceeds the total number of samples")
        
        Group <- as.factor(rmlowReactive()$samples$group)
        bx.p<-boxplot(cpm(rmlowReactive(), log=T)[,])
        bxp(bx.p, boxfill=as.numeric(Group)+1, 
            cex.axis=1.5, whisklwd=3, outcol=as.numeric(Group)+1, 
            main="Normalized sample distribution", cex.main=2)
      }) 
  })
  
  output$sampleMDS <- renderPlot({ 
    if(input$rmlow) 
      isolate({
        if (as.numeric(input$gThreshold) > length(colnames(datareactive()$counts)) ) 
          stop("Cutoff sample number exceeds the total number of samples")
        
        Group <- as.factor(rmlowReactive()$samples$group)
        observeEvent(input$rmlow, { 
          progress$time$set(value = 1, detail = "processing 100%")
        })
        par(mar=c(5,5,4,2))
        plotMDS(rmlowReactive(), col=as.numeric(Group)+1, 
                cex=2, main="MDS plot", ndim=3, gene.selection="common",
                xlab = "logFC dim 1", ylab="logFC dim 2", cex.lab=2, cex.main=2, cex.axis=1.5)
      })
  }) 
  
  
  output$errorFiltering <- renderText({
    input$rmlow 
    isolate({
      if (as.numeric(input$gThreshold) > length(colnames(datareactive()$counts)) ) {
        stop("Cutoff sample number exceeds the total number of samples")
        #paste("ERROR: cutoff sample number exceeds the total number of samples!")
      }
    })
  })  
  ##End Data summary tab Panel
  ################################################
  
  ################################################
  ###tab panel edgeR
  output$edgerGroupLevel <- renderText({ 
    paste("The available group levels are: " , paste(as.character(levels((datareactive())$samples$group)), collapse=", "), ".", sep="")
  })
  
  ##Estimation dispersion BCV
  output$edgerBCV <- renderPlot({ 
    if (input$edgerdeAnalysis)
      isolate({ 
        
        par(mar=c(5,5,2,2))
        plotBCV(edgerDispersionEst(), cex=0.5, cex.lab=1.8, cex.axis=1.5)
      })
  })
  
  output$edgerCommonDisp <- renderText({
    if(input$edgerdeAnalysis)
      isolate({
        paste("Esitmated biological coefficient of variation (BCV) is ", round(sqrt(edgerDispersionEst()$common.dispersion)*100, 4), "%", sep="")
      })
    
  })
  
  output$edgerTagwiseDispExp <- renderText({
    paste("Esitmated tagwise dispersion can be summarized as below:")
  })
  
  output$edgerTagwiseDisp <- renderTable({ 
    
    if(input$edgerdeAnalysis)
      isolate({       
        res<- as.table(t(summary(edgerDispersionEst()$tagwise.dispersion)))
        res <- t(res)
        colnames(res) <- "Tagwise"
        res    
      })
  }, digits=3, align = "l|c")
  ##############
  ##tab panel for DE analysis results   
  output$erroredgeR <- renderText({
    if(input$edgerdeAnalysis)
      isolate({
        Group <- edgerDispersionEst()$samples$group
        if (as.character(trim(input$edgercompGroup2)) =="" | as.character(trim(input$edgercompGroup1)) == "") {
          stop("Please select 2 group levels for DE analysis!")
        } else {
          if (as.character(trim(input$edgercompGroup2)) == as.character(trim(input$edgercompGroup1)) ) {
            stop("Group level 1 and level 2 are the same, please select 2 different group levels for DE analysis!")
          } else if( (! as.character(trim(input$edgercompGroup2)) %in% levels(Group)) | (! as.character(trim(input$edgercompGroup1)) %in% levels(Group)) ) {
            stop("Group level 1 or level 2 are not in the available group levels!")
          }
        }        
      })
  })
  
  output$edgerRes <- DT::renderDataTable({ 
    if(input$edgerdeAnalysis)
      isolate({ 
        tp <- topTags(edgerDEres(), n=Inf, adjust.method="BH", sort.by="PValue")
        res <- tp$table[,c("logFC", "PValue", "FDR")]
        #res$FC <- 2^res$logFC
        res$logFC <- round(res$logFC, digits = 3)
        res
      })
  }, 
  options = list(order = list(2, 'asc'), pageLength = 20, 
                 columnDefs = list(list(className = 'dt-center', targets = c(1,2,3)))),  
  colnames = c("log2FC"=2, "p"=3, "FDR"=4, "Tag/Gene Name"=1)
  )
  
  output$edgerTestDGEtitle <- renderText({
    if(input$edgerdeAnalysis)
      isolate({
        Group <- edgerDispersionEst()$samples$group
        if (as.character(trim(input$edgercompGroup2)) =="" | as.character(trim(input$edgercompGroup1)) == "") {
          stop("Please select 2 group levels for DE analysis!")
        } else {
          if (as.character(trim(input$edgercompGroup2)) == as.character(trim(input$edgercompGroup1)) ) {
            stop("Group level 1 and level 2 are the same, please select 2 different group levels for DE analysis!")
          } else if( (! as.character(trim(input$edgercompGroup2)) %in% levels(Group)) | (! as.character(trim(input$edgercompGroup1)) %in% levels(Group)) ) {
            stop("Group level 1 or level 2 are not in the available group levels!")
          }
        } 
        paste(paste(as.character(trim(input$edgercompGroup2)), as.character(trim(input$edgercompGroup1)), sep="-"), " DE analysis", sep="")
      })   
  })
  
  output$edgerTestDGE <- renderTable({ 
    if(input$edgerdeAnalysis)
      isolate({ 
        Group <- edgerDispersionEst()$samples$group
        if (as.character(trim(input$edgercompGroup2)) =="" | as.character(trim(input$edgercompGroup1)) == "") {
          stop("Please select 2 group levels for DE analysis!")
        } else {
          if (as.character(trim(input$edgercompGroup2)) == as.character(trim(input$edgercompGroup1)) ) {
            stop("Group level 1 and level 2 are the same, please select 2 different group levels for DE analysis!")
          } else if( (! as.character(trim(input$edgercompGroup2)) %in% levels(Group)) | (! as.character(trim(input$edgercompGroup1)) %in% levels(Group)) ) {
            stop("Group level 1 or level 2 are not in the available group levels!")
          }
        }
        
        edgerDGEsummary <- as.data.frame(summary((edgerDEfilter()$table)$filter))      
        #colnames(edgerDGEsummary) <- paste("No. of gene (",paste(as.character(input$edgercompGroup2), as.character(input$edgercompGroup1), sep="-"), ")", sep="")
        colnames(edgerDGEsummary) <- "Number"
        if (!is.null(rownames(edgerDGEsummary)==-1)) rownames(edgerDGEsummary)[rownames(edgerDGEsummary)==-1] <- "down-regulated DEG"
        if (!is.null(rownames(edgerDGEsummary)==1)) rownames(edgerDGEsummary)[rownames(edgerDGEsummary)==1] <- "up-regulated DEG"
        if (!is.null(rownames(edgerDGEsummary)==0)) rownames(edgerDGEsummary)[rownames(edgerDGEsummary)==0] <- "Non DEG"
        
        edgerDGEsummary
      })
  }, align="l|c", digits=0)
  
  output$edgerVolcano <- renderPlot({
    if(input$edgerdeAnalysis)
      isolate({
        fcval <- as.numeric(input$edgerfc)
        fdr <- as.numeric(input$edgerfdr)
        plotres <- edgerDEfilter()$table
        levels(plotres$filter) <- c("Down-regulated", "Non-DE", "Up-regulated")
        par(mar=c(4,4,2,2))
        if (input$edgerP == 'normp') {
          g <- ggplot(data=plotres, aes(x=logFC, y=-log10(PValue), colour=filter, shape=filter))
          g <- g + labs(x="log2 FC", y="-log10 (nominal p)")
        } else if (input$edgerP == 'fdrp') {
          g <- ggplot(data=plotres, aes(x=logFC, y=-log10(FDR), colour=filter, shape=filter))
          g <- g + labs(x="log2 FC", y="-log10 (FDR adjusted-p)")
        }
        g <- g + geom_point(size=4) + geom_hline(yintercept=-log10(as.numeric(input$edgerfdr))) 
        g <- g + geom_vline(xintercept=log2(fcval)) + geom_vline(xintercept=-log2(fcval))
        g <- g + theme(legend.title=element_blank(),legend.position = "top", legend.direction="vertical")
        g <- g + theme(axis.title.x=element_text(size=rel(1.5))) + theme(axis.title.y=element_text(size=rel(1.5),vjust=1.5))
        g <- g + theme(axis.text.x=element_text(size=rel(1.5))) + theme(axis.text.y=element_text(size=rel(1.5)))
        g <- g + theme(legend.text=element_text(size=rel(1.2)))
        g
      })
  })
  
  output$edgerDownload <- downloadHandler(  
    filename = function() {paste("edger-DEG-res-", Sys.Date(), ".txt", sep="")},
    content = function(file) {write.table(rbind(subset(edgerDEfilter()$table, filter==1),
                                                subset(edgerDEfilter()$table, filter==-1)), 
                                          file, row.names=T, col.names=NA, quote=F, sep="\t"
    )},
    contentType = "text"
  )
  ###end tab panel edgeR
  ################################################
  
  ################################################
  ###Limma-voom panel DE analysis
  output$errorVoom <- renderText({
    if(input$voomdeAnalysis)
      isolate({
        Group <- datareactive()$samples$group
        if (as.character(trim(input$voomcompGroup2)) =="" | as.character(trim(input$voomcompGroup1)) == "") {
          stop("Please select 2 groups levels for DE analysis!")
        } else if (as.character(trim(input$voomcompGroup1)) == as.character(trim(input$voomcompGroup2)) ) {
          stop("Group level 1 and level 2 are the same, please select 2 different group levels for DE analysis!")
        } else if( (! as.character(trim(input$voomcompGroup2)) %in% levels(Group)) | (! as.character(trim(input$voomcompGroup1)) %in% levels(Group)) ) {
          stop("Group level 1 or level 2 are not in the available group levels!")
        } 
      })    
  })  
  
  output$voomGroupLevel <- renderText({ 
    paste("The available group levels are: " ,paste(as.character(levels((datareactive())$samples$group)), collapse=", "), sep="")     
  })
  
  output$voomBCV <- renderPlot({
    if(input$voomdeAnalysis)
      isolate({
        Group <- rmlowReactive()$samples$group
        f <- factor(Group, levels=levels(Group))
        design <- model.matrix(~0+f)
        rownames(design) <- rownames(rmlowReactive()$samples)
        colnames(design) <- levels(Group)
        par(mar=c(4,4,2,2))
        voom(rmlowReactive(), design=design, plot=T, normalize="quantile")
      })
  })
  
  output$voomRes <- DT::renderDataTable({
    if (input$voomdeAnalysis)
      isolate({ 
        tpvoom <- topTable(voomDEres(), n=Inf, adjust.method="BH", sort.by="p")
        res <- tpvoom[,c("logFC","P.Value", "adj.P.Val")]
        res$logFC <- round(res$logFC, digits=3)
        observeEvent(input$voomdeAnalysis, { 
          progress$time$set(value = 1, detail = "processing 100%")
        })
        res
      })
  }, 
  options = list(order = list(2, 'asc'), columnDefs = list(list(className = 'dt-center', targets = c(1,2,3))), pageLength = 20),  
  colnames = c("log2FC"=2, "p"=3, "FDR"=4, "Tag/Gene Name"=1)
  )
  
  output$voomTestDGEtitle <- renderText({
    if(input$voomdeAnalysis)
      isolate({
        Group <- datareactive()$samples$group
        if (as.character(trim(input$voomcompGroup2)) =="" | as.character(trim(input$voomcompGroup1)) == "") {
          stop("Please select 2 groups levels for DE analysis!")
        } else if (as.character(trim(input$voomcompGroup1)) == as.character(trim(input$voomcompGroup2)) ) {
          stop("Group level 1 and level 2 are the same, please select 2 different group levels for DE analysis!")
        } else if( (! as.character(trim(input$voomcompGroup2)) %in% levels(Group)) | (! as.character(trim(input$voomcompGroup1)) %in% levels(Group)) ) {
          stop("Group level 1 or level 2 are not in the available group levels!")
        } 
        paste(paste(as.character(trim(input$voomcompGroup2)), as.character(trim(input$voomcompGroup1)), sep="-"), " DE analysis", sep="")
      })   
  })
  
  output$voomTestDGE <- renderTable({
    if(input$voomdeAnalysis)
      isolate({
        Group <- datareactive()$samples$group
        if (as.character(trim(input$voomcompGroup2)) =="" | as.character(trim(input$voomcompGroup1)) == "") {
          stop("Please select 2 groups levels for DE analysis!")
        } else if (as.character(trim(input$voomcompGroup1)) == as.character(trim(input$voomcompGroup2)) ) {
          stop("Group level 1 and level 2 are the same, please select 2 different group levels for DE analysis!")
        } else if( (! as.character(trim(input$voomcompGroup2)) %in% levels(Group)) | (! as.character(trim(input$voomcompGroup1)) %in% levels(Group)) ) {
          stop("Group level 1 or level 2 are not in the available group levels!")
        } 
        voomDGEsummary <- as.data.frame(summary(voomDEfilter()$filter))
        colnames(voomDGEsummary) <- "Number"
        #colnames(voomDGEsummary) <- paste("No. of gene (",paste(as.character(input$voomcompGroup2), as.character(input$voomcompGroup1), sep="-"), ")", sep="")
        #rownames(voomDGEsummary) <- c("down-regulated DEG", "Non DEG", "up-regulated DEG")
        if (!is.null(rownames(voomDGEsummary)==-1)) rownames(voomDGEsummary)[rownames(voomDGEsummary)==-1] <- "down-regulated DEG"
        if (!is.null(rownames(voomDGEsummary)==1)) rownames(voomDGEsummary)[rownames(voomDGEsummary)==1] <- "up-regulated DEG"
        if (!is.null(rownames(voomDGEsummary)==0)) rownames(voomDGEsummary)[rownames(voomDGEsummary)==0] <- "Non DEG"
        
        voomDGEsummary
      })
  }, align="l|c", digits=0)
  
  output$voomVolcano <- renderPlot({
    if (input$voomdeAnalysis)
      isolate({
        plotres <- voomDEfilter()
        levels(plotres$filter) <- c("Down-regulated", "Non-DE", "Up-regulated")
        if (input$voomP == 'normp') {
          g <- ggplot(data=plotres, aes(x=logFC, y=-log10(P.Value), colour=filter, shape=filter))
          g <- g + labs(x="log2 FC", y="-log10 (nominal p)")
        } else if (input$voomP == 'fdrp') {
          g <- ggplot(data=plotres, aes(x=logFC, y=-log10(adj.P.Val), colour=filter, shape=filter))
          g <- g + labs(x="log2 FC", y="-log10 (FDR adjusted-p)")
        }
        g <- g + geom_point(size=4) + geom_hline(yintercept=-log10(as.numeric(input$voomfdr))) 
        g <- g + geom_vline(xintercept=log2(as.numeric(input$voomfc))) + geom_vline(xintercept=-log2(as.numeric(input$voomfc)))
        g <- g + theme(legend.title=element_blank(),legend.position = "top", legend.direction="vertical")
        g <- g + theme(axis.title.x=element_text(size=rel(1.5))) + theme(axis.title.y=element_text(size=rel(1.5),vjust=1.5))
        g <- g + theme(axis.text.x=element_text(size=rel(1.5))) + theme(axis.text.y=element_text(size=rel(1.5)))
        g <- g + theme(legend.text=element_text(size=rel(1.2)))
        g
      })
  })
  
  output$voomDownload <- downloadHandler(  
    filename = function() {paste("voom-DEG-res-", Sys.Date(), ".txt", sep="")},
    content = function(file) {write.table(rbind(subset(voomDEfilter(), filter==1),
                                                subset(voomDEfilter(), filter==-1)), 
                                          file, row.names=T, col.names=NA, quote=F, sep="\t"
    )},
    contentType = "text"
  )
  
  ###end tab panel limma-voom
  ################################################
  
  ################################################
  ###DEseq2 panel DE analysis
  output$errorDeseq2 <- renderText({
    if(input$deseq2deAnalysis)
      isolate({
        Group <- datareactive()$samples$group
        if (as.character(trim(input$deseq2compGroup2)) =="" | as.character(trim(input$deseq2compGroup1)) == "") {
          stop("Please select 2 groups levels for DE analysis!")
        } else if (as.character(trim(input$deseq2compGroup1)) == as.character(trim(input$deseq2compGroup2)) ) {
          stop("Group level 1 and level 2 are the same, please select 2 different group levels for DE analysis!")
        } else if( (! as.character(trim(input$deseq2compGroup2)) %in% levels(Group)) | (! as.character(trim(input$deseq2compGroup1)) %in% levels(Group)) ) {
          stop("Group level 1 or level 2 are not in the available group levels!")
        }        
      })    
  }) 
  
  output$deseq2GroupLevel <- renderText({  
    paste("The available group levels are: " ,paste(as.character(levels((datareactive())$samples$group)), collapse=", "), sep="")  
  })
  
  output$deseq2BCV <- renderPlot({
    if (input$deseq2deAnalysis)
      isolate({
        observeEvent(input$deseq2deAnalysis, { 
          progress$time$set(value = 0.3, detail = "processing 30%")
        })
        par(mar=c(4,5,2,2))
        plotDispEsts(deseq2Res(), cex.lab=1.8, cex.axis=1.5)
      })
  })
  
  output$deseq2Res <- DT::renderDataTable({
    if (input$deseq2deAnalysis)
      isolate({
        res.pre <- deseq2DEres()[,c("log2FoldChange", "pvalue", "padj")]
        res <- as.data.frame(res.pre) 
        res$log2FoldChange <- round(res$log2FoldChange, digits=3)
        observeEvent(input$deseq2deAnalysis, { 
          progress$time$set(value = 1, detail = "processing 100%")
        })
        res
      })
  }, 
  options = list(order = list(2, 'asc'), columnDefs = list(list(className = 'dt-center', targets = c(1,2,3))), pageLength = 20),  
  colnames = c("log2FC"=2, "p"=3, "FDR"=4,"Tag/Gene Name"=1)
  )
  
  output$deseq2TestDGEtitle <- renderText({
    if (input$deseq2deAnalysis)
      isolate({
        Group <- datareactive()$samples$group
        if (as.character(trim(input$deseq2compGroup2)) =="" | as.character(trim(input$deseq2compGroup1)) == "") {
          stop("Please select 2 groups levels for DE analysis!")
        } else if (as.character(trim(input$deseq2compGroup1)) == as.character(trim(input$deseq2compGroup2)) ) {
          stop("Group level 1 and level 2 are the same, please select 2 different group levels for DE analysis!")
        } else if( (! as.character(trim(input$deseq2compGroup2)) %in% levels(Group)) | (! as.character(trim(input$deseq2compGroup1)) %in% levels(Group)) ) {
          stop("Group level 1 or level 2 are not in the available group levels!")
        } 
        paste(paste(as.character(trim(input$deseq2compGroup2)), as.character(trim(input$deseq2compGroup1)), sep="-"), " DE analysis", sep="")
      })
  })
  
  output$deseq2TestDGE <- renderTable({
    if (input$deseq2deAnalysis)
      isolate({
        Group <- datareactive()$samples$group
        if (as.character(trim(input$deseq2compGroup2)) =="" | as.character(trim(input$deseq2compGroup1)) == "") {
          stop("Please select 2 groups levels for DE analysis!")
        } else if (as.character(trim(input$deseq2compGroup1)) == as.character(trim(input$deseq2compGroup2)) ) {
          stop("Group level 1 and level 2 are the same, please select 2 different group levels for DE analysis!")
        } else if( (! as.character(trim(input$deseq2compGroup2)) %in% levels(Group)) | (! as.character(trim(input$deseq2compGroup1)) %in% levels(Group)) ) {
          stop("Group level 1 or level 2 are not in the available group levels!")
        } 
        deseq2DGEsummary <- as.data.frame(summary(deseq2DEfilter()$filter))
        #colnames(deseq2DGEsummary) <- paste("No. of gene (",paste(as.character(input$deseq2compGroup2), as.character(input$deseq2compGroup1), sep="-"), ")", sep="")
        colnames(deseq2DGEsummary) <- "Number"
        if (!is.null(rownames(deseq2DGEsummary)==-1)) rownames(deseq2DGEsummary)[rownames(deseq2DGEsummary)==-1] <- "down-regulated DEG"
        if (!is.null(rownames(deseq2DGEsummary)==1)) rownames(deseq2DGEsummary)[rownames(deseq2DGEsummary)==1] <- "up-regulated DEG"
        if (!is.null(rownames(deseq2DGEsummary)==0)) rownames(deseq2DGEsummary)[rownames(deseq2DGEsummary)==0] <- "Non DEG"
        
        deseq2DGEsummary
      })
  }, align="l|c", digits=0)
  
  output$deseq2Volcano <- renderPlot({
    if (input$deseq2deAnalysis)
      isolate({
        plotres <- deseq2DEfilter()
        levels(plotres$filter) <- c("Down-regulated", "Non-DE", "Up-regulated")
        if (input$deseq2P == 'normp') {
          g <- ggplot(data=plotres, aes(x=log2FoldChange, y=-log10(pvalue), colour=filter, shape=filter))
          g <- g + labs(x="log2 FC", y="-log10 (nominal p)")
        } else if (input$deseq2P == 'fdrp') {
          g <- ggplot(data=plotres, aes(x=log2FoldChange, y=-log10(padj), colour=filter, shape=filter))
          g <- g + labs(x="log2 FC", y="-log10 (FDR adjusted-p)")
        }           
        g <- g + geom_point(size=4) + geom_hline(yintercept=-log10(as.numeric(input$deseq2fdr))) 
        g <- g + geom_vline(xintercept=log2(as.numeric(input$deseq2fc))) + geom_vline(xintercept=-log2(as.numeric(input$deseq2fc)))
        g <- g + theme(legend.title=element_blank(),legend.position = "top", legend.direction="vertical")
        g <- g + theme(axis.title.x=element_text(size=rel(1.5))) + theme(axis.title.y=element_text(size=rel(1.5),vjust=1.5))
        g <- g + theme(axis.text.x=element_text(size=rel(1.5))) + theme(axis.text.y=element_text(size=rel(1.5)))
        g <- g + theme(legend.text=element_text(size=rel(1.2)))
        g
      })
  })
  
  output$deseq2Download <- downloadHandler(  
    filename = function() {paste("DESeq2-DEG-res-", Sys.Date(), ".txt", sep="")},
    content = function(file) {write.table(rbind(subset(deseq2DEfilter(), filter==1),
                                                subset(deseq2DEfilter(), filter==-1)), 
                                          file, row.names=T, col.names=NA, quote=F, sep="\t"
    )},
    contentType = "text"
  )
  ###End DEseq2 panel DE analysis
  ################################################
  
  ################################################
  ##DE comparison results
  
  output$errorComp <- renderText({
    if(input$decompAnalysis)
      isolate({
        Group <- datareactive()$samples$group
        if (as.character(trim(input$decompGroup2)) =="" | as.character(trim(input$decompGroup1)) == "") {
          stop("Please select 2 groups levels for DE analysis!")
        } else if (as.character(trim(input$decompGroup1)) == as.character(trim(input$decompGroup2)) ) {
          stop("Group level 1 and level 2 are the same, please select 2 different group levels for DE analysis!")
        } else if( (! as.character(trim(input$decompGroup2)) %in% levels(Group)) | (! as.character(trim(input$decompGroup1)) %in% levels(Group)) ) {
          stop("Group level 1 or level 2 are not in the available group levels!")
        } 
      })
    
  })
  
  output$compGroupLevel <- renderText({ 
    
    paste("The available group levels are: " ,paste(as.character(levels((datareactive())$samples$group)), collapse=", "), sep="")
    
  })
  
  output$decomp <- renderPlot({
    if (input$decompAnalysis)
      isolate({      
        Group <- datareactive()$samples$group
        if (as.character(trim(input$decompGroup2)) =="" | as.character(trim(input$decompGroup1)) == "") {
          stop("Please select 2 groups levels for DE analysis!")
        } else if (as.character(trim(input$decompGroup1)) == as.character(trim(input$decompGroup2)) ) {
          stop("Group level 1 and level 2 are the same, please select 2 different group levels for DE analysis!")
        } else if( (! as.character(trim(input$decompGroup2)) %in% levels(Group)) | (! as.character(trim(input$decompGroup1)) %in% levels(Group)) ) {
          stop("Group level 1 or level 2 are not in the available group levels!")
        } 
        if (as.character("edger")%in%input$decompMethods & as.character("voom")%in%input$decompMethods & as.character("deseq2")%in%input$decompMethods) {
          venncount <- vennCounts(as.matrix(cbind(cbind(edgerDecomp()$table$filter, voomDecomp()$filter),deseq2Decomp()$filter)))
          vennDiagram(venncount, names=c("edgeR","voom", "DESeq2"),
                      mar=c(0,0,0,0), circle.col=c("red", "green", "blue"), 
                      counts.col=c(1, "red"), lwd=3)
        } else if (as.character("edger")%in%input$decompMethods & as.character("voom")%in%input$decompMethods) {
          vennDiagram(cbind(edgerDecomp()$table$filter, voomDecomp()$filter), include=c("both"), 
                      names=c("edgeR", "voom"),mar=c(0,0,0,0), circle.col=c("red", "green"), 
                      counts.col=c(1, "red"), lwd=3)
        } else if (as.character("edger")%in%input$decompMethods & as.character("deseq2")%in%input$decompMethods) {
          vennDiagram(cbind(edgerDecomp()$table$filter, deseq2Decomp()$filter), include=c("both"), 
                      names=c("edgeR", "DESeq2"),mar=c(0,0,0,0), circle.col=c("red", "green"), 
                      counts.col=c(1, "red"), lwd=3)
        } else if (as.character("voom")%in%input$decompMethods & as.character("deseq2")%in%input$decompMethods) {
          vennDiagram(cbind(voomDecomp()$filter, deseq2Decomp()$filter), include=c("both"), 
                      names=c("voom", "DESeq2"),mar=c(0,0,0,0), circle.col=c("red", "green"), 
                      counts.col=c(1, "red"), lwd=3)
        }
      })
  })
  
  output$decompTitle <- renderText({
    if (input$decompAnalysis)
      isolate({
        Group <- datareactive()$samples$group
        if (as.character(trim(input$decompGroup2)) =="" | as.character(trim(input$decompGroup1)) == "") {
          stop("Please select 2 groups levels for DE analysis!")
        } else if (as.character(trim(input$decompGroup1)) == as.character(trim(input$decompGroup2)) ) {
          stop("Group level 1 and level 2 are the same, please select 2 different group levels for DE analysis!")
        } else if( (! as.character(trim(input$decompGroup2)) %in% levels(Group)) | (! as.character(trim(input$decompGroup1)) %in% levels(Group)) ) {
          stop("Group level 1 or level 2 are not in the available group levels!")
        } 
        
        if (input$decompP == 'normp') {
          title <- paste("Below results are based on the nominal p ", 
                         "with filtering level of norminal p = ", as.character(input$decompfdr),
                         " and FC = ", as.character(input$decompfc), sep="")
        }
        if (input$decompP == 'fdrp') {
          title <- paste("Below results are based on the FDR-adjusted p ", 
                         "with filtering level of FDR-adjusted p = ", as.character(input$decompfdr),
                         " and FC = ", as.character(input$decompfc), sep="")
        }
        title
      })
  })
  
  output$decompText <- renderText({
    if (input$decompAnalysis)
      isolate({
        res.matrix <- decompRes()
        if (as.character("edger")%in%input$decompMethods & as.character("voom")%in%input$decompMethods & as.character("deseq2")%in%input$decompMethods) {
          res.text.edgeR <- round(res.matrix[7,1]/res.matrix[1,1]*100, digits=2)
          res.text.voom <- round(res.matrix[7,1]/res.matrix[2,1]*100, digits=2)
          res.text.deseq2 <- round(res.matrix[7,1]/res.matrix[3,1]*100, digits=2)
          res.text <- paste(as.character(res.text.edgeR),"% identified DEGs with edgeR were identified by all three methods.\n",
                            as.character(res.text.voom),"% identified DEGs with limma-voom were identified by all three methods.\n", 
                            as.character(res.text.deseq2),"% identified DEG withe DESeq2 were identified by all three methods.\n",
                            sep="")     
        } else if (as.character("edger")%in%input$decompMethods & as.character("voom")%in%input$decompMethods) {
          res.text.edgeR <- round(res.matrix[4,1]/res.matrix[1,1]*100, digits=2)
          res.text.voom <- round(res.matrix[4,1]/res.matrix[2,1]*100, digits=2)
          res.text <- paste(as.character(res.text.edgeR),"% identified DEGs with edgeR were identified by both edgeR and limma-voom.\n",
                            as.character(res.text.voom),"% identified DEGs with limma-voom were identified by both edgeR and limma-voom.\n",
                            sep="")        
        } else if (as.character("edger")%in%input$decompMethods & as.character("deseq2")%in%input$decompMethods) {
          res.text.edgeR <- round(res.matrix[5,1]/res.matrix[1,1]*100, digits=2)
          res.text.deseq2 <- round(res.matrix[5,1]/res.matrix[3,1]*100, digits=2)
          res.text <- paste(as.character(res.text.edgeR),"% identified DEGs with edgeR were identified by both edgeR and DESeq2.\n",
                            as.character(res.text.deseq2),"% identified DEGs with DESeq2 were identified by both edgeR and DESeq2.\n",
                            sep="")  
        } else if (as.character("voom")%in%input$decompMethods & as.character("deseq2")%in%input$decompMethods) {
          res.text.voom <- round(res.matrix[6,1]/res.matrix[2,1]*100, digits=2)
          res.text.deseq2 <- round(res.matrix[6,1]/res.matrix[3,1]*100, digits=2)
          res.text <- paste(as.character(res.text.voom),"% identified DEGs with limma-voom were identified by both limma-voom and DESeq2.\n",
                            as.character(res.text.deseq2),"% identified DEGs with DESeq2 were identified by both limma-voom and DESeq2.\n",
                            sep="")  
        }
        res.text
      })
  })
  
  output$decompTab <- renderTable({
    if (input$decompAnalysis)
      isolate({
        res.matrix <- decompRes()
        if (as.character("edger")%in%input$decompMethods & as.character("voom")%in%input$decompMethods & as.character("deseq2")%in%input$decompMethods) {
          res <- res.matrix
        } else if (as.character("edger")%in%input$decompMethods & as.character("voom")%in%input$decompMethods) {
          res <- data.frame(x <- res.matrix[c(1,2,4),])
          rownames(res) <- rownames(res.matrix)[c(1,2,4)]   
        } else if (as.character("edger")%in%input$decompMethods & as.character("deseq2")%in%input$decompMethods) {
          res <- data.frame(x <- res.matrix[c(1,3,5),])
          rownames(res) <- rownames(res.matrix)[c(1,3,5)]
        } else if (as.character("voom")%in%input$decompMethods & as.character("deseq2")%in%input$decompMethods) {
          res <- data.frame(x <- res.matrix[c(2,3,6),])
          rownames(res) <- rownames(res.matrix)[c(2,3,6)]
        }
        colnames(res) <- "No. identified DEGs"
        observeEvent(input$decompAnalysis, { 
          progress$time$set(message = "Comparison analysis", value = 0.7)
          progress$time$set(value = 1, detail = "processing 100%")
        })
        res
      })
  },digits = 0, align="l|c")
  
}) 


