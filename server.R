## This App is used for RNA-seq DE analysis with different methods together with method comparison
## Developed by Yan Li, last update on June, 2015

shinyServer(function(input, output, session) {
  
  ##Reactive expression object for original row count
  datareactive <- reactive ({  
    #inputDataFile <- input$countFile       
    if ( !is.null(input$countFile) & is.null(input$countFileMulti) ) {
      inputDataFile <- input$countFile 
      org.counts <- read.delim(inputDataFile$datapath, header=T, sep=input$coutFileSep, row.names=1 )
    } else if ( !is.null(input$countFileMulti) & is.null(input$countFile) ) {
      inputDataFile <- input$countFileMulti
      org.counts <- read.delim(inputDataFile$datapath, header=T, sep=input$coutFileSepMulti, row.names=1 )
    } else if (is.null(input$countFile) & is.null(input$countFileMulti) ) {
      org.counts <- read.delim(paste(getwd(),"data/TestData-feature-count-res.txt",sep="/"), header=T, row.names=1)
    }  
    
    #inputMetatab <- input$metaTab
    if ( !is.null(input$metaTab) & is.null(input$metaTabMulti) ) {
      inputMetatab <- input$metaTab
      metadata <- read.delim(inputMetatab$datapath, header=T, sep=input$metaSep)
    } else if ( !is.null(input$metaTabMulti) & is.null(input$metaTab) ) {
      inputMetatab <- input$metaTabMulti
      metadata <- read.delim(inputMetatab$datapath, header=T, sep=input$metaSepMulti)
    } else if ( is.null(input$metaTab) & is.null(input$metaTabMulti) ) {
      metadata <- read.delim(paste(getwd(),"/data/metatable.txt",sep=""), header=T)
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
    dge.count <- calcNormFactors(datareactive())
    cpm.count <- cpm(dge.count)
    keep = rowSums(cpm.count >= as.numeric(input$cpmVal) ) >= as.numeric(input$gThreshold)
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
    comp <- makeContrasts(contrasts=paste(as.character(input$edgercompGroup2), as.character(input$edgercompGroup1), sep="-"), levels=design )
    test <- glmLRT(edgerglmFit(), contrast=comp)
    test
  })
  
  #Reactive expression object for edgeR desideTestsDEG results
  edgerDEfilter <- reactive({
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
    Group <- datareactive()$samples$group
    f <- factor(Group, levels=levels(Group))
    design <- model.matrix(~0+f)
    rownames(design) <- rownames(rmlowReactive()$samples)
    colnames(design) <- levels(Group)
    par(mar=c(0,0,0,0))
    v <- voom(rmlowReactive(), design=design, plot=T, normalize="quantile")
    fit <- lmFit(v, design)
    fit
  })
  
  ##Reactive expression object for limma-voom DE analysis after contrast.fit
  voomDEres <- reactive({
    Group <- datareactive()$samples$group
    f <- factor(Group, levels=levels(Group))
    design <- model.matrix(~0+f)
    rownames(design) <- rownames(rmlowReactive()$samples)
    colnames(design) <- levels(Group)
    
    comp <- makeContrasts(contrasts=paste(as.character(input$edgercompGroup2), as.character(input$edgercompGroup1), sep="-"), levels=design )
    
    contrast.fit <- contrasts.fit(voomRes(), contrasts=comp)
    contrast.fit <- eBayes(contrast.fit)  
    contrast.fit
  })
  
  ##Reactive expression object for limma-voom decideTests results
  voomDEfilter <- reactive({
    voomfcval <- as.numeric(input$voomfc)
    voomfdr <- as.numeric(input$voomfdr)
    if (input$voomP == 'normp') {
      tpvoom <- topTable(voomDEres(), number=Inf, adjust="BH", sort.by="none") 
      filtervoom <- decideTests(voomDEres(), adjust.method="none", p.value=voomfdr, lfc=log2(voomfcval))
    } else if (input$voomP == 'fdrp') {
      tpvoom <- topTable(voomDEres(), number=Inf, adjust="BH", sort.by="none") 
      filtervoom <- decideTests(voomDEres(), adjust.method="BH", p.value=voomfdr, lfc=log2(voomfcval))
    }
    tpvoom$filter <- as.factor(filtervoom)
    tpvoom
  })
  
  ##Reactive expression object for DEseq2 DESeq analysis results
  deseq2Res <- reactive({
    Group <- datareactive()$samples$group
    colData <- data.frame(Group)
    dds <- DESeqDataSetFromMatrix(rmlowReactive()$count, colData=colData, design=formula(~Group) )
    dds <- DESeq(dds, test="Wald")  
    dds
  })
  
  ##Reactive expression object for DESeq2 results()
  deseq2DEres <- reactive({
    res <- results(deseq2Res(), contrast=c("Group", as.character(input$deseq2compGroup2), as.character(input$deseq2compGroup1)), format="DataFrame")
    as.data.frame(res)
  })
  
  ##Reactive expression object for limma-voom decideTests results
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
    
    comp <- makeContrasts(contrasts=paste(as.character(input$decompGroup2), as.character(input$decompGroup1), sep="-"), levels=design )
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
    
    comp <- makeContrasts(contrasts=paste(as.character(input$decompGroup2), as.character(input$decompGroup1), sep="-"), levels=design )
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
    deseq2res  <- results(deseq2Res(), contrast=c("Group", as.character(input$decompGroup2), as.character(input$decompGroup1)), format="DataFrame")
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
  
  output$countTabSamp1 <- renderTable({     
    org.counts <- read.delim(paste(getwd(),"www/dataInput1-exp.txt",sep="/"), header=T, row.names=1)
    org.counts
  }, align="l|cccccc")
  
  output$metaTabSamp <- renderTable({
    metadata <- read.delim(paste(getwd(),"/www/dataInput2-exp0.txt",sep=""), header=T)
    metadata
  }, include.rownames=F, align="llc")
  
  output$metaTabSamp1 <- renderTable({
    metadata <- read.delim(paste(getwd(),"/www/dataInput2-exp0.txt",sep=""), header=T)
    metadata
  }, align="ll|c", include.rownames=F, caption = "Single-factor",
  caption.placement = getOption("xtable.caption.placement", "top"), 
  caption.width = getOption("xtable.caption.width", NULL))
  
  output$countTabSampMulti1 <- renderTable({     
    org.counts <- read.delim(paste(getwd(),"www/Multi-dataInput1-exp.txt",sep="/"), header=T, row.names=1)
    org.counts
  }, align="l|cccccc")
  
  output$multimetaTabSamp22 <- renderTable({
    metadata <- read.delim(paste(getwd(),"/www/Multi-dataInput2-exp2.txt",sep=""), header=T)
    metadata
  }, include.rownames=F, align="llcccc")
  
  output$multimetaTabSamp1 <- renderTable({
    metadata <- read.delim(paste(getwd(),"/www/Multi-dataInput2-exp1.txt",sep=""), header=T)
    metadata
  }, align="ll|ccc", include.rownames=F, caption = "Multi-factor",
  caption.placement = getOption("xtable.caption.placement", "top"), 
  caption.width = getOption("xtable.caption.width", NULL))
  
  ################################################################################################
  ################################################################################################
  
  output$overallDataSummary <- renderTable({ 
    input$dataSubmit
    isolate({
      no.samples <- length(colnames(datareactive()$counts))
      no.gene <- dim((datareactive())$counts)[1]
      res.summary <- rbind(no.samples, no.gene)
      rownames(res.summary) <- c("Samples", "Tags")
      colnames(res.summary) <- "Number"
      res.summary
    })
    
  },digits=0, align="l|c")
  
  output$sampleGroup <- renderTable({ 
    input$dataSubmit
    isolate({
      groupinfo <- as.matrix(summary((datareactive())$samples$group))
      colnames(groupinfo) <- "No. in each group"
      groupinfo
    })
  },digits=0, align="l|c")
  
  output$sampleInfo <- renderText({ 
    input$dataSubmit
    isolate({
      paste(as.character(rownames((datareactive())$samples)), collapse=", " )
    })
  })
  
  output$sampleTitle <- renderText({ 
    input$dataSubmit
    isolate({
      no.samples <- length(colnames(datareactive()$counts))
      paste("A total of ", as.character(no.samples), " samples in the experiment, they are:", sep="")
    })  
  })
  
  output$expDesign <- renderText({
    input$dataSubmit
    isolate({
      inputMetatab <- input$metaTab
      
      if (is.null(inputMetatab)) metadata <- read.delim(paste(getwd(),"/data/metatable.txt",sep=""), header=T)
      else metadata <- read.delim(inputMetatab$datapath, header=T, sep=input$metaSep)
      
      if (dim(metadata)[2]>2) {
        paste("ERROR: Your input data is multi-factor experiment, please use 'Multi-factor Experiment' tab to input your data.")
      } else {
        paste("This is a single-factor experiment with factor - '", as.character(colnames(metadata)[2]), "', the level of this factor is",sep="")
      }  
    })
  })
  
  output$GroupLevel <- renderText({ 
    input$dataSubmit
    isolate({
      inputMetatab <- input$metaTab
      
      if (is.null(inputMetatab)) metadata <- read.delim(paste(getwd(),"/data/metatable.txt",sep=""), header=T)
      else metadata <- read.delim(inputMetatab$datapath, header=T, sep=input$metaSep)
      
      if (dim(metadata)[2]>2) {
        paste("ERROR!!!")
      } else {
        paste(as.character(levels(as.factor(metadata[,2]))), collapse=", ") 
      }
    })
  })
  
  #################################################
  ##Multi-factor Exp Res Summary
  output$overallDataSummaryMulti <- renderTable({ 
    if ( input$dataSubmitMulti ) {
      isolate({
        no.samples <- length(colnames(datareactive()$counts))
        no.gene <- dim((datareactive())$counts)[1]
        res.summary <- rbind(no.samples, no.gene)
        rownames(res.summary) <- c("Samples", "Tags")
        colnames(res.summary) <- "Number"
        res.summary
      })
    }    
  },digits=0, align="l|c")
  
  output$sampleGroupMulti <- renderTable({ 
    if ( input$dataSubmitMulti ) {
      isolate({
        groupinfo <- as.matrix(summary((datareactive())$samples$group))
        colnames(groupinfo) <- "No. in each group"
        groupinfo
      })
    }
    
  },digits=0, align="l|c")
  
  output$sampleInfoMulti <- renderText({ 
    if ( input$dataSubmitMulti ) {
      isolate({
        paste(as.character(rownames((datareactive())$samples)), collapse=", " )
      })
    }
  })
  
  output$sampleTitleMulti <- renderText({ 
    if ( input$dataSubmitMulti ) {
      isolate({
        no.samples <- length(colnames(datareactive()$counts))
        paste("A total of ", as.character(no.samples), " samples in the experiment, they are:", sep="")
      }) 
    }     
  })
  
  output$expDesignMulti <- renderText({
    if ( input$dataSubmitMulti ) {
      isolate({        
        inputMetatab <- input$metaTabMulti
        metadata <- read.delim(inputMetatab$datapath, header=T, sep=input$metaSepMulti)
        if (dim(metadata)[2]>2) {
          paste("This is a multi-factor experiment with ", (dim(metadata)[2]-1), " factors levels, they are", sep="" )
        } else {
          paste("ERROR: Your data input is a single-factor experiment, please use 'Single-facotr Experiment' tab to input your data.")
        }       
      })
    }
    
  })
  
  output$GroupLevelMulti <- renderText({ 
    if ( input$dataSubmitMulti ) {
      isolate({ 
        inputMetatab <- input$metaTabMulti
        metadata <- read.delim(inputMetatab$datapath, header=T, sep=input$metaSepMulti)
        
        if (dim(metadata)[2]>2) {
          factor <- paste(as.character(colnames(metadata)[-1]), collapse=", ")        
          group <- paste(as.character(levels((datareactive())$samples$group)), collapse=", ") 
          
          factor.group <- paste("\nFactor - ", as.character(colnames(metadata)[2]), " includes factore levels of " , paste(as.character(levels(as.factor(metadata[,2]))), collapse=", "), ".",sep="")          
          for (i in 3:length(metadata[1,])) factor.group <- paste(factor.group, "\n", paste("Factor - ", as.character(colnames(metadata)[i]), " include factor levels of " ,paste(as.character(levels(as.factor(metadata[,i]))), collapse=", "), ".", sep=""), sep="")
          
          paste(as.character(factor), ".", as.character(factor.group), "\n\nThe combined factor levels are ", as.character(group), ".", sep="")
        } else {
          paste("ERROR!") 
        }
        
      })
    }
    
  })
  
  ##End data input tab panel
  ################################################
  
  ################################################
  ###tab panel Data summary
  output$orgLibsizeNormfactor <- renderTable({ 
    tab <- datareactive()$samples[,-1]
    colnames(tab) <- c("Library sizes", "Normalization factors")
    tab
  }, align="l|cc", digits=c(0,0,2), display=c("s", "e", "f"))
  
  output$rmlowLibsizeNormfactor <- renderTable({ 
    input$rmlow | input$dataSubmitMulti
    isolate({ 
      tab <- (rmlowReactive())$samples[, -1]
      colnames(tab) <- c("Library sizes", "Normalization factors")
      tab
    })
  }, align="l|cc", digits=c(0,0,2), display=c("s", "e", "f"))
  
  output$orgSamplesize <- renderTable({ 
    
    no.samples <- length(colnames(datareactive()$counts))
    no.gene <- dim((datareactive())$counts)[1]
    res.summary <- rbind(no.samples, no.gene)
    rownames(res.summary) <- c("Samples", "Tags")
    colnames(res.summary) <- "Number"
    res.summary
    
  },digits=0, align="l|c")
  
  output$rmlowSamplesize <- renderTable({ 
    input$rmlow | input$dataSubmitMulti
    isolate({ 
      no.samples <- length(colnames(rmlowReactive()$counts))
      no.gene <- dim(rmlowReactive()$counts)[1]
      res.summary <- rbind(no.samples, no.gene)
      rownames(res.summary) <- c("Samples", "Tags")
      colnames(res.summary) <- "Number"
      res.summary
    })
  },digits=0, align="l|c")
  
  output$sampleBoxplot <- renderPlot({ 
    input$rmlow | input$dataSubmitMulti
    isolate({
      Group <- as.factor(rmlowReactive()$samples$group)
      bx.p<-boxplot(cpm(rmlowReactive(), log=T)[,])
      bxp(bx.p, boxfill=as.numeric(Group)+1, 
          cex.axis=1.5, whisklwd=3, outcol=as.numeric(Group)+1, 
          main="Normalized sample distribution", cex.main=2)
    }) 
  })
  
  output$sampleMDS <- renderPlot({ 
    input$rmlow | input$dataSubmitMulti
    isolate({
      Group <- as.factor(rmlowReactive()$samples$group)
      par(mar=c(5,5,4,2))
      plotMDS(rmlowReactive(), col=as.numeric(Group)+1, 
              cex=2, main="MDS plot", ndim=3, gene.selection="common",
              xlab = "logFC dim 1", ylab="logFC dim 2", cex.lab=2, cex.main=2, cex.axis=1.5)
    })
  }) 
  ##End Data summary tab Panel
  ################################################
  
  ################################################
  ###tab panel edgeR
  output$edgerGroupLevel <- renderText({ 
    paste("The available group levels are: " ,paste(as.character(levels((datareactive())$samples$group)), collapse=", "), sep="")
  })
  
  ##Estimation dispersion BCV
  output$edgerBCV <- renderPlot({ 
    input$edgerdeAnalysis
    isolate({ 
      par(mar=c(5,5,2,2))
      plotBCV(edgerDispersionEst(), cex=0.5, cex.lab=1.8, cex.axis=1.5)
    })
  })
  
  output$edgerCommonDisp <- renderText({
    paste("Esitmated biological coefficient of variation (BCV) is ", round(sqrt(edgerDispersionEst()$common.dispersion)*100, 4), "%", sep="")
  })
  
  output$edgerTagwiseDispExp <- renderText({
    paste("Esitmated tagwise dispersion can be summarized as below:")
  })
  
  output$edgerTagwiseDisp <- renderTable({ 
    input$edgerdeAnalysis
    isolate({ 
      res<- as.table(t(summary(sqrt(edgerDispersionEst()$tagwise.dispersion))))
      res <- t(res)
      colnames(res) <- "Tagwise"
      res
    })
  }, digits=3, align = "l|c")
  ##############
  ##tab panel for DE analysis results 
  output$edgerRes <- DT::renderDataTable({ 
    input$edgerdeAnalysis
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
    paste(paste(as.character(input$edgercompGroup2), as.character(input$edgercompGroup1), sep="-"), " DE analysis", sep="")
  })
  
  output$edgerTestDGE <- renderTable({ 
    input$edgerdeAnalysis
    isolate({ 
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
    input$edgerdeAnalysis
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
  output$voomGroupLevel <- renderText({ 
    paste("The available group levels are: " ,paste(as.character(levels((datareactive())$samples$group)), collapse=", "), sep="")     
  })
  
  output$voomBCV <- renderPlot({
    Group <- rmlowReactive()$samples$group
    f <- factor(Group, levels=levels(Group))
    design <- model.matrix(~0+f)
    rownames(design) <- rownames(rmlowReactive()$samples)
    colnames(design) <- levels(Group)
    par(mar=c(4,4,2,2))
    voom(rmlowReactive(), design=design, plot=T, normalize="quantile")
  })
  
  output$voomRes <- DT::renderDataTable({
    input$voomdeAnalysis
    isolate({ 
      tpvoom <- topTable(voomDEres(), n=Inf, adjust.method="BH", sort.by="p")
      res <- tpvoom[,c("logFC","P.Value", "adj.P.Val")]
      res$logFC <- round(res$logFC, digits=3)
      res
    })
  }, 
  options = list(order = list(2, 'asc'), columnDefs = list(list(className = 'dt-center', targets = c(1,2,3))), pageLength = 20),  
  colnames = c("log2FC"=2, "p"=3, "FDR"=4, "Tag/Gene Name"=1)
  )
  
  output$voomTestDGEtitle <- renderText({
    paste(paste(as.character(input$voomcompGroup2), as.character(input$voomcompGroup1), sep="-"), " DE analysis", sep="")
  })
  
  output$voomTestDGE <- renderTable({
    input$voomdeAnalysis
    isolate({
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
    input$voomdeAnalysis
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
  output$deseq2GroupLevel <- renderText({ 
    
    paste("The available group levels are: " ,paste(as.character(levels((datareactive())$samples$group)), collapse=", "), sep="")
    
  })
  
  output$deseq2BCV <- renderPlot({
    input$deseq2deAnalysis
    isolate({
      par(mar=c(4,5,2,2))
      plotDispEsts(deseq2Res(), cex.lab=1.8, cex.axis=1.5)
    })
  })
  
  output$deseq2Res <- DT::renderDataTable({
    input$deseq2deAnalysis
    isolate({
      res.pre <- deseq2DEres()[,c("log2FoldChange", "pvalue", "padj")]
      res <- as.data.frame(res.pre) 
      res$log2FoldChange <- round(res$log2FoldChange, digits=3)
      res
    })
  }, 
  options = list(order = list(2, 'asc'), columnDefs = list(list(className = 'dt-center', targets = c(1,2,3))), pageLength = 20),  
  colnames = c("log2FC"=2, "p"=3, "FDR"=4,"Tag/Gene Name"=1)
  )
  
  output$deseq2TestDGEtitle <- renderText({
    paste(paste(as.character(input$deseq2compGroup2), as.character(input$deseq2compGroup1), sep="-"), " DE analysis", sep="")
  })
  
  output$deseq2TestDGE <- renderTable({
    input$deseq2deAnalysis
    isolate({
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
    input$deseq2deAnalysis
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
  output$compGroupLevel <- renderText({ 
    
    paste("The available group levels are: " ,paste(as.character(levels((datareactive())$samples$group)), collapse=", "), sep="")
        
  })
  
  output$decomp <- renderPlot({
    input$decompAnalysis
    isolate({      
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
    input$decompAnalysis
    isolate({
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
    input$decompAnalysis
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
    input$decompAnalysis
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
      res
    })
  },digits = 0, align="l|c")
  
}) 

