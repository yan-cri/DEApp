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


header <- dashboardHeader(
  title = "DE analysis App"
)

sidebar <- dashboardSidebar( 
  sidebarMenu(id="menu1",
              menuItem("Introduction", icon = icon("user"),
                       menuSubItem("Introduction", tabName = "intro"),
                       menuSubItem("Analysis Workflow", tabName = "workflow"),
                       menuSubItem("Q & A", tabName ="QuestionsAndAnswers")
                       ),
              menuItem("Data Input", icon = icon("file"),
                       menuSubItem("Single-factor Experiment", tabName = "dataInputSingle"),
                       menuSubItem("Multi-factor Experiment", tabName ="dataInputMulti")
                       ),
              menuItem("Data Summarization", tabName="dataSummary", icon = icon("table")),
              menuItem("DE analysis", icon = icon("area-chart"), 
                       menuSubItem("edgeR", tabName="edger"),
                       menuSubItem("Limma-Voom", tabName="limmavoom"),
                       menuSubItem("DEseq2", tabName="deseq2")
              ),
              menuItem("Methods Comparison", icon = icon("bar-chart-o"),tabName="decomp"),
              menuItem("Feedback", icon = icon("comment", lib="glyphicon"), tabName="feedback")
              
  )
  ,helpText("Developed by bioinformatics core, Center for Research Informatics (CRI), University of Chicago", style="padding-left:1em; padding-right:1em;position:absolute; bottom:1em; ")
  #,helpText("Developed by bioinformatics core, Center for Research Informatics (CRI), University of Chicago", style="padding-left:1em; padding-right:1em;position:absolute; bottom:5em; ")
  #,img(src="CRI_Logo_Text.png", style="height:2%;margin-left:3em;position:absolute; bottom:1em;")
  
)


body <- dashboardBody(
  tabItems(
    #########################################
    ## Introduction tab panel
    tabItem(tabName="intro",
            #img(src="cri.png"),
            h2("Introduction"),
            p("This interactive web application (DEApp) is developed in R with Shiny to 1). conduct differential 
              expression (DE) analysis with ",
              a("edgeR, ", href="http://bioconductor.org/packages/release/bioc/html/edgeR.html"),
              a("limma-voom, ", href="http://www.genomebiology.com/2014/15/2/R29"), "and ",
              a("DESeq2 ", href="http://bioconductor.org/packages/release/bioc/html/DESeq2.html"),
              " based on the provided input data (raw count results and experimental design information); 
              2). cross-validate the DE analysis results with these 3 different DE analysis methods."
              , style="padding-left: 0em"),
            p("The goal of this App is to provide biologists with an easy way to 
              conduct and cross-validate DE analysis with 3 different methods on their own data."
              , style="padding-left: 0em"),
            #####################################################
            ##Data input section
            h3("1. Data input", style="padding-left: 1em"),
            p("The input data of this App is 2 files in '.txt' or '.csv' format, 
              they are named as 'Raw Count Data' and 'Meta-data Table'."
              , style="padding-left: 2em"),
            h4("1.1 Raw Count Data", style="padding-left: 3em"),            
            p("This input file should contain summarized count results of all samples in the experiment, 
              an example of the expected input data format is presented as below:"
              , style="padding-left: 5em"),
            tableOutput("countTabSamp"),
            tags$style("#countTabSamp table {border: 1px solid black; align: left; margin-left: 6em}","#countTabSamp th {border: 1px solid black;}","#countTabSamp td {border: 1px solid black;}"),

            p("Where, columns correspond to samples, rows correspond to mapped genomic features 
              (e.g. genes, exons, transcript, miRNAs etc.)."
              , style="padding-left: 5em"),
            p("An example of demo 'Raw Count Data' input text file for single-factor experiment 
              used in this App is provided in the 'data' folder named as 'pnas-count_singleFactor.txt', 
              it is also accessible ", 
              a("here.", href=as.character(paste("file://~",getwd(),"/data/pnas-count_singleFactor.txt", sep="")))
              , style="padding-left: 5em"),
            
            h4("1.2 Meta-data Table", style="padding-left: 3em"),
            p("This input data contains summarized experimental design information for each sample. 
              This App is able to conduct DE analysis of both single-factor and multi-factor experiments, 
              and the experiment design information is illustrated on the below 'Meta-data Table':"
              , style="padding-left: 5em"),
            h5("1.2.1 Single-factor Experiment", style="padding-left: 5em; font-weight: bold"),
            p("If the experiment has one single experimental factor, such as 'Group', 
              the input 'Meta-data Table' file should be prepared as below:"
              , style="padding-left: 5em"),
            tableOutput("metaTabSamp"),
            tags$style("#metaTabSamp table {border: 1px solid black; align: left;margin-left: 6em}","#metaTabSamp th {border: 1px solid black;}","#metaTabSamp td {border: 1px solid black;}"),
            p("An example of corresponding demo 'Meta-data Table' text file for single-factor experiment 
              used in this App is provided in the 'data' folder named as 'pnas-count_singleFactor-meta.txt', it is also available "
              , a("here.", href=as.character(paste("file://~",getwd(),"/data/pnas-count_singleFactor-meta.txt", sep=""))) 
              , style="padding-left: 5em"),
            
            h5("1.2.2 Multi-factor Experiment", style="padding-left: 5em; font-weight: bold"),
            p("For multi-factor experiment, the 'Meta-data Table' should look as below:"   
              , style="padding-left: 5em"),
            tableOutput("multimetaTabSamp22"),
            tags$style("#multimetaTabSamp22 table {border: 1px solid black; align: left;margin-left: 6em}","#multimetaTabSamp22 th {border: 1px solid black;}","#multimetaTabSamp22 td {border: 1px solid black;}"),
            
            p("An example of the 'Meta-data Table' csv file for multi-factor experiment used in this App 
              is provided in the 'data' folder named as 'ReadCounts-Chen-edgeRSpringer-multiFactor-meta.csv',
              it is also available "
              , a("here.", href=as.character(paste("file://~",getwd(),"/data/ReadCounts-Chen-edgeRSpringer-multiFactor-meta.csv", sep=""))) 
              , "The corresponding 'Raw Count Data' csv file is also provided in the 'data' folder named as 
              'ReadCounts-Chen-edgeRSpringer-multiFactor.csv', it is accessible "
              , a("here.", href=as.character(paste("file://~",getwd(),"/data/ReadCounts-Chen-edgeRSpringer-multiFactor.csv", sep=""))) 
              , style="padding-left: 5em"),
            
            #####################################################
            ##Filter low expression tags section
            h3("2. Filter low expression tags", style="padding-left: 1em"),
            p("For the 'Data Summarization' section in this App, it is aiming to filter out genetic features with
              very low counts. The guideline of this step is to keep genetic features which are expressed in 
              at least one sample out of each factor level."
              , style="padding-left: 2em"),
            p("For example, if there are 3 combined factor levels in the experiment, 
              then at least 3 samples should have expression level 
              above the expression cutoff value presented in count per million (CPM)."
              , style="padding-left: 2em"),
            withMathJax(),
            p("The expression cutoff (CPM) value is determined according to the library size 
              and normalization factors with formula $$\\text{CPM} = \\frac{\\text{counts}}{\\text{library size} * \\text{normalization factors} * 10^{-6}}$$ 
              For example, if the expression cutoff CPM value is 10, 
              the library size and normalization factors are estimated approximately equal to \\(\\ 2 \\text{ x} 10 ^ 6\\) and 1 for majority samples, 
              then 10 CPM expression cutoff corresponds to about 20 read counts. 
              Therefore, in this example genetic features in more than 3 samples have less than 
              20 read counts (10 CPM) is classified as low expression genetic features 
              and are removed for further downstream DE analysis."
              , style="padding-left: 2em"),
            #####################################################
            ##DE analysis methods section
            h3("3. DE analysis methods", style="padding-left: 1em"),
            p("This App implements 3 different methods to conduct DE analysis, 
              they are ",
              a("edgeR, ", href="http://bioconductor.org/packages/release/bioc/html/edgeR.html"),
              a("limma-voom, ", href="http://www.genomebiology.com/2014/15/2/R29"), "and ",
              a("DESeq2 ", href="http://bioconductor.org/packages/release/bioc/html/DESeq2.html"),
              ", which are implemented in the 'DE analysis' section. 
              For the single-factor experiment, the DE analysis could be conducted between any 
              2 levels of that single-factor; 
              for the multi-factor experiment, the DE analysis could be conducted in a way 
              to combine all the experimental factors into combined factors, 
              so that DE analysis can be conducted in any 2 selected combined factor levels."
              , style="padding-left: 2em")
            #h4("1. edgeR", style="padding-left: 3em"),
            #p("The multi-factor experiment DE analysis is conducted in a way to combine all the experimental 
            #  factors into tone combined factor, followed with the DE analysis in any chosen combined factor levels."
            #  , style="padding-left: 4em"),
            #h4("2. Limma-Voom", style="padding-left: 3em"),
            #h4("3. DESeq2", style="padding-left: 3em"),
            
            #p("Due to the algorithm differences of different DE analysis methods, the method comparison is 
            #      conducted based on the single-factor experimental design and multi-factor experimental design.
            #      Resutls of method comparison can be seen with the 'Methods Comparison' tab button"),
            
            ),
    
    tabItem(tabName="workflow",
            ############################################
            ##Workflow section
            h3("Analysis Workflow"),

            p(strong("Step 1: "), "Upload your input data ('Raw Count Table' and 'Meta-data Table')
              via 'Data Input' section panel for single-factor or multi-factor experiment, 
              a summary of your input data will be presented."
              , style="padding-left: 2em; padding-top: 1em"),
            p(strong("Step 2: "), "Filter out the low expression genetic features via 'Data Summarization' section panel, 
              summarized count results after filtering will be presented here."
              , style="padding-left: 2em"),
            p(strong("Step 3: "), "Conduct DE analysis on the 'Data Analysis' section."
              , style="padding-left: 2em"),
            p(strong("Step 4: "), "Compare/cross-validate DE analysis results via 'Methods Comparison' section panel."
              , style="padding-left: 2em"),
            #includeHTML(path=paste(getwd(),"www/Introduction.html",sep="/"))
            #imageOutput(paste(getwd(),"www/Introduction.pdf",sep="/"))
            #p("The analysis workflow is summarized as below:"
            #  , style="padding-left: 0em, padding-bottom: 4em"),
          
            #img(src="analysis-workflow.png", width=800, style="display: block; margin-left: auto; margin-right: auto;"),
            p("The analysis execution workflow of this App is illustrated in a pdf file, which can be downloaded from  ",
              a("here.", href=as.character(paste("file://~",getwd(),"/www/workflow.pdf", sep="")))
              , style="padding-left: 1em; padding-top: 1em")
            ),
    
    tabItem(tabName="QuestionsAndAnswers",
            ############################################
            ##Q&A section
            h3("Frequently Asked Questions and Answers"),
            br(),
            p("Q1. Can I use this App to analyze RPKM data, such as quantified results from cuffquant?"
              , style="padding-left: 2em; font-weight: bold"),
            p("A: No, this App can only be used to analyze raw count data via edgeR, limma-voom, and DESeq2. 
              The recommended 'Raw Count Data' input is usually obtained from ",
              a("'HTSeq'", href="http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html"), " or ", 
              a("'featureCounts'", href="http://bioinf.wehi.edu.au/featureCounts/"), "program."
              , style="padding-left: 4em"),
            
            p("Q2. What kinds of statistical analysis methods are implemented in this App for DE analysis?"
              , style="padding-left: 2em; font-weight: bold"),
            p("A: This App applied 3 different methods including ",
              a("edgeR, ", href="http://bioconductor.org/packages/release/bioc/html/edgeR.html"),
              a("limma-voom, ", href="http://www.genomebiology.com/2014/15/2/R29"), "and ",
              a("DESeq2 ", href="http://bioconductor.org/packages/release/bioc/html/DESeq2.html"),
              "to conduct DE analysis."
              , style="padding-left: 4em"),
            
            p("Q3. Can I use this App to conduct DE analysis from a multi-factor experiment?"
              , style="padding-left: 2em; font-weight: bold"),
            p("A: Yes, this App supports both single factor and multi-factor experiment's DE analysis, 
              which reflects in the 'Meta-data Table' input. 
              For the multi-factor experiment, the DE analysis is conducted in a way 
              to combine all the experimental factors into one combined factor, 
              so that DE analysis can be conducted in any 2 chosen combined factor levels.
              Additionally, please make sure there is more than 1 biological replicate (>=2 samples) 
              for each factor/combined-factor levels."
              , style="padding-left: 4em")
            
    ),
    
    tabItem(tabName='feedback',
            ############################################
            ##Feedback section
            h2("Feedback"),
            p("This App is developed and maintained by Yan Li at the bioinformatics core, Center for Research Informatics (CRI), 
              Biological Science Division (BSD), University of Chicago."),
            p("As a bioinformatics core, we are actively improving and expanding our NGS analysis services and analysis products.
              If you have any questions, comments, or suggestions, feel free to contact our core at bioinformatics@bsd.uchicago.edu or the developer at yli22@bsd.uchicago.edu"),
            br()#,
            #img(src="cri.png")
            ),
    
    ## End introduction tab panel
    #########################################    
    
    #########################################
    ## First tab content for data input
    tabItem(tabName="dataInputSingle",
            
            fluidRow(
              box(title = "Input data",
                  solidHeader = T, status = "info",
                  collapsible = T, collapsed = F,
                  width = 12,
                  
                  fluidRow(
                    ##Raw count input box under Data Input tab panel
                    box(title = "Input 1: Raw Count Data",
                        solidHeader = T, status = "info",
                        width = 6,
                        helpText("Upload your 'Raw Count Data' here, if no file is selected, 
                                 the demo file for single-factor experiment will be used and displayed."
                                 ,style="color:black; padding-right:0em;"),
                        
                        fileInput(inputId="countFile", 
                                  label=NULL, 
                                  accept=c('text/tab-separated-values',
                                           'text/csv',
                                           'text/comma-separated-values',
                                           'text/tab-separated-values',
                                           '.txt',
                                           '.csv',
                                           '.tsv')
                        ),
                        radioButtons(inputId="coutFileSep", 
                                     label="Separator",
                                     choices=c(Comma=',',
                                               Semicolon=';',
                                               Tab='\t'
                                     ),
                                     selected='\t'
                        ),
                        
                        helpText("The demo file of 'Raw Count Data' for the single-factor experiment is available ", 
                                 a("here", href=as.character(paste("file://~",getwd(),"/data/pnas-count_singleFactor.txt", sep=""))),
                                 style="color: black")
                    ),
                    ##Meta-data input box under Data Input tab panel
                    box(title = "Input 2: Meta-data Table",
                        solidHeader = T, status = "info",
                        width = 6,
                        helpText("Upload your 'Meta-data Table' here, if no file is selected, 
                                 the corresponding demo file for single-factor experiment 
                                 will be used and displayed."
                                 ,style="color:black; padding-right:0em;"),
                        
                        fileInput(inputId="metaTab", 
                                  label=NULL,
                                  accept = c('text/tab-separated-values',
                                             'text/csv',
                                             'text/comma-separated-values',
                                             'text/tab-separated-values',
                                             '.txt',
                                             '.csv',
                                             '.tsv')
                        ), 
                        radioButtons(inputId="metaSep", 
                                     label="Separator",
                                     choices=c(Comma=',',
                                               Semicolon=';',
                                               Tab='\t'
                                     ),
                                     selected='\t'
                        ),
                        
                        helpText("The corresponding 'Meta-data Table' of the demo file for the single factor experiment is accessible ", 
                                 a("here", href=as.character(paste("file://~",getwd(),"/data/pnas-count_singleFactor-meta.txt", sep="")))
                                 ,style="color:black;")
                    )
                  ),

                  actionButton("dataSubmit", label = "Submit"),
                  tags$style("button#dataSubmit {margin-left:auto;margin-right:auto;display:block; background-color:#00CCFF; padding: 5px 25px; font-family:Andika, Arial, sans-serif; font-size:1.5em;  letter-spacing:0.05em; text-transform:uppercase ;color:#fff; text-shadow: 0px 1px 10px #000;border-radius: 15px;box-shadow: rgba(0, 0, 0, .55) 0 1px 6px;}"),
                  helpText(textOutput("errorInputSingle"), style="color:red;")
                  
                    )
                    ),
            
            ##Input information summary under Data Input Panel
            fluidRow(
              box(title = "Input Information Summary",
                  solidHeader = T, status = "success",
                  #footer = "Note: There must be a biological replicates from the experiment for DE analysis",
                  width = 12, 
                  fluidRow(
                    column(6, 
                           h4(textOutput("sampleTitle")),
                           verbatimTextOutput("sampleInfo")
                    ),
                    column(6,
                           h4(textOutput("expDesign")),
                           verbatimTextOutput("GroupLevel")
                    )
                  )
              )
            ),
            
            fluidRow(
              
              ##Input Raw Counts Summary box under Data Input Panel
              box(title = "Input Raw Count Summary",
                  solidHeader = T, status = "success",
                  #background = "green",
                  width = 6,
                  tableOutput("overallDataSummary"),
                  tags$style("#overallDataSummary table {border: 1px solid black; align: center; margin:auto;}","#overallDataSummary th {border: 1px solid black;}","#overallDataSummary td {border: 1px solid black;}")
              ),
              
              ##Sample Group Information Summary box under Data Input Panel
              box(title = "Sample Group Information Summary",
                  solidHeader = T, status = "success",
                  #background = "green",
                  width = 6,
                  tableOutput("sampleGroup"),
                  tags$style("#sampleGroup table {border: 1px solid black; align: center; margin:auto;}","#sampleGroup th {border: 1px solid black;}","#sampleGroup td {border: 1px solid black;}")
              )
            )
              ),
    
    tabItem(tabName="dataInputMulti",
            fluidRow(
              box(title = "Input data",
                  solidHeader = T, status = "info",
                  collapsible = T, collapsed = F,
                  width = 12,
                  
                  fluidRow(
                    ##Raw count input box under Data Input tab panel
                    box(title = "Input 1: Raw Count Data",
                        solidHeader = T, status = "info",
                        width = 6,
                        
                        helpText("Upload your 'Raw Count Data' here, if no file is selected, 
                                 the demo file for multi-factor experiment will be used and displayed."
                                 ,style="color:black; padding-right:0em;"),
                        
                        fileInput(inputId="countFileMulti", 
                                  label=NULL, 
                                  accept=c('text/tab-separated-values',
                                           'text/csv',
                                           'text/comma-separated-values',
                                           'text/tab-separated-values',
                                           '.txt',
                                           '.csv',
                                           '.tsv')
                        ),
                        radioButtons(inputId="coutFileSepMulti", 
                                     label="Separator",
                                     choices=c(Comma=',',
                                               Semicolon=';',
                                               Tab='\t'
                                     ),
                                     selected=','
                        ),
                        
                        helpText(HTML("<div style=\"color: black \">An example of full 'Raw Count Data' csv file for the multi-factor experiment is accessible "), 
                                 a(HTML("here</div>"), href=as.character(paste("file://~",getwd(),"/data/ReadCounts-Chen-edgeRSpringer-multiFactor.csv", sep="")))
                        )
                    ),
                    ##Meta-data input box under Data Input tab panel
                    box(title = "Input 2: Meta-data Table",
                        solidHeader = T, status = "info",
                        width = 6,
                        
                        helpText("Upload your 'Meta-data Table' here, if no file is selected, 
                                 the corresponding demo file for multi-factor experiment 
                                 will be used and displayed."
                                 ,style="color:black; padding-right:0em;"),
                        
                        fileInput(inputId="metaTabMulti", 
                                  label=NULL,
                                  accept = c('text/tab-separated-values',
                                             'text/csv',
                                             'text/comma-separated-values',
                                             'text/tab-separated-values',
                                             '.txt',
                                             '.csv',
                                             '.tsv')
                        ), 
                        radioButtons(inputId="metaSepMulti", 
                                     label="Separator",
                                     choices=c(Comma=',',
                                               Semicolon=';',
                                               Tab='\t'
                                     ),
                                     selected=','
                        ),
                        
                        helpText("An example of 'Meta-data Table' csv file corresponding to the example of input 1 - 'Raw Count Data' for the multi-factor experiment is accessible ",
                                 a("here", href=as.character(paste("file://~",getwd(),"/data/ReadCounts-Chen-edgeRSpringer-multiFactor-meta.csv", sep="")))
                                 , style="color: black")
                        
                        #helpText(HTML("<div style=\"color: black; padding-left: 1em; padding-right: 1em; margin-top: 1cm \"> An example of 'Meta-data Table' csv file corresponding to the example of input 1 - 'Raw Count Data' for the multi-factor experiment is accessible ")
                        #         , a(HTML("here</div>"), href=as.character(paste("file://~",getwd(),"/data/ReadCounts-Chen-edgeRSpringer-multiFactor-meta.csv", sep="")))
                        #)
                    )
                    
                  ),

                  actionButton("MultiSubmit", label = "Submit"),
                  tags$style("button#MultiSubmit {margin-left:auto;margin-right:auto;display:block; background-color:#00CCFF; padding: 5px 25px; font-family:Andika, Arial, sans-serif; font-size:1.5em;  letter-spacing:0.05em; text-transform:uppercase ;color:#fff; text-shadow: 0px 1px 10px #000;border-radius: 15px;box-shadow: rgba(0, 0, 0, .55) 0 1px 6px;}"),

                  helpText(textOutput("errorInputMulti"), style="color:red;")

              )
            ),    
            
            fluidRow(
              box(title = "Input Information Summary",
                  solidHeader = T, status = "success",
                  #footer = "Note: There must be a biological replicates from the experiment for DE analysis",
                  width = 12, 
                  fluidRow(
                    column(6, 
                           h4(textOutput("sampleTitleMulti")),
                           verbatimTextOutput("sampleInfoMulti")
                    ),
                    column(6,
                           h4(textOutput("expDesignMulti")),
                           verbatimTextOutput("GroupLevelMulti")
                    )
                  )
              )
            ),

            
            fluidRow(
              
              ##Input Raw Counts Summary box under Data Input Panel
              box(title = "Input Raw Count Summary",
                  solidHeader = T, status = "success",
                  #background = "green",
                  width = 6,
                  tableOutput("overallDataSummaryMulti"),
                  tags$style("#overallDataSummaryMulti table {border: 1px solid black; align: center; margin:auto;}","#overallDataSummaryMulti th {border: 1px solid black;}","#overallDataSummaryMulti td {border: 1px solid black;}")
              ),
              
              ##Sample Group Information Summary box under Data Input Panel
              box(title = "Sample Group Information Summary",
                  solidHeader = T, status = "success",
                  #background = "green",
                  width = 6,
                  tableOutput("sampleGroupMulti"),
                  tags$style("#sampleGroupMulti table {border: 1px solid black; align: center; margin:auto;}","#sampleGroupMulti th {border: 1px solid black;}","#sampleGroupMulti td {border: 1px solid black;}")
              )
            )
            
              ),
            
    
    ## End first tab content for data input panel
    #########################################
    
    #########################################
    ## Second tab content for data summarization panel
    tabItem(tabName="dataSummary",
            
            ##First 3 box under Data Summarization panel
            fluidRow(
              box(title = "Low Expression Removal",
                  solidHeader = T, status = "info",
                  collapsible = T, collapsed = F,
                  width = 12,
                  
                  fluidRow(
                    ##Raw Count Summary box under the Data Summarization panel
                    box(title = "Raw Count Summary",
                        #background = "yellow",
                        solidHeader = T, status = "warning",
                        width = 5,
                        fluidRow(
                          column(7,
                                 tableOutput("orgLibsizeNormfactor"),
                                 tags$style("#orgLibsizeNormfactor table {border: 1px solid black; float: center; position:relative;}","#orgLibsizeNormfactor th {border: 1px solid black;}","#orgLibsizeNormfactor td {border: 1px solid black;}")
                          ),
                          column(5,
                                 tableOutput("orgSamplesize"),
                                 tags$style("#orgSamplesize table {border: 1px solid black; float: center; position:relative;}","#orgSamplesize th {border: 1px solid black;}","#orgSamplesize td {border: 1px solid black;}")
                          )
                        )               
                    ),
                    
                    ##Low expression removal threshold chosing box under Data Summarization panel
                    box(title = "Low Expression Removal Options",
                        solidHeader = T, status = "info",
                        width = 2,
                        h6("Low expression mapped genetic features will be removed with"),
                        fluidRow(
                          column(6,
                                 textInput(inputId="cpmVal", 
                                           label="CPM value", 
                                           value="1")
                          ),
                          column(6,
                                 textInput(inputId="gThreshold", 
                                           label="at least in No. Samples", 
                                           value="2")
                          )
                        ),
                        actionButton(inputId="rmlow", label="Submit"),
                        tags$style("button#rmlow {margin-left:auto;margin-right:auto;display:block;background-color:#00CCFF; font-family:Andika, Arial, sans-serif; font-size:0.8em; letter-spacing:0.05em; text-transform:uppercase ;color:#fff; text-shadow: 0px 1px 10px #000;box-shadow: rgba(0, 0, 0, .55) 0 1px 6px;}")
                        
                    ),
                    
                    ##After Low Expression Removal box under Data Summarization panel
                    box(title = "After Low Expression Removal",
                        #background = "green",
                        solidHeader = T, status = "success",
                        width = 5,
                        fluidRow(
                          column(7,
                                 tableOutput("rmlowLibsizeNormfactor"),
                                 tags$style("#rmlowLibsizeNormfactor table {border: 1px solid black; float: center; position:relative;}","#rmlowLibsizeNormfactor th {border: 1px solid black;}","#rmlowLibsizeNormfactor td {border: 1px solid black;}")
                          ), 
                          column(5,
                                 tableOutput("rmlowSamplesize"),
                                 tags$style("#rmlowSamplesize table {border: 1px solid black; float: center; position:relative;}","#rmlowSamplesize th {border: 1px solid black;}","#rmlowSamplesize td {border: 1px solid black;}")
                          )
                        )
                    )
                  ),
                  helpText(textOutput("errorFiltering"), style="color:red;")
              )
    ),
    
    ##Second row with 2 box under the Data Summarization panel
    fluidRow(
      
      ## Sample Normalization box plot under Data Summarization panel
      box(title = "Sample Normalization Results",
          solidHeader = T, status = "success",
          width = 6,
          plotOutput("sampleBoxplot")
      ),
      
      ## Sample MDS Exploration plot under Data Summarization panel
      box(title = "Sample MDS Exploration",
          solidHeader = T, status = "success",
          width = 6,
          plotOutput("sampleMDS")
      )
    )
    ),
    ## End second tab content for data summarization panel
    #########################################
    
    #########################################
    ## DE analysis tab panel content for edgeR analysis under DE analysis menu tab
    tabItem(tabName="edger",
            
            ## 1st row with 2 boxes under edgeR tab panel
            fluidRow(          
              box(title = "edgeR DE Analysis Options",
                  solidHeader = T, status = "info",
                  collapsible = T, collapsed = F,
                  width = 12,
                  fluidRow( 
                    ## GLE DE analysis group chosing box under edgeR tab panel
                    box(title = "DE Analysis Group Levels",
                        #footer = "DE analysis is conducted between any 2 specified levels out of available group levels",
                        width = 4,
                        solidHeader = T,
                        status = "info",
                        
                        textOutput("edgerGroupLevel"),
                        tags$style("#edgerGroupLevel{ font-weight: bold; color: #0033FF; padding-top: .3cm; padding-bottom: .3cm;}"),
                        p("Please select any 2 levels from the above available group levels for DE analysis."
                          , style="font-weight: bold"),
                        
                        
                        fluidRow(
                          column(6,
                                 textInput(inputId="edgercompGroup1", 
                                           label="Level 1", 
                                           value="")
                          ),
                          column(6,
                                 textInput(inputId="edgercompGroup2", 
                                           label="Level 2", 
                                           value="")
                          )
                          #helpText("Please specify the factor group level based on your experimental design for DE analysis, 
                          #for the single-factor experiment, it can be any 2 levels of the factor level,
                          #for the multi-factor experiment, it can be any 2 levels of the combined factor levels. 
                          #Please provide correct comparison factor levels for DE analysis, which are summarized on the 'Data Input' tab.")
                        )
                        
                    ),
                    
                    ## DE filtering criteria option box under edgeR tab panel
                    box(title = "DE Analysis Filtering Criteria",
                        width = 8,
                        solidHeader = T,
                        status = "info",
                        fluidRow(
                          column(4,
                                 radioButtons(inputId="edgerP", 
                                              label="DE Analysis is based on", 
                                              choices=c("Nominal p-value" ='normp',
                                                        "FDR adjusted p-value" ='fdrp'
                                              ),
                                              selected='fdrp'
                                 )
                          ),
                          column(4,
                                 textInput(inputId="edgerfdr", 
                                           label="p-value or FDR adjusted p-value", 
                                           value="0.05")
                          ),
                          column(4,
                                 textInput(inputId="edgerfc", 
                                           label="Fold Change (FC)", 
                                           value="1.5")
                          )
                        )
                        
                    ),
                    
                    actionButton(inputId="edgerdeAnalysis", label="Submit"),
                    tags$style("button#edgerdeAnalysis {margin-top:0.5em;float:right; margin-right: 1em; background-color:#00CCFF; padding: 5px 25px; font-family:Andika, Arial, sans-serif; font-size:1.5em;  letter-spacing:0.05em; text-transform:uppercase ;color:#fff; text-shadow: 0px 1px 10px #000;border-radius: 15px;box-shadow: rgba(0, 0, 0, .55) 0 1px 6px;}")
                  ),
                  helpText(textOutput("erroredgeR"), style="color:red;")
              )
            ),
            
            ## 2nd row under edgeR tab panel with 3 columns box contents
            fluidRow(
              
              ## Column 1 with 2 boxes under it
              column(3,
                     
                     ## Estimated BCV summay box under edgeR tab panel
                     box(title = "Estimated BCV Summary",
                         width = NULL,
                         #background = "green",
                         solidHeader = T, status = "success",
                         verbatimTextOutput("edgerCommonDisp"),
                         h5(textOutput("edgerTagwiseDispExp")),
                         tableOutput("edgerTagwiseDisp"),
                         tags$style("#edgerTagwiseDisp table {border: 1px solid black; align: center; margin:auto;}","#edgerTagwiseDisp th {border: 1px solid black;}","#edgerTagwiseDisp td {border: 1px solid black;}")
                     ),
                     
                     ## Estimated dispersion under edgeR tab panel
                     box(title = "Estimated Dispersion",
                         width = NULL,
                         solidHeader = T, status = "success",
                         plotOutput("edgerBCV")
                     )
              ),
              
              ## Column 2 with only 1 box for DE analysis results from edgeR
              ## The results include all gene analysis results
              column(6,
                     box(title = "DE Analysis Results",
                         width = NULL,
                         solidHeader = T, status = "success",
                         DT::dataTableOutput("edgerRes")
                     )
              ),
              
              ## Column 3 with 2 boxes under it
              column(3,
                     
                     ## DEGs Summary to summarize the up- and down-regulated DEGs
                     box(title = "DE Results Summary",
                         width = NULL,
                         #background = "green",
                         solidHeader = T, status = "success",
                         h4(textOutput("edgerTestDGEtitle"), align="center"),
                         tableOutput("edgerTestDGE"),
                         tags$style("#edgerTestDGE table {border: 1px solid black; align: center; margin:auto;}","#edgerTestDGE th {border: 1px solid black;}","#edgerTestDGE td {border: 1px solid black;}"),
                         br(),
                         downloadButton("edgerDownload", 
                                        label = "Download",
                                        class = NULL),
                         tags$style("#edgerDownload {float:right; }")
                         #helpText("Download DEGs results", align="right")
                     ),
                     
                     ## Volcano plot of DE analysis results
                     box(title = "Volcano Exploration",
                         width = NULL,
                         solidHeader = T, status = "success",
                         plotOutput("edgerVolcano")
                     )
              )
            )            
    ),
    ## End edgeR DE analysis tab panel under DE analysis menu tab
    #########################################
    
    #########################################
    ## DE analysis tab content panel for limma-voom
    tabItem(tabName="limmavoom",            
            ## 1st row with 2 boxes under limma-voom tab panel
            fluidRow(
              box(title = "Limma-voom DE Analysis Options",
                  solidHeader = T, status = "info",
                  collapsible = T, collapsed = F,
                  width = 12,
                  fluidRow(
                    ## DE analysis group chosing box under limma-voom tab panel
                    box(title = "DE Analysis Group Levels",
                        #footer = "DE analysis is conducted between any 2 specified levels out of available group levels",
                        width = 4,
                        solidHeader = T,
                        status = "info",
                        
                        textOutput("voomGroupLevel"),
                        tags$style("#voomGroupLevel{ font-weight: bold; color: #0033FF; padding-top: .3cm; padding-bottom: .3cm;}"),
                        p("Please select any 2 levels from the above available group levels for DE analysis."
                          , style="font-weight: bold"),
                        
                        fluidRow(
                          column(6,
                                 textInput(inputId="voomcompGroup1", 
                                           label="Group 1", 
                                           value="")
                          ),
                          column(6,
                                 textInput(inputId="voomcompGroup2", 
                                           label="Group 2", 
                                           value="")
                          )
                        )
                    ),
                    
                    ## DE filtering criteria option box under limma-voom tab panel
                    box(title = "DE Analysis Filtering Criteria",
                        width = 8,
                        solidHeader = T,
                        status = "info",
                        fluidRow(
                          column(4,
                                 radioButtons(inputId="voomP", 
                                              label="DE Analysis is based on", 
                                              choices=c("Nominal p-value" ='normp',
                                                        "FDR adjusted p-value" ='fdrp'
                                              ),
                                              selected='fdrp'
                                 )
                          ),
                          column(4,
                                 textInput(inputId="voomfdr", 
                                           label="p-value or FDR adjusted p-value", 
                                           value="0.05")
                          ),
                          column(4,
                                 textInput(inputId="voomfc", 
                                           label="Fold Change (FC)", 
                                           value="1.5")
                          )
                        )
                        
                    ),
                    actionButton(inputId="voomdeAnalysis", label="Submit"),
                    tags$style("button#voomdeAnalysis {margin-top:0.5em;float:right; margin-right: 1em; background-color:#00CCFF; padding: 5px 25px; font-family:Andika, Arial, sans-serif; font-size:1.5em;  letter-spacing:0.05em; text-transform:uppercase ;color:#fff; text-shadow: 0px 1px 10px #000;border-radius: 15px;box-shadow: rgba(0, 0, 0, .55) 0 1px 6px;}")
                  ),
                  helpText(textOutput("errorVoom"), style="color:red;")
              )
            ),
            
            ## 2nd row under edgeR tab panel with 3 columns box contents
            fluidRow(
              
              ## Column 1 with Estimated Dispersion box under it
              column(3,                   
                     ## Estimated dispersion under limma-voom tab panel
                     box(title = "Estimated Dispersion",
                         width = NULL,
                         solidHeader = T, status = "success",
                         plotOutput("voomBCV")
                     )
              ),
              
              ## Column 2 with only 1 box for DE analysis results from voom
              ## The results include all gene analysis results
              column(6,
                     box(title = "DE Analysis Results",
                         width = NULL,
                         solidHeader = T, status = "success",
                         DT::dataTableOutput("voomRes")
                     )
              ),
              
              ## Column 3 with 2 boxes under it
              column(3,
                     
                     ## DEGs Summary to summarize the up- and down-regulated DEGs
                     box(title = "DE Results summary",
                         width = NULL,
                         #background = "green",
                         solidHeader = T, status = "success",
                         h4(textOutput("voomTestDGEtitle"), align="center" ),
                         tableOutput("voomTestDGE"),
                         tags$style("#voomTestDGE table {border: 1px solid black; align: center; margin:auto;}","#voomTestDGE th {border: 1px solid black;}","#voomTestDGE td {border: 1px solid black;}"),
                         br(),
                         downloadButton("voomDownload", 
                                        label = "Download",
                                        class = NULL),
                         tags$style("#voomDownload {float:right; }")
                         #helpText("Download DEGs results")
                     ),
                     
                     ## Volcano plot of DE analysis results
                     box(title = "Volcano Exploration",
                         width = NULL,
                         solidHeader = T, status = "success",
                         plotOutput("voomVolcano")
                     )
              )
            )            
    ),
    ## End DE analysis tab panel for limma-voom
    #########################################
    
    #########################################
    ## DE analysis tab content panel for deseq2 analysis
    tabItem(tabName="deseq2",             
            ## 1st row with 2 boxes under deseq2 tab panel
            fluidRow(
              box(title = "DESeq2 DE Analysis Options",
                  solidHeader = T, status = "info",
                  collapsible = T, collapsed = F,
                  width = 12,
                  fluidRow(
                    ## DE analysis group chosing box under deseq2 tab panel
                    box(title = "DE Analysis Group Levels",
                        #footer = "DE analysis is conducted between any 2 specified levels out of available group levels",
                        width = 4,
                        solidHeader = T,
                        status = "info",
                        
                        textOutput("deseq2GroupLevel"),
                        tags$style("#deseq2GroupLevel{ font-weight: bold; color: #0033FF; padding-top: .3cm; padding-bottom: .3cm;}"),
                        p("Please select any 2 levels from the above available group levels for DE analysis."
                          , style="font-weight: bold"),
                        
                        fluidRow(
                          column(6,
                                 textInput(inputId="deseq2compGroup1", 
                                           label="Group 1", 
                                           value="")
                          ),
                          column(6,
                                 textInput(inputId="deseq2compGroup2", 
                                           label="Group 2", 
                                           value="")
                          )
                        )
                    ),
                    
                    ## DE filtering criteria option box under deseq2 tab panel
                    box(title = "DE Analysis Filtering Criteria",
                        width = 8,
                        solidHeader = T,
                        status = "info",
                        fluidRow(
                          column(4,
                                 radioButtons(inputId="deseq2P", 
                                              label="DE Analysis is based on", 
                                              choices=c("Nominal p-value" ='normp',
                                                        "FDR adjusted p-value" ='fdrp'
                                              ),
                                              selected='fdrp'
                                 )
                          ),
                          column(4,
                                 textInput(inputId="deseq2fdr", 
                                           label="p-value or FDR adjusted p-value", 
                                           value="0.05")
                          ),
                          column(4,
                                 textInput(inputId="deseq2fc", 
                                           label="Fold Change (FC)", 
                                           value="1.5")
                          )
                        )
                        
                    ),
                    actionButton(inputId="deseq2deAnalysis", label="Submit"),
                    tags$style("button#deseq2deAnalysis {margin-top:0.5em;float:right; margin-right: 1em; background-color:#00CCFF; padding: 5px 25px; font-family:Andika, Arial, sans-serif; font-size:1.5em;  letter-spacing:0.05em; text-transform:uppercase ;color:#fff; text-shadow: 0px 1px 10px #000;border-radius: 15px;box-shadow: rgba(0, 0, 0, .55) 0 1px 6px;}")
                  ),
                  helpText(textOutput("errorDeseq2"), style="color:red;")
              )
            ),
            
            ## 2nd row under DEseq2 tab panel with 3 columns box contents
            fluidRow(             
              ## Column 1 with Estimated dispersion boxes under it
              column(3,                  
                     ## Estimated dispersion under deseq2 tab panel
                     box(title = "Estimated Dispersion",
                         width = NULL,
                         solidHeader = T, status = "success",
                         plotOutput("deseq2BCV")
                     )
              ),
              
              ## Column 2 with only 1 box for DE analysis results from deseq2
              ## The results include all gene analysis results
              column(6,
                     box(title = "DE Analysis Results",
                         width = NULL,
                         solidHeader = T, status = "success",
                         DT::dataTableOutput("deseq2Res")
                     )
              ),
              
              ## Column 3 with 2 boxes under it
              column(3,
                     
                     ## DEGs Summary to summarize the up- and down-regulated DEGs
                     box(title = "DE Results summary",
                         width = NULL,
                         #background = "green",
                         solidHeader = T, status = "success",
                         h4(textOutput("deseq2TestDGEtitle"), align="center"),
                         tableOutput("deseq2TestDGE"),
                         tags$style("#deseq2TestDGE table {border: 1px solid black; align: center; margin:auto;}","#deseq2TestDGE th {border: 1px solid black;}","#deseq2TestDGE td {border: 1px solid black;}"),
                         br(),
                         downloadButton("deseq2Download", 
                                        label = "Download",
                                        class = NULL),
                         tags$style("#deseq2Download {float:right; }")
                         
                         #helpText("Download DEGs results", style="float:right")
                     ),
                     
                     ## Volcano plot of DE analysis results
                     box(title = "Volcano Exploration",
                         width = NULL,
                         solidHeader = T, status = "success",
                         plotOutput("deseq2Volcano")
                     )
              )
            )            
    ),
    ## End DE analysis tab panel for deseq2 analysis under DE analysis menu tab
    #########################################
    
    #########################################
    ## DE analysis comparison tab panel for single factor experiment
    tabItem(tabName="decomp", 
            fluidRow(
              box(title = "DE Analysis Comparison Options",
                  solidHeader = T, status = "info",
                  collapsible = T, collapsed = F,
                  width = 12,
                  fluidRow(
                    box(title = "Methods for Comparison",
                        width = 3,
                        solidHeader = T,
                        status = "info",
                        checkboxGroupInput(inputId="decompMethods",
                                           label = "DE analysis method selection",
                                           choices = list("edgeR" = 'edger',
                                                          "limma-voom" = 'voom',
                                                          "DESeq2" = 'deseq2'),
                                           selected = c('edger', 'voom')
                        )
                    ),
                    box(title = "DE Analysis Group Levels",
                        width = 3,
                        solidHeader = T,
                        status = "info",
                        
                        textOutput("compGroupLevel"),
                        tags$style("#compGroupLevel{ font-weight: bold; color: #0033FF; padding-top: .3cm; padding-bottom: .3cm;}"),
                        p("Please select any 2 levels from the above available group levels for DE analysis."
                          , style="font-weight: bold"),
                        
                        fluidRow(
                          column(6,
                                 textInput(inputId="decompGroup1", 
                                           label="Group 1", 
                                           value="")
                          ),
                          column(6,
                                 textInput(inputId="decompGroup2", 
                                           label="Group 2", 
                                           value="")
                          )
                        )
                    ),
                    box(title = "DE Analysis Filtering Criteria",
                        width = 6,
                        solidHeader = T,
                        status = "info",
                        fluidRow(
                          column(3,
                                 radioButtons(inputId="decompP", 
                                              label="DE Analysis is based on", 
                                              choices=c("Nominal p-value" ='normp',
                                                        "FDR adjusted p-value" ='fdrp'
                                              ),
                                              selected='fdrp'
                                 )
                          ),
                          column(5,
                                 textInput(inputId="decompfdr", 
                                           label="Nominal p-value or FDR adjusted p-value", 
                                           value="0.05")
                          ),
                          column(4,
                                 textInput(inputId="decompfc", 
                                           label="Fold Change (FC)", 
                                           value="1.5")
                          )
                        )
                    )
                  ),
                  actionButton(inputId="decompAnalysis", label="Submit"),
                  tags$style("button#decompAnalysis {margin-left:auto;margin-right:auto;display:block; background-color:#00CCFF; padding: 5px 25px; font-family:Andika, Arial, sans-serif; font-size:1.5em;  letter-spacing:0.05em; text-transform:uppercase ;color:#fff; text-shadow: 0px 1px 10px #000;border-radius: 15px;box-shadow: rgba(0, 0, 0, .55) 0 1px 6px;}"),
                  helpText(textOutput("errorComp"), style="color:red;")
                  #helpText("Note 1). Please make sure at least 2 methods are selected for DE analysis comparison; 2). DE analysis is conducted between any 2 specified levels out of available group levels.")
              )
            ),
            fluidRow( 
              column(6,
                     box(title = "Comparison Summary",
                         width = NULL,
                         #background = "green",
                         solidHeader = T, status = "success",
                         h5(textOutput("decompTitle")),
                         tableOutput("decompTab"),
                         tags$style("#decompTab table {border: 1px solid black; align: center; margin:auto;}","#decompTab th {border: 1px solid black;}","#decompTab td {border: 1px solid black;}"),
                         br(),
                         verbatimTextOutput("decompText")
                     )
              ),
              column(6,
                     box(title = "Comparison Venn-Diagram",
                         width = NULL,
                         solidHeader = T, status = "success",
                         plotOutput("decomp")
                     )
              )
            )
    )
    ## End DE analysis comparison tab panel for single factor experiment
    #########################################
    
    
)
            )


ui <- dashboardPage(header, sidebar, body, skin = "red")
