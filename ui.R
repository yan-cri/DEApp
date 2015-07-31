## This App is used for RNA-seq DE analysis with different methods together with method comparison
## Developed by Yan Li, last update on June, 2015

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
              menuItem("Feedback", tabName="feedback")
  )
)


body <- dashboardBody(
  tabItems(
    #########################################
    ## Introduction tab panel
    tabItem(tabName="intro",
            img(src="cri.png"),
            h2("Introduction"),
            p("This interactive web application (DEApp) is developed in R with Shiny, it is aiming to 1). conduct differential 
              expression (DE) analysis with edgeR, limma-voom, and DESeq2 based on 
              the provided input data (raw count results and experimental design information); 
              2). cross-validate the DE analysis results with these 3 different DE analysis methods."
              , style="padding-left: 0em"),
            p("The goal of this App development is that biologists could easily 
              conduct and cross-validate DE analysis with 3 different methods with their own data."
              , style="padding-left: 0em"),
            #####################################################
            ##Data input section
            h3("1. Data input", style="padding-left: 1em"),
            p("The input data of this App is 2 files in '.txt' format or '.csv' format, 
              they are named as 'Raw Count Data' and 'Meta-data Table'."
              , style="padding-left: 2em"),
            h4("1.1 Raw Count Data", style="padding-left: 3em"),            
            p("This input data contains summarized count results of all samples in the experiment, 
              an example of the top part of this input data is presented as below:"
              , style="padding-left: 5em"),
            tableOutput("countTabSamp"),
            tags$style("#countTabSamp table {border: 1px solid black; align: left; margin-left: 6em}","#countTabSamp th {border: 1px solid black;}","#countTabSamp td {border: 1px solid black;}"),

            p("Where, columns correspond to the sample, rows correspond to the tag. 
              The 1st column and row indicate the sample ID and tag ID respectively."
              , style="padding-left: 5em"),
            p("An example of full 'Raw Count Data' input plain text file tested in this App is available ", 
              a("here", href=as.character(paste("file://~",getwd(),"/data/TestData-featureCount.txt", sep="")))
              , style="padding-left: 5em"),
            
            h4("1.2 Meta-data Table", style="padding-left: 3em"),
            p("This input data contains summarized experimental design information for each sample. 
              This App is able to conduct DE analysis of both single-factor and multi-factor experiment, 
              and the experiment design information is illustrated in this 'Meta-data Table' input."
              , style="padding-left: 5em"),
            h5("1.2.1 Single-factor Experiment", style="padding-left: 5em; font-weight: bold"),
            p("If the experiment only has only one single experimental factor, such as 'Group'. 
              To illustrate such experimental design, the input 'Meta-data Table' file
              should be prepared as below."
              , style="padding-left: 5em"),
            tableOutput("metaTabSamp"),
            tags$style("#metaTabSamp table {border: 1px solid black; align: left;margin-left: 6em}","#metaTabSamp th {border: 1px solid black;}","#metaTabSamp td {border: 1px solid black;}"),
            
            p("Where, the 1st column corresponds to the sample name, and the 2nd column corresponds 
              to the single experimental factor information - 'Group' for each sample."
              , style="padding-left: 5em"),
            p("An example of the 'Meta-data Table' plain text file for single-factor experiment tested in this App is accessible "
              , a("here", href=as.character(paste("file://~",getwd(),"/data/TestData-featureCount-meta.txt", sep=""))) 
              , style="padding-left: 5em"),
            h5("1.2.2 Multi-factor Experiment", style="padding-left: 5em; font-weight: bold"),
            p("If there are more than one experimental factors used in the experiment, it is a multi-factor experiment,
              such multi-factor information can be listed in the following columns one by one."   
              , style="padding-left: 5em"),
            p("Such as in addition to the 'Group' factor, 
              I have another experimental factor 'Time' used in the experiment, 
              this factor information with respect to each sample can be listed to 
              the 3rd column of the 'Meta-data table'. Similarly, other factor information
              can be summarized in the 'Meta-data table' in the same way. "   
              , style="padding-left: 5em"),
            tableOutput("multimetaTabSamp22"),
            tags$style("#multimetaTabSamp22 table {border: 1px solid black; align: left;margin-left: 6em}","#multimetaTabSamp22 th {border: 1px solid black;}","#multimetaTabSamp22 td {border: 1px solid black;}"),
            
            p("An example of the 'Meta-data Table' text file for multi-factor experiment tested in this App is available "
              , a("here", href=as.character(paste("file://~",getwd(),"/data/ReadCounts-Chen-edgeRSpringer-multiFactor-meta.csv", sep=""))) 
              , style="padding-left: 5em"),
            
            #####################################################
            ##Filter low expression tags section
            h3("2. Filter low expression tags", style="padding-left: 1em"),
            p("This step is implemented with the 'Data Summarization' tab in this App, aiming to filter out tags/genes with
              very low counts. The guideline of this step is to keep tags which are expressed in 
              at least one sample out of each factor levels."
              , style="padding-left: 2em"),
            p("Such as there are 3 factor/combined-factor levels in the experiment, 
              then at least 3 samples should better have expression level 
              above the expression cutoff value presented in count per million (CPM)."
              , style="padding-left: 2em"),
            withMathJax(),
            p("The expression cutoff (CPM) value is determined according to the library size 
              and normalization factors with formula $$\\text{CPM} = \\frac{\\text{counts}}{\\text{library size} * \\text{normalization factors} * 1e^{-6}}$$ 
              Such as the expression cutoff CPM value is 10, 
              the library size and normalization factors are estimated approximately equal to 2e6 and 1 for majority samples, 
              then 10 CPM expression cutoff corresponds to about 20 read counts. 
              Therefore, in this example tags/genes in more than 3 samples have less than 
              20 read counts (10 CPM) is classified as low expression tags and are removed for further downstream DE analysis."
              , style="padding-left: 2em"),
            #####################################################
            ##DE analysis methods section
            h3("3. DE analysis methods", style="padding-left: 1em"),
            p("This App implements 3 different methods to conduct DE analysis, 
              they are ",
              a("edgeR, ", href="http://bioconductor.org/packages/release/bioc/html/edgeR.html"),
              a("limma-voom, ", href="http://www.genomebiology.com/2014/15/2/R29"), "and ",
              a("DESeq2 ", href="http://bioconductor.org/packages/release/bioc/html/DESeq2.html"),
              ". The DE analysis with respect to different methods is implemented  
              with sub-tab button under the 'DE analysis' tab button. 
              For the single-factor experiment, the DE analysis is conducted between any 
              2 levels of that single-factor; 
              for the multi-factor experiment, the DE analysis is conducted in a way 
              to combine all the experimental factors into one combined factor, 
              so that DE analysis can be conducted in any 2 chosen combined factor levels."
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
            p("In this App, different color box is used to distinguish the analysis results, 
              the user input options, the example and pre-processing results."
              , style="padding-left: 0em"),
            
            p("In general, the analysis results are presented in 'Green' color box; 
              the user input options are provided in 'Blue' color box which allows to set up 
              various DE analysis criteria, and it is collapsible by clicking 
              the right top corner '+/-' sign; the example and pre-processing results 
              are illustrated in 'Yellow' color box, which provides example and 
              helps you to understand your original input data."
              , style="padding-left: 0em"),
            
            p("Therefore, to conduct DE analysis, you just need to follow the tabs on the left black dashboard panel 
              from uploading your input data to the end of DE analysis results comparison, 
              this four steps analysis workflow is summarized and illustrated as below."
              , style="padding-left: 0em"),
            
            p("Step 1: Uploading your input data ('Raw Count Table' and 'Meta-data Table')
              via sub-tabs under 'Data Input' tab on the black dashboard panel 
              for single-factor and multi-factor experiment respectively, 
              and visualizing the summarized results of your input data."
              , style="padding-left: 2em"),
            p("Step 2: Filtering out the low expression tags 
              via 'Data Summarization' tab on the black dashboard panel, 
              and visualizing summarized count results after filtering. 
              The low expression filtering criteria (CPM cutoff and sample number) 
              can be set up in the collapsible 'Blue' color box."
              , style="padding-left: 2em"),
            p("Step 3: Conducting DE analysis via sub-tabs under the 'Data Analysis' tab on the black dashboard panel 
              corresponding to 3 different analysis methods. The same as low expression cut off options, 
              the DE analysis options can be set up in the collapsible 'Blue' color box."
              , style="padding-left: 2em"),
            p("Step 4: Comparing/cross-validating DE analysis results 
              via 'Methods Comparison' tab on the black dashboard panel. 
              There are up to 3 different DE analysis methods to choose 
              for comparison in the collapsible 'Blue' color box together with other comparison options."
              , style="padding-left: 2em"),
            #includeHTML(path=paste(getwd(),"www/Introduction.html",sep="/"))
            #imageOutput(paste(getwd(),"www/Introduction.pdf",sep="/"))
            #p("The analysis workflow is summarized as below:"
            #  , style="padding-left: 0em, padding-bottom: 4em"),
          
            img(src="analysis-workflow.png", width=800, style="display: block; margin-left: auto; margin-right: auto;"),
            p("The analysis execution workflow with this App is illustrated in a pdf file, which can be downloaded from  ",
              a("here.", href=as.character(paste("file://~",getwd(),"/www/workflow.pdf", sep="")))
              , style="padding-left: 0em; padding-top: 2em")
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
                        fileInput(inputId="countFile", 
                                  label=h4("Raw Count Table"), 
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
                        )
                    ),
                    ##Meta-data input box under Data Input tab panel
                    box(title = "Input 2: Meta-data Table",
                        solidHeader = T, status = "info",
                        width = 6,
                        fileInput(inputId="metaTab", 
                                  label=h4("Metadata Table"),
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
                        br(),
                        actionButton("dataSubmit", label = "Submit")
                        #submitButton("dataSubmit")
                    )
                  ),
                  
                  helpText(HTML("<div style=\" color: black; font-size: 20px;\">Note: Please make sure the format of you input data are the same as examples of input data shown as below (collapsible by clicking top right corner '+' sign).</div>")),
                  
                  fluidRow(
                    box(title = "Example of input 1",
                        solidHeader = T, status = "warning",
                        collapsible = T, collapsed = T,
                        width = 6,
                        tableOutput("countTabSamp1"),
                        tags$style("#countTabSamp1 table {border: 1px solid black; align: center; margin:auto; margin-bottom: 1em}","#countTabSamp1 th {border: 1px solid black;}","#countTabSamp1 td {border: 1px solid black;}"),
                        helpText(HTML("<div style=\"color: black \">Example of the input data1 - 'raw count data' for the single-factor experiment, 
                                      the 1st column corresponds to the tag ID, such as gene ID, transcript ID, or miRNA ID; 
                                      and 1st row corresponds to the sample ID. 
                                      <br>An example of full 'Raw Count Data' text file for the single-factor experiment is available "), 
                                 a(HTML("here</div>"), href=as.character(paste("file://~",getwd(),"/data/TestData-featureCount.txt", sep="")))
                        )
                        ),
                    box(title = "Example of input 2",
                        solidHeader = T, status = "warning",
                        collapsible = T, collapsed = T,
                        width = 6,
                        fluidRow(
                          tableOutput("metaTabSamp1"),
                          tags$style("#metaTabSamp1 table {border: 1px solid black; align: center; margin:auto; margin-bottom: 1em}","#metaTabSamp1 th {border: 1px solid black;}","#metaTabSamp1 td {border: 1px solid black;}", "#metaTabSamp1 caption {font-size: 16px; color: black; text-align: center;}"),
                          #column(6,
                          #       tableOutput("metaTabSamp1"),
                          #       tags$style("#metaTabSamp1 table {border: 1px solid black; align: center; margin:auto; margin-bottom: 1em}","#metaTabSamp1 th {border: 1px solid black;}","#metaTabSamp1 td {border: 1px solid black;}", "#metaTabSamp1 caption {font-size: 16px; color: black; text-align: center;}")
                          
                          #),
                          #column(6,
                          #       tableOutput("multimetaTabSamp1"),
                          #       tags$style("#multimetaTabSamp1 table {border: 1px solid black; align: center; margin:auto; margin-bottom: 1em}","#multimetaTabSamp1 th {border: 1px solid black;}","#multimetaTabSamp1 td {border: 1px solid black;}", "#multimetaTabSamp1 caption {font-size: 16px; color: black; text-align: center;}")
                          #),
                          #column(4,
                          #       tableOutput("multimetaTabSamp2"),
                          #       tags$style("#multimetaTabSamp2 table {border: 1px solid black; align: center; margin:auto; margin-bottom: 1em}","#multimetaTabSamp2 th {border: 1px solid black;}","#multimetaTabSamp2 td {border: 1px solid black;}", "#multimetaTabSamp2 caption {font-size: 16px; color: black; text-align: center;}")
                          #),
                          br(),
                          br(),
                          helpText(HTML("<div style=\"color: black; padding-left: 1em; padding-right: 1em; margin-top: 1cm \"> Example of the meta-data table for single-factor experiment, 
                                        the 1st column corresponds to the sample name, the 2nd columns corresponds to the single experimental factor - 'Group'. 
                                        <br>An example of 'Meta-data Table' text file corresponding to the example of input 1 - 'Raw Count Data' for the single factor experiment is accessible "), 
                                   a(HTML("here</div>"), 
                                     href=as.character(paste("file://~",getwd(),"/data/TestData-featureCount-meta.txt", sep=""))
                                     )
                          )
                          
                          #helpText("Please specify the factor group level based on your experimental design for DE analysis, 
                          #for the single-factor experiment, it can be any 2 levels of the factor level,
                          #for the multi-factor experiment, it can be any 2 levels of the combined factor levels. 
                          #Please provide correct comparison factor levels for DE analysis, which are summarized on the 'Data Input' tab.")
                          )
                          )
                        )
                  
                    )
                    ),
            
            ##Input information summary under Data Input Panel
            fluidRow(
              box(title = "Input Information Summary",
                  solidHeader = T, status = "success",
                  footer = "Note: There must be a biological replicates from the experiment for DE analysis",
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
                        fileInput(inputId="countFileMulti", 
                                  label=h4("Raw Count Table"), 
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
                                     selected='\t'
                        )
                    ),
                    ##Meta-data input box under Data Input tab panel
                    box(title = "Input 2: Meta-data Table",
                        solidHeader = T, status = "info",
                        width = 6,
                        fileInput(inputId="metaTabMulti", 
                                  label=h4("Metadata Table"),
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
                                     selected='\t'
                        ),
                        br(),
                        actionButton("MultiSubmit", label = "Submit")
                        #submitButton("MultiSubmit")
                    )
                    
                  ),
                  
                  helpText(HTML("<div style=\" color: black; font-size: 20px;\">Note: Please make sure the format of you input data are the same as examples of input data shown as below (collapsible by clicking top right corner '+' sign).</div>")),
                  
                  fluidRow(
                    box(title = "Example of input 1",
                        solidHeader = T, status = "warning",
                        collapsible = T, collapsed = T,
                        width = 6,
                        tableOutput("countTabSampMulti1"),
                        tags$style("#countTabSampMulti1 table {border: 1px solid black; align: center; margin:auto; margin-bottom: 1em}","#countTabSampMulti1 th {border: 1px solid black;}","#countTabSampMulti1 td {border: 1px solid black;}"),
                        helpText(HTML("<div style=\"color: black \">Example of the input data1 - 'Raw Count Data' for the multi-factor experiment, 
                                the 1st column corresponds to the tag ID, such as gene ID, transcript ID, or miRNA ID; 
                                and 1st row corresponds to the sample ID. The format of this 'Raw Count Data' for multi-factor experiment is the same as the single-factor experiment.
                                <br>An example of full 'Raw Count Data' csv file for the multi-factor experiment is accessible "), 
                                 a(HTML("here</div>"), href=as.character(paste("file://~",getwd(),"/data/ReadCounts-Chen-edgeRSpringer-multiFactor.csv", sep="")))
                        )
                    ),
                    
                    box(title = "Example of input 2",
                        solidHeader = T, status = "warning",
                        collapsible = T, collapsed = T,
                        width = 6,
                        fluidRow(tableOutput("multimetaTabSamp1"),
                                 tags$style("#multimetaTabSamp1 table {border: 1px solid black; align: center; margin:auto; margin-bottom: 1em}","#multimetaTabSamp1 th {border: 1px solid black;}","#multimetaTabSamp1 td {border: 1px solid black;}", "#multimetaTabSamp1 caption {font-size: 16px; color: black; text-align: center;}"),
                                 br(),
                                 br(),
                                 helpText(HTML("<div style=\"color: black; padding-left: 1em; padding-right: 1em; margin-top: 1cm \"> Example of the meta-data table for multi-factor experiment, 
                                        the 1st column corresponds to the sample name, the rest columns correspond to the different experimental factors, which are listed in columns one by one.
                                        <br>An example of 'Meta-data Table' csv file corresponding to the example of input 1 - 'Raw Count Data' for the multi-factor experiment is accessible ")
                                          , a(HTML("here</div>"), href=as.character(paste("file://~",getwd(),"/data/ReadCounts-Chen-edgeRSpringer-multiFactor-meta.csv", sep="")))
                                 )
                        )
                    )
                  )
              )
            ),            
            
            fluidRow(
              box(title = "Input Information Summary",
                  solidHeader = T, status = "success",
                  footer = "Note: There must be a biological replicates from the experiment for DE analysis",
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
              box(title = "Low Expression Tags/Genes Removal",
                  solidHeader = T, status = "info",
                  collapsible = T, collapsed = F,
                  width = 12,
                  footer = "Note: Contents in the left dark yellow window are calculated from original raw count data, and contents in 
                  the right green window are results calculated from removing the low expression tags/genes based on the filtering 
                  options specified in the middle blue window, where CPM is count per million, it corresponds 
                  to the read count of (CPM * Library sizes * Normalization factors * 1e-6). 
                  For detailed information regarding filtering low expression tag, 
                  please check the 'Introduction' tab.",
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
                        h6("Low expression tag/gene will be removed with"),
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
                        actionButton(inputId="rmlow", label="Submit")
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
                  )
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
                                           value="control")
                          ),
                          column(6,
                                 textInput(inputId="edgercompGroup2", 
                                           label="Level 2", 
                                           value="treated")
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
                                              choices=c("Nominal-p" ='normp',
                                                        "FDR adjusted-p" ='fdrp'
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
                        ),
                        actionButton(inputId="edgerdeAnalysis", label="Submit")
                    )
                  )
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
                         helpText("Download DEGs results", align="right")
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
                                           value="control")
                          ),
                          column(6,
                                 textInput(inputId="voomcompGroup2", 
                                           label="Group 2", 
                                           value="treated")
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
                                              choices=c("Nominal-p" ='normp',
                                                        "FDR adjusted-p" ='fdrp'
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
                        ),
                        actionButton(inputId="voomdeAnalysis", label="Submit")
                    )
                  )
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
                         helpText("Download DEGs results")
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
                                           value="control")
                          ),
                          column(6,
                                 textInput(inputId="deseq2compGroup2", 
                                           label="Group 2", 
                                           value="treated")
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
                                              choices=c("Nominal-p" ='normp',
                                                        "FDR adjusted-p" ='fdrp'
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
                        ),
                        actionButton(inputId="deseq2deAnalysis", label="Submit")
                    )
                  )
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
                         helpText("Download DEGs results")
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
                                           label = "DE analysis method choices",
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
                                           value="control")
                          ),
                          column(6,
                                 textInput(inputId="decompGroup2", 
                                           label="Group 2", 
                                           value="treated")
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
                                              choices=c("Nominal-p" ='normp',
                                                        "FDR adjusted-p" ='fdrp'
                                              ),
                                              selected='fdrp'
                                 )
                          ),
                          column(5,
                                 textInput(inputId="decompfdr", 
                                           label="p-value or FDR adjusted p-value", 
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
                  helpText("Note 1). Please make sure at least 2 methods are selected for DE analysis comparison; 2). DE analysis is conducted between any 2 specified levels out of available group levels.")
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
