## This App is used for RNA-seq DE analysis with different methods together with method comparison
## Developed by Yan Li, last update on June, 2015

header <- dashboardHeader(
  title = "DE analysis App"
)

sidebar <- dashboardSidebar( 
  sidebarMenu(id="menu1",
              menuItem("Introduction", tabName="intro" , icon = icon("user")),
              menuItem("Data Input", tabName="dataInput", icon = icon("file")),
              menuItem("Data Summarization", tabName="dataSummary", icon = icon("table")),
              menuItem("DE analysis", icon = icon("area-chart"), 
                       menuSubItem("edgeR", tabName="edger"),
                       menuSubItem("Limma-Voom", tabName="limmavoom"),
                       menuSubItem("DEseq2", tabName="deseq2")
              ),
              menuItem("Methods Comparison", icon = icon("bar-chart-o"),tabName="decomp")        
  )
)


body <- dashboardBody(
  tabItems(
    #########################################
    ## Introduction tab panel
    tabItem(tabName="intro",
            h2("Introduction"),
            p("This Application is developed in R with Shiny, it is developed to 1). conduct differential 
              expression (DE) analysis with edgeR, limma-voom, and DESeq2 based on 
              the user provided input data (raw count results and experimental design information); 
              2). cross-validate the DE analysis results with these 3 different DE analysis methods. 
              The goal of this App development is that biologists can use this App to conduct DE analysis
              on their own data and further compare the DE analysis results with 3 different methods."
              , style="padding-left: 0em"),
            #####################################################
            ##Data input section
            h3("1. Data input", style="padding-left: 1em"),
            p("The input data of this App is 2 plain text files, 
              they are named as 'Raw Count Data' and 'Meta-data Table'."
              , style="padding-left: 2em"),
            h4("1. Raw Count Data", style="padding-left: 3em"),            
            p("It contains summarized count results for each sample, the part of this input table can be seen as below:"
              , style="padding-left: 5em"),
            tableOutput("countTabSamp"),
            tags$style("#countTabSamp table {border: 1px solid black; align: left; margin-left: 6em}","#countTabSamp th {border: 1px solid black;}","#countTabSamp td {border: 1px solid black;}"),
            p("Where, the 1st column corresponds to the tag ID, such as gene ID, transcript ID, or miRNA ID; 
              and 1st row corresponds to the sample ID. 
              The full 'Raw Count Data' text file of the test data used in the App can be seen ", 
              a("here", href=as.character(paste("file://~",getwd(),"/data/TestData-feature-count-res.txt", sep="")))
              , style="padding-left: 5em"),
            
            h4("2. Meta-data Table", style="padding-left: 3em"),
            p("It contains summarized experimental design information for each sample. 
              This App supports both single-factor and multi-factor experiment's DE analysis, 
              and the 'Meat-data Table' input is used to reflect the experiment design."
              , style="padding-left: 5em"),
            h5("2.1 Single-factor Experiment", style="padding-left: 5em; font-weight: bold"),
            p("If the experiment only has only one single experimental factor, such as 'Group'. 
              The 'Meta-data Table' input to summarize such experimental design can be prepared 
              in the way presented as below."
              , style="padding-left: 5em"),
            tableOutput("metaTabSamp"),
            tags$style("#metaTabSamp table {border: 1px solid black; align: left;margin-left: 6em}","#metaTabSamp th {border: 1px solid black;}","#metaTabSamp td {border: 1px solid black;}"),
            
            p("Where, the 1st column corresponds to the sample name, and the 2nd column corresponds 
              to the single factor information - 'Group' for each sample."
              , style="padding-left: 5em"),
            
            h5("2.2 Multi-factor Experiment", style="padding-left: 5em; font-weight: bold"),
            #tableOutput("multimetaTabSamp"),
            #tags$style("#multimetaTabSamp table {border: 1px solid black; align: left;margin-left: 6em}","#multimetaTabSamp th {border: 1px solid black;}","#multimetaTabSamp td {border: 1px solid black;}"),
            p("If it is a multi-factor experiment, and there are more than one experimental factors used in the experiment, 
              such multi-factor information can be listed in the following columns one by one."   
              , style="padding-left: 5em"),
            p("Such as in addition to the 'Group' factor, 
              I have two more experimental factors 'Time' and 'Drug' used in the experiment, 
              these factor information for each sample can be listed to 
              the 3rd and 4th column in the 'Meta-data table' 
              after the factor 'Group' in the 2nd column as below."   
              , style="padding-left: 5em"),
            tableOutput("multimetaTabSamp22"),
            tags$style("#multimetaTabSamp22 table {border: 1px solid black; align: left;margin-left: 6em}","#multimetaTabSamp22 th {border: 1px solid black;}","#multimetaTabSamp22 td {border: 1px solid black;}"),
            
            p("The 'Meta-data Table' text file of the test data used in this App is available "
              , a("here", href=as.character(paste("file://~",getwd(),"/data/metatable.txt", sep=""))) 
              , style="padding-left: 5em"),
            #####################################################
            ##Filter low expression tags section
            h3("2. Filter low expression tags", style="padding-left: 1em"),
            p("This step is implemented in the 'Data Summarization' tab, aiming to filter out tags/genes with
              very low counts. The guideline for this step is to keep tags which are expressed in 
              at least one sample out of each factor levels."
              , style="padding-left: 2em"),
            p("Such as there are 3 factor levels in the experiment, 
              then at least 3 samples should better to have expression level 
              above the expression cutoff in count per million (CPM)."
              , style="padding-left: 2em"),
            withMathJax(),
            p("The expression cutoff (CPM) value is determined based on the library size 
              and normalization factors with formula $$\\text{CPM} = \\frac{\\text{counts}}{\\text{library size} * \\text{normalization factors} * 1e^{-6}}$$ 
              Such as the expression cutoff CPM is set to be 10, 
              and the library size and normalization factors are estimated approximately equal to 2e6 and 1 for each sample, 
              then 10 CPM expression cutoff corresponds to about 20 read counts. 
              Therefore, in this example tags/genes with more than 3 samples having less than 
              20 read counts (10 CPM) is classified as low expression tags and will be removed for further DE analysis."
              , style="padding-left: 2em"),
            #####################################################
            ##DE analysis methods section
            h3("3. DE analysis methods", style="padding-left: 1em"),
            p("This App implements 3 different methods to conduct DE analysis, 
              they are edgeR, limma-voom, and DEseq2. 
              The analysis results with respect to different methods can be seen 
              under the 'DE analysis' tab button. 
              For the single-factor experiment, the DE analysis is conducted between any 
              2 levels of that single-factor; 
              for the multi-factor experiment, the DE analysis is conducted in a way 
              to combine all the experimental factors into one combined factor, 
              so that DE analysis can be conducted in any 2 chosen combined factor levels."
              , style="padding-left: 2em"),
            #h4("1. edgeR", style="padding-left: 3em"),
            #p("The multi-factor experiment DE analysis is conducted in a way to combine all the experimental 
            #  factors into tone combined factor, followed with the DE analysis in any chosen combined factor levels."
            #  , style="padding-left: 4em"),
            #h4("2. Limma-Voom", style="padding-left: 3em"),
            #h4("3. DESeq2", style="padding-left: 3em"),
            
            #p("Due to the algorithm differences of different DE analysis methods, the method comparison is 
            #      conducted based on the single-factor experimental design and multi-factor experimental design.
            #      Resutls of method comparison can be seen with the 'Methods Comparison' tab button"),
            
            ############################################
            ##Workflow section
            h2("Workflow"),
            p("To distinguish the user input, analysis results, example and pre-processing results, 
              various color box is used in this App to clarify them."
              , style="padding-left: 0em"),
            p("In general, the analysis results are presented in the 'Green Box', 
              the user input is provided in the 'Blue Box' which is usually collapsible by the right top color '+/-' sign, 
              the example and pre-processing results is usually summarized in the 'Yellow Box', 
              which helps you to understand the input data structure and set up the CPM expression cutoff. 
              Therefore, you just need to follow the tabs in the left black dashboard panel 
              from uploading your own input data to the end DE analysis results comparison, 
              these four steps usage is also summarized and illustrated as below."
              , style="padding-left: 0em"),
            p("Step 1: Uploading your own input data ('Raw Count Table' and 'Meta-data Table')
              via 'Data Input' Tab on the black dashboard panel, and visualizing the input data summary."
              , style="padding-left: 2em"),
            p("Step 2: Filtering out the low expression tags 
              via 'Data Summarization' Tab on the black dashboard panel, 
              and visualizing the count results summary after filtering. 
              The low expression filtering criteria (CPM cutoff and sample number) 
              can be set up in the collapsible blue box."
              , style="padding-left: 2em"),
            p("Step 3: Conducting DE analysis via 'Data Analysis' Tab on the black dashboard panel 
              with 3 different methods based on your own filtering options (in the collapsible blue box)."
              , style="padding-left: 2em"),
            p("Step 4: Comparing/cross-validating DE analysis results 
              via 'Methods Comparison' Tab on the black dashboard panel. 
              There are up to 3 different DE analysis methods to chose (in the collapsible blue box) 
              for cross-validating."
              , style="padding-left: 2em"),
            #includeHTML(path=paste(getwd(),"www/Introduction.html",sep="/"))
            #imageOutput(paste(getwd(),"www/Introduction.pdf",sep="/"))
            img(src="workflow.png", width=1000),
            p("The above workflow can also be downloaded ",
              a("here.", href=as.character(paste("file://~",getwd(),"/www/workflow.pdf", sep="")))
              , style="padding-left: 0em"),
            
            ############################################
            ##Q&A section
            
            h2("Q&A"),
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
              , style="padding-left: 4em"),
            
            ############################################
            ##Feedback section
            h2("Feedback"),
            p("This App is developed and will be maintained by Yan Li at the bioinformatics core, Center for Research Informatics (CRI), 
              Biological Science Division (BSD), University of Chicago."),
            p("As a bioinformatics core, we are actively improving and expanding our NGS analysis services and products.
              If you have any questions, comments, or suggestions, feel free to contact our core at bioinformatics@bsd.uchicago.edu or the developer at yli22@bsd.uchicago.edu"),
            br(),
            img(src="cri.png")
            
            ),
    
    ## End introduction tab panel
    #########################################    
    
    #########################################
    ## First tab content for data input
    tabItem(tabName="dataInput",
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
                    )
                  ),
                  
                  helpText(HTML("<div style=\" color: black; font-size: 20px;\">Note: Please make sure the format of you input data are the same as examples of input data shown as below.</div>")),
                  
                  fluidRow(
                    box(title = "Example of input 1",
                        solidHeader = T, status = "warning",
                        width = 6,
                        tableOutput("countTabSamp1"),
                        tags$style("#countTabSamp1 table {border: 1px solid black; align: center; margin:auto; margin-bottom: 1em}","#countTabSamp1 th {border: 1px solid black;}","#countTabSamp1 td {border: 1px solid black;}"),
                        helpText(HTML("<div style=\"color: black \">The example of the input data1 - 'raw count data', 
                                      the 1st column corresponds to the tag ID, such as gene ID, transcript ID, or miRNA ID; 
                                      and 1st row corresponds to the sample ID. 
                                      The full 'Raw Count Data' text file can be seen "), 
                                 a(HTML("here</div>"), href=as.character(paste("file://~",getwd(),"/data/TestData-feature-count-res.txt", sep="")))
                        )
                        ),
                    box(title = "Example of input 2",
                        solidHeader = T, status = "warning",
                        width = 6,
                        fluidRow(
                          column(6,
                                 tableOutput("metaTabSamp1"),
                                 tags$style("#metaTabSamp1 table {border: 1px solid black; align: center; margin:auto; margin-bottom: 1em}","#metaTabSamp1 th {border: 1px solid black;}","#metaTabSamp1 td {border: 1px solid black;}", "#metaTabSamp1 caption {font-size: 16px; color: black; text-align: center;}")
                          ),
                          column(6,
                                 tableOutput("multimetaTabSamp1"),
                                 tags$style("#multimetaTabSamp1 table {border: 1px solid black; align: center; margin:auto; margin-bottom: 1em}","#multimetaTabSamp1 th {border: 1px solid black;}","#multimetaTabSamp1 td {border: 1px solid black;}", "#multimetaTabSamp1 caption {font-size: 16px; color: black; text-align: center;}")
                          ),
                          #column(4,
                          #       tableOutput("multimetaTabSamp2"),
                          #       tags$style("#multimetaTabSamp2 table {border: 1px solid black; align: center; margin:auto; margin-bottom: 1em}","#multimetaTabSamp2 th {border: 1px solid black;}","#multimetaTabSamp2 td {border: 1px solid black;}", "#multimetaTabSamp2 caption {font-size: 16px; color: black; text-align: center;}")
                          #),
                          br(),
                          br(),
                          helpText(HTML("<div style=\"color: black; padding-left: 1em; padding-right: 1em; margin-top: 1cm \"> Example of the meta-data table for single-factor experiment (left) and multi-factor experiment (right), 
                                        the 1st column corresponds to the sample name, the rest columns correspond to the experimental factors. 
                                        <br>For single-factor experiment (left table), the experimental factor - 'Group' is listed in the 2nd column.
                                        <br>For multi-factor experiment (right table), experimental factors are listed in columns one by one.
                                        <br>Since this App is presenting a single-factor experiment, the 'Meta-data Table' is as the single-factor experimental meta-data table on the left, and the full text file of this test data can be seen"), a(HTML("here</div>"), 
                                                                                                                                                                                                                                                       href=as.character(paste("file://~",getwd(),"/data/metatable.txt", sep=""))
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
                        footer = "DE analysis is conducted between any 2 specified levels out of available group levels",
                        width = 4,
                        solidHeader = T,
                        status = "info",
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
                        ),
                        verbatimTextOutput("edgerGroupLevel")
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
                        footer = "DE analysis is conducted between any 2 specified levels out of available group levels",
                        width = 4,
                        solidHeader = T,
                        status = "info",
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
                        ),
                        verbatimTextOutput("voomGroupLevel")
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
                        footer = "DE analysis is conducted between any 2 specified levels out of available group levels",
                        width = 4,
                        solidHeader = T,
                        status = "info",
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
                        ),
                        verbatimTextOutput("deseq2GroupLevel")
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
                        ),
                        verbatimTextOutput("compGroupLevel")
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
