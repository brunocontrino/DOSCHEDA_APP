# https://bruniec.shinyapps.io/doscheda/

library(hexbin)
library(MASS)
library(shiny)
library(shinydashboard)
library(stringr)
library(affy)
library(limma)
library(DT)
library(ggplot2)
library(vsn)
library(gridExtra)
library(lattice)
library(corrgram)
library(calibrate)
library(reshape2)
library(readxl)
library(MASS)
library(lazyeval)
library(drc)
library(httr)
library(jsonlite)
library(rmarkdown)
library(dplyr)
library(d3heatmap)
library(prodlim)


################### Load Function

loadingLogo <- function( src, loadingsrc, height = NULL, width = NULL, alt = NULL) {
  tagList(
    tags$head(
      tags$script(
        "setInterval(function(){
        if ($('html').attr('class')=='shiny-busy') {
        $('div.busy').show();
        $('div.notbusy').hide();
        } else {
        $('div.busy').hide();
        $('div.notbusy').show();
        }
},100)")
  ),
  tags$a(
    div(class = "busy",  
        img(src=loadingsrc,height = height, width = width, alt = alt)),
    div(class = 'notbusy',
        img(src = src, height = height, width = width, alt = alt))
  )
    )
  }


# plotting functions for sigmoidal
shape_for_ggplot_pred<-function(df_ordered,conc,pred.names){
  cols_to_keep_pred<-c(pred.names,"GeneID","Accession")
  
  forggplot_pred<-vector(mode = "list",length = length(df_ordered$GeneID))
  
  for(i in 1:length(df_ordered$GeneID)){
    tmp_pred<-df_ordered[,cols_to_keep_pred]
    forggplot_pred[[i]]<-melt(tmp_pred[i,], id = c("GeneID", "Accession"), na.rm = TRUE)
  }
  
  forggplot_pred_1<-do.call(rbind, forggplot_pred)
  forggplot_pred_1<-data.frame(forggplot_pred_1,"x"=conc)
  #return(forggplot_pred_1)
}

shape_for_ggplot_perc<-function(df_ordered,conc,finalNames){
  cols_to_keep_perc<-c(finalNames, "GeneID","Accession")
  
  forggplot_perc<- vector(mode = "list",length = length(df_ordered$GeneID))
  for(i in 1:length(df_ordered$GeneID)){
    tmp_perc<-df_ordered[,cols_to_keep_perc]
    forggplot_perc[[i]]<-melt(tmp_perc[i,], id = c("GeneID", "Accession"), na.rm = TRUE)
  }
  
  forggplot_perc_1<-do.call(rbind, forggplot_perc)
  forggplot_perc_1<-data.frame(forggplot_perc_1,"x"=conc)
  #return(forggplot_perc_1)
}


pie_chart <- function(df, main, labels = NULL, condition = NULL) {
  
  # convert the data into percentages. group by conditional variable if needed
  df <- group_by_(df, .dots = c(condition, main)) %>%
    summarise(counts = n()) %>%
    mutate(perc = counts / sum(counts)) %>%
    arrange(desc(perc)) %>%
    mutate(label_pos = cumsum(perc) - perc / 2,
           perc_text = paste0(round(perc * 100), "%", "\n","(",counts, ")"))
  
  # reorder the category factor levels to order the legend
  df[[main]] <- factor(df[[main]], levels = unique(df[[main]]))
  
  # if labels haven't been specified, use what's already there
  if (is.null(labels)) labels <- as.character(df[[main]])
  
  p <- ggplot(data = df, aes_string(x = factor(1), y = "perc", fill = main)) +
    
    # make stacked bar chart with black border
    geom_bar(stat = "identity", color = "black", width = 1) +
    
    # add the percents to the interior of the chart
    geom_text(aes(x = 1.25, y = label_pos, label = perc_text), size = 4) +
    
    # add the category labels to the chart
    # increase x / play with label strings if labels aren't pretty
    geom_text(aes(x = 1.82, y = label_pos, label = labels), size = 4) +
    
    # convert to polar coordinates
    coord_polar(theta = "y") +
    
    # formatting
    scale_y_continuous(breaks = NULL) +
    scale_fill_discrete(name = "", labels = unique(labels)) +
    theme(text = element_text(size = 12),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank())
  
  # facet wrap if that's happening
  if (!is.null(condition)) p <- p + facet_wrap(condition)
  
  return(p)
}


##upload kinase file 
kinome1<-read.csv(file="KinaseFile.csv",  header=TRUE, sep=",")
kinome<-as.vector(toupper(kinome1$new.Symbol))

logvc<- function(x) sqrt(exp(x*log(2))^2 - 1 )

# Peptide matching function 

peptide.match<- function(dr1,dr2,nchan){
  maxrow <- max(nrow(dr1),nrow(dr2))
  minrow <- min(nrow(dr1),nrow(dr2))
  
  adVal <- maxrow - minrow
  
  if(nrow(dr1) == maxrow){
    dr2$addedVals <- adVal
    big.pep <- dr1
    small.pep <- dr2
  } else {
    dr1$addedVals <- adVal
    big.pep <- dr2
    small.pep <- dr1
  }
  
  
  newframe <- big.pep
  colnames(newframe) <- colnames(small.pep)
  newframe[1:minrow,] <- small.pep
  
  for(i in  1:nchan){
    
    newframe[(minrow+1):maxrow,i]<- mean(small.pep[,i])
  }
  
  
  if (all.equal(dim(big.pep),dim(dr1)) == TRUE){
    
    dr2 <- newframe
  } else {
    dr1 <- newframe
  }
  
  list(dr1 = dr1, dr2 = dr2 )
}

panel.shadeNtext <- function (x, y, corr = NULL, col.regions, ...)
{
  if (is.null(corr))
    corr <- cor(x, y, use = "
                pair")
  ncol <- 14
  pal <- col.regions(ncol)
  col.ind <- as.numeric(cut(corr, breaks = seq(from = -1, to = 1,
                                               length = ncol + 1), include.lowest = TRUE))
  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], col = pal[col.ind],
       border = NA)
  box(col = "lightgray")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- formatC(corr, digits = 2, format = "f")
  cex.cor <- .8/strwidth("-X.xx")
  text(0.5, 0.5, r, cex = cex.cor)
}

#########################################################################
## UI
#########################################################################

ui <- shinyUI(dashboardPage(
  dashboardHeader(title = loadingLogo('titleogo.png','logogifopt.gif',height = 55,width = 165)),
  dashboardSidebar(
    sidebarMenu(id = "menu1",
                
                menuItem(text = "Introduction", tabName = "intro"),
                #menuItem(text = "Step by Step Guide", tabName = 'stbyst'),
                conditionalPanel(
                  condition = "input.menu1 == 'stbyst'",
                  radioButtons(inputId = "stage",label = NULL, choices = 
                                 c("Stage 1: Data Upload" = "stage1",
                                   "Stage 2: Navigating Plots" = "stage2",
                                   "Stage 3: Download Results" = "stage3"))
                ),
                
                menuItem(text = "Data Upload", tabName = "dataload"),
                menuItem(text = "Box and Density Plots", tabName = "box"),
                menuItem(text = "MeanSD Plot", tabName = "meansd"),
                menuItem(text = "Mean vs Difference", tabName = "meandiff"),
                conditionalPanel(
                  condition = "input.menu1 == 'meandiff'",
                  uiOutput("meandiff1")
                ),
                menuItem(text = "Corrgram",tabName = "corrgram"),
                menuItem(text = "Compare Replicates",tabName = "repvsrep"),
                conditionalPanel(
                  condition = "input.menu1 == 'repvsrep'",
                  uiOutput("repvsrep1") ),
                menuItem(text = "PCA", tabName = "pca"),
                menuItem(text = "Heatmap", tabName = "testsigmoid"),
                menuItem(text = "Linear Model", tabName = "volcano"),
                conditionalPanel(
                  condition = "input.menu1 == 'volcano' & input.modtyp == 'lin'",
                  sliderInput(inputId = "pvalsli", label = "Change Pvalue threshold",min = 0,max = 0.1,value = 0.05),
                  sliderInput(inputId = "avthrssli", label = "Change avfold threshold",min = 0,max = 1,value = 0.2)
                ),
                menuItem(text = "Sigmoidal fit",tabName = "sigfit"),
                menuItem(text = "Summary", tabName = "summary"),
                menuItem(text = "Download", tabName = "download")
                # menuItem(text = "Test", tabName = "tst")
                
                
                
    )
  ),
  
  dashboardBody(
    tabItems(
      # tabItem(tabName = 'tst',box( DT::dataTableOutput('testkd'))),
      tabItem(tabName = "testsigmoid", 
              box(width = 12, height = 800,
                  d3heatmapOutput("plot9",height = 750)
              )),
      
      
      tabItem(tabName = "intro",  includeHTML(path ="introduction.html")),
      # The step by step has been removed from the app but can be re instated by including the menuItem at the start of the UI 
      tabItem(tabName = "stbyst", 
              conditionalPanel(condition = 'input.stage == "stage1"',
                               tabBox(id = 'stg1tabs', width = 12,
                                      tabPanel(title = "Step 1",
                                               h4("Go to the tab on the sidebar"),
                                               br(),br()
                                               ,HTML('<p><img src="stage1step1.gif" height = 500 /></p>')
                                      ),
                                      tabPanel(title = "Step 2",
                                               h4("Select your data type, please note that intensities will take a few minutes to run."),
                                               br(),HTML('<p><img src="stage1step2.gif" height = 450 /></p>')
                                               
                                      ),
                                      tabPanel(title = "Step 3",
                                               h4("Select the number of channels and replicates."),
                                               br(), 
                                               HTML('<p><img src="stage1step3.gif" height = 450 /></p>')
                                      ),
                                      tabPanel(title = "Step 4",
                                               h4("Choose the columns from your data set which indicate the channels you would like to analyse. Note that you can search for these column names by typing in the input box."),
                                               br(),h4("Remember to check that the uploaded column names correspond with the correct standardised name by looking at the table. You can change these positions by dragging and dropping the selected column names in the desired order in the input box."),
                                               br(),HTML('<p><img src="stage1step4.gif" height = 450 /></p>')
                                      ),
                                      tabPanel(title = "Step 5",
                                               h4("(optional depending on input) If you would like to apply a sigmoidal fit and have an adequate amount of data, use the radio button to change the fit and input your concentrations, separating them with a comma
                                                  . Please Note that this may take a few minutes to run"),
                                               br(),HTML('<p><img src="doschedasigmoid.gif" height = 450 /></p>')
                                      ),
                                      tabPanel(title = "Step 6",
                                               h4("(optional depending on input)  If you are using intensities, you required to input the peptide quality score and  sequence."),
                                               br(),HTML('<p><img src="doschedaintensity.gif" height = 450 /></p>')),
                                      tabPanel(title = "Step 7",
                                               h4("Change to one of the plot tabs and you will see the loading sign appear in place of the logo, now cycle through all the tabs.")
                                      )
                                      
                               )),
              conditionalPanel(condition = 'input.stage == "stage2"'
                               , 
                               tabBox(id = 'stg2tabs', width = 12,
                                      tabPanel(title = "Stage 2",
                                               h4(" Cycle through plots, this is done by clicking on the sidebar tabs and the tabs within the main panel.")
                                      ))
                               
              ),
              conditionalPanel(condition = 'input.stage == "stage3"'
                               , 
                               tabBox(id = 'stg3tabs', width = 12,
                                      tabPanel(title = "Download Data",
                                               h4("Go to the Download tab and press the Download Data button to download the processed data in csv format. and use the â€˜Download Reportâ€™ to download a html file with all your plots with descriptions with them to save your results from the pipeline.")
                                      ),
                                      tabPanel(title = "Download Report",
                                               h4("Go to the Download tab and press the â€˜Download Reportâ€™ button to download a html file with all your plots with description, allowing you to save your results from the pipeline.")
                                               
                                      )
                                      
                               )
                               
                               
              )),
      # End of the step by step 
      
      tabItem(tabName = "dataload",
              box(width = 4,
                  
                  
                  fluidRow(column(5,radioButtons("datype", "Data Type:",
                                                 c("Intensities" = "intensity",
                                                   "Fold Change" = "FC",
                                                   "Log-Fold Change" = "lFC"))),
                           
                           column(5,radioButtons("filetype", "File Type:",
                                                 c(".csv" = "CSV",
                                                   ".txt" = "TXT",
                                                   ".xlsx" = "XLSX")))
                  ),
                  
                  conditionalPanel(condition = 'input.datype == "intensity" & input.modtyp != "sigmoid"',
                                   radioButtons("dorem", "Do removal:",
                                                c("Yes" = "yes",
                                                  "No" = "no"),selected = 'no')
                                
                  )
                  ,
                  
                  
                  
                  # DATA file upload
                  fileInput('file2', 'Choose file',
                            accept = c(".txt",".csv",".xlsx")),
                  
                  fluidRow(column(5, numericInput("chans","# Channels",value = 4, min = 1)),
                           column(5, numericInput("reps","# Replicates",value = 1,min = 1))
                  ),
                  
                  actionButton("changenames", label = "Change Names"),
                  
                  conditionalPanel(condition = "(input.changenames % 2) == 1",
                                   
                                   uiOutput("setNames")
                                   
                  ),
                  uiOutput("ui_choice"),
                  conditionalPanel(condition = 'input.datype == "intensity"',
                                   
                                   uiOutput("ui_sequence"),
                                   uiOutput("ui_qual")
                                   
                                   
                  ),
                  conditionalPanel(condition = 'input.toacc == true || input.datype == "intensity"',
                                   uiOutput("ui_accession")
                                   
                  ),
                  conditionalPanel(condition = 'input.toacc == true & input.datype != "intensity"',
                                   uiOutput("ui_uniqpep")
                                   
                  ),
                  
                  conditionalPanel(condition = "input.modtyp == 'sigmoid'",
                    checkboxInput('incpd', 'Insert pulldown of pull down', value = FALSE)
                    
                  ),
                  conditionalPanel(condition = 'input.incpd == true',
                                   uiOutput("ui_pdofpd")
                                   
                  ),
                  conditionalPanel( condition =  "input.datype != 'intensity'",
                                    
                                    checkboxInput("toacc", "Data is NOT Proteome Discoverer 2.1", FALSE)
                  ),
                  
                  checkboxInput('genefile', 'Upload an Accession to Gene ID file', value = FALSE)

                  
                  # uiOutput("ui_sequence"),
                  ## remove column obselete 
                  
                  # fluidRow(
                  # column(5,conditionalPanel(
                  #   condition = "input.conf != 0",
                  #   actionButton(inputId = "removechan", label = "Remove Column")
                  # ))),
                  # 
                  # conditionalPanel(
                  #   condition = "(input.removechan % 2) == 1",
                  #   uiOutput("remchan")
                  # )

              ),
              box(width = 8,
                  DT::dataTableOutput("test")
              ),
              box(width = 8,
                  conditionalPanel(condition = 'input.chans >= 5 & input.reps == 1',
                                   radioButtons("modtyp", "Fit model:",
                                                c("Sigmoidal" = "sigmoid",
                                                  "Linear" = "lin"
                                                ),selected = "lin")),
                  conditionalPanel(condition = 'input.chans < 5',
                                   h1("Less than 5 channels, only a linear model can be applied")
                  ),
                  conditionalPanel(condition = "input.modtyp == 'sigmoid'",
                                   textInput("concsig", "Enter vector of concentrations from low to high [comma delimited]. Ensure these are not log concentrations"),
                                   textOutput("concsig")
                                   
                                   
                  )
              ),
              
              box(width = 8,
                  selectizeInput( 'organism',label = 'Select your organism:', choice = c('H.sapiens','R.norvegicus','Canis familiaris','Mus musculus','C. elegans', 'S.cerevisiae'), selected =  'H.sapiens'),
                  conditionalPanel(condition = 'input.genefile == true',
                                   
                                   fluidRow(column(5,fileInput('geneF', 'Choose you Accession to Gene ID file',
                                                               accept = c(".txt",".csv",".xlsx"))),
                                            column(5,radioButtons("generadio", "File Type:",
                                                                  c(".csv" = "CSV",
                                                                    ".txt" = "TXT",
                                                                    ".xlsx" = "XLSX"))
                                                   
                                            ))         
                  )
                  )
      ),
      
      tabItem(tabName = "box",
              
              fluidRow(
                tabBox(id = "tb1", width = 9, height = 700,
                       tabPanel("Box",plotOutput(height = 700, "bar")),
                       tabPanel("Density", plotOutput(height = 700, "plot2")),
                       tabPanel("Venn", plotOutput(height = 700, "venn"))
                       
                       
                ),
                box(width = 3,
                    conditionalPanel(condition = "input.tb1 == 'Venn'",
                                     
                                     checkboxInput(inputId = 'venninput' ,label = 'Include file', value = FALSE)
                    ),
                    conditionalPanel('input.venninput == true',
                                     radioButtons("filetype2", "File Type:",
                                                  c(".csv" = "CSV",
                                                    ".txt" = "TXT",
                                                    ".xlsx" = "XLSX")),
                                     
                                     fileInput(inputId = 'venninp',label = 'Choose input file',
                                               accept = c(".txt",".csv",".xlsx")                                                )
                    )
                )
              )
      ),
      
      tabItem(tabName = "meansd",
              fluidRow(box(width = 9, height = 700, plotOutput(height = 600, "plot3")))),
      tabItem(tabName = "meandiff",
              fluidRow(
                box(width = 9, height = 700, plotOutput(height = 600, "plot6"))
              )),
      tabItem("corrgram", fluidRow(box(width = 9, height = 800, plotOutput(height = 700,"plot5")))
      ),
      tabItem("repvsrep", fluidRow(box(width = 9, height = 800, plotOutput(height = 700,"repvsrep")))),
      
      tabItem(tabName = "pca", fluidRow( box(width = 9, height = 800, plotOutput(height = 700,"plot7")))),
      
      tabItem(tabName = "volcano", fluidRow( tabBox(id = "volplots", width = 9, height = 800,
                                                    tabPanel(title = "P-Value Distribution",plotOutput(height = 700,"plot4")),
                                                    tabPanel(title = "Slope", plotOutput(height = 700,"plot8")),
                                                    tabPanel(title = "Intercept", plotOutput(height = 700,"volcanoint")),
                                                    tabPanel(title = "Quadratic", plotOutput(height = 700,"volcanoquad"))
      ))),
      tabItem(tabName = "sigfit", fluidRow( tabBox(id = "sigplots", width = 10, height = 800,
                                                   tabPanel(title = "Difference Top-Bottom",plotOutput(height = 700,"DiffTopBottom")),
                                                   tabPanel(title = "RB50", plotOutput(height = 700,"RB50")),
                                                   tabPanel(title = "Slope", plotOutput(height = 700,"Slope_pl"))
                                                 
      ))),
      
      
      tabItem(tabName = "summary", fluidRow(
        tabBox(id = "sumtab", width = 8,
               tabPanel("Data p-values",div(DT::dataTableOutput("testmerge",width = '100%'),style = "font-size:90%")),
               tabPanel("Kinases", DT::dataTableOutput("kintab"))
        ),box(title = "Corrgram QC", width = 4,
              infoBoxOutput(outputId = "corrinfo",width = 12)
              
        ),box(title = "P Value QC", width = 4,
              
              conditionalPanel(condition = "input.modtyp == 'sigmoid'",
                               infoBoxOutput("siginfodt",width = 12),
                               infoBoxOutput("siginfoslop",width = 12),
                               infoBoxOutput("siginfodiff",width = 12)
                               ),
              
              conditionalPanel(condition = "input.modtyp != 'sigmoid'",
                               infoBoxOutput(outputId = "infopvalslo",width = 12),
                               br(),br(),br(), br(),br(), br(),
                               
                               infoBoxOutput(outputId = "infopvalint",width = 12),
                               br(),br(),br(),  br(),br(), br(),
                               infoBoxOutput(outputId = "infopvalquad",width = 12)
              )
        )
      )),
      tabItem(tabName = "download",
              box(width = 12,    
                  textInput("dataset", "Filename for Data", value = "Doscheda"),
                  
                  downloadButton('downloadData', 'Download Data'),
                  
                  downloadButton("report", "Generate report")
                  
                  
                  
              ),
              
              box(width = 12,
                  conditionalPanel(condition = "input.datype == 'intensity'",
                                   actionButton(inputId = 'intcalc',"Calculate Removed Peptides")
                                   
                  ), br(),
                  
                  conditionalPanel(condition = "input.intcalc != 0",
                                   
                                   downloadButton('peprmv', "Download Removed Peptides")
                                   
                  )
               )
              
      )
      
      
    )
  )
  
))

#########################################################
## SERVER
#########################################################

server <- shinyServer(function(input, output) {
  
  options(shiny.maxRequestSize=100*1024^2)
  ##############################################
  
  ## uniprot gene load in any changes to uni prot file will affect this section 
  
  uniprotGene <- reactive({
    
    if(input$organism == 'H.sapiens'){
      uniprot <- read.delim(file = 'uniprot/HumanSwissP_TrEMBL_Jan2017.tab',stringsAsFactors = FALSE)
    } else if(input$organism == 'Canis familiaris'){
      
      uniprot <- read.delim(file = 'uniprot/DogSwissP_TrEMBL_Jan2017.tab',stringsAsFactors = FALSE)
      
    }else if(input$organism == 'Mus musculus'){
      uniprot <- read.delim(file = 'uniprot/MouseSwissP_TrEMBL_Jan2017.tab',stringsAsFactors = FALSE)
      
      
    }else if (input$organism == 'Rattus norvegicus'){
      
      uniprot <- read.delim(file = 'uniprot/RatSwissP_TrEMBL_Jan2017.tab',stringsAsFactors = FALSE)
      uniprot[,2] <- gsub('_RAT','',uniprot[,2])
      colnames(uniprot)[2] <- 'Gene.names'
    }else{
      uniprot <- read.delim(file = 'uniprot/YeastSwissP_TrEMBL_Jan2017.tab',stringsAsFactors = FALSE)
      
    }
    
    uniprot[,1:2]
  })
  
  ## upload gene file 
  
  uploadGene <- reactive({
    
    inFile <- input$geneF
    
    if(is.null(inFile)){
      return(NULL)
    }
    if( input$generadio == "TXT"){
      file.rename(inFile$datapath,
                  paste(inFile$datapath, ".txt", sep=""))
      # read.xlsx(paste(inFile$datapath, ".xlsx", sep=""),
      #           sheetIndex = 1, header = TRUE, stringsAsFactors = FALSE)
      conv<- read.delim(paste(inFile$datapath, ".txt", sep=""),
                        header = TRUE,stringsAsFactors = FALSE)
      
    } else if (input$generadio == "XLSX"){
      file.rename(inFile$datapath,
                  paste(inFile$datapath, ".xlsx", sep=""))
      # paste(inFile$datapath, ".xlsx", sep=""))
      # read.xlsx(paste(inFile$datapath, ".xlsx", sep=""),
      #           sheetIndex = 1, header = TRUE, stringsAsFactors = FALSE)
      conv<- read_excel(paste(inFile$datapath, ".xlsx", sep=""),
                        sheet = 1)
      
    } else {
      file.rename(inFile$datapath,
                  paste(inFile$datapath, ".csv", sep=""))
      # paste(inFile$datapath, ".xlsx", sep=""))
      # read.xlsx(paste(inFile$datapath, ".xlsx", sep=""),
      #           sheetIndex = 1, header = TRUE, stringsAsFactors = FALSE)
      conv<- read.csv(paste(inFile$datapath, ".csv", sep=""), header = TRUE, stringsAsFactors = FALSE)
      conv<- conv[,-1]
    }
    
    conv <- conv[,1:2]
    colnames(conv) <- c('Entry', 'Gene.names')
    conv
  })
  
## venn diagram intersect file upload
  
uploadVenn <- reactive({

    if(input$venninput == FALSE){
      return(NULL)
    } else{ 
    inFile <- input$venninp
    
    if(is.null(inFile)){
      return(NULL)
    }
    if( input$filetype2 == "TXT"){
      file.rename(inFile$datapath,
                  paste(inFile$datapath, ".txt", sep=""))
      # read.xlsx(paste(inFile$datapath, ".xlsx", sep=""),
      #           sheetIndex = 1, header = TRUE, stringsAsFactors = FALSE)
      conv<- read.delim(paste(inFile$datapath, ".txt", sep=""),
                        header = TRUE,stringsAsFactors = FALSE)
      
    } else if (input$filetype2 == "XLSX"){
      file.rename(inFile$datapath,
                  paste(inFile$datapath, ".xlsx", sep=""))
      # paste(inFile$datapath, ".xlsx", sep=""))
      # read.xlsx(paste(inFile$datapath, ".xlsx", sep=""),
      #           sheetIndex = 1, header = TRUE, stringsAsFactors = FALSE)
      conv<- read_excel(paste(inFile$datapath, ".xlsx", sep=""),
                        sheet = 1)
      
    } else {
      file.rename(inFile$datapath,
                  paste(inFile$datapath, ".csv", sep=""))
      # paste(inFile$datapath, ".xlsx", sep=""))
      # read.xlsx(paste(inFile$datapath, ".xlsx", sep=""),
      #           sheetIndex = 1, header = TRUE, stringsAsFactors = FALSE)
      conv<- read.csv(paste(inFile$datapath, ".csv", sep=""), header = TRUE, stringsAsFactors = FALSE)
      conv<- conv[,-1]
    }
    
    conv
    }
  })
  
  ### crapome file upload
  
  
  
  crapome<- reactive({
    if(input$organism == 'S.cerevisiae'){
      a <- read.csv('CRAPOME_YEAST',stringsAsFactors = FALSE)
    } else{
      a <- read.csv('CRAPOME_HUMAN',stringsAsFactors = FALSE)
    }
    
    print(head(a))
    a
  })
  
  ### DATA file upload
  
    data <- reactive({
    
    
    inFile <- input$file2
    
    if(is.null(inFile)){
      return(NULL)
    }
    if( input$filetype == "TXT"){
      file.rename(inFile$datapath,
                  paste(inFile$datapath, ".txt", sep=""))
      # read.xlsx(paste(inFile$datapath, ".xlsx", sep=""),
      #           sheetIndex = 1, header = TRUE, stringsAsFactors = FALSE)
      conv<- read.delim(paste(inFile$datapath, ".txt", sep=""),
                        header = TRUE,stringsAsFactors = FALSE)
      
    } else if (input$filetype == "XLSX"){
      file.rename(inFile$datapath,
                  paste(inFile$datapath, ".xlsx", sep=""))
      # paste(inFile$datapath, ".xlsx", sep=""))
      # read.xlsx(paste(inFile$datapath, ".xlsx", sep=""),
      #           sheetIndex = 1, header = TRUE, stringsAsFactors = FALSE)
      conv<- read_excel(paste(inFile$datapath, ".xlsx", sep=""),
                        sheet = 1)
      
    } else {
      file.rename(inFile$datapath,
                  paste(inFile$datapath, ".csv", sep=""))
      # paste(inFile$datapath, ".xlsx", sep=""))
      # read.xlsx(paste(inFile$datapath, ".xlsx", sep=""),
      #           sheetIndex = 1, header = TRUE, stringsAsFactors = FALSE)
      conv<- read.csv(paste(inFile$datapath, ".csv", sep=""), header = TRUE, stringsAsFactors = FALSE)
    }
    
    
    data.frame(conv)
  })
  

    ## ui elements 
  output$ui_accession<- renderUI({
    
    temp <- data()
    
    ## using selectizeInput with drag_drop and DT
    selectizeInput("accession", "Choose Accession", choices  = colnames(temp) ,
                   selected = NULL,
                   multiple = FALSE,
                   options = list(plugins = list('remove_button')))
  })  
  
  output$ui_uniqpep<- renderUI({
    
    temp <- data()
    
    ## using selectizeInput with drag_drop and DT
    selectizeInput("unipeps", "Choose Unique Peptides", choices  = colnames(temp) ,
                   selected = NULL,
                   multiple = FALSE,
                   options = list(plugins = list('remove_button')))
  }) 
  
  output$ui_pdofpd<- renderUI({
    
    temp <- data()
    
    ## using selectizeInput with drag_drop and DT
    selectizeInput("pdofpd", "Choose Pulldown of pulldown", choices  = colnames(temp) ,
                   selected = NULL,
                   multiple = FALSE,
                   options = list(plugins = list('remove_button')))
    
  })  
  
  
  
  output$ui_sequence<- renderUI({
    
    temp <- data()
    
    ## using selectizeInput with drag_drop and DT
    selectizeInput("sequence", "Choose Sequence", choices  = colnames(temp) ,
                   selected = NULL,
                   multiple = FALSE,
                   options = list(plugins = list('remove_button')))
  })  
  
  output$ui_choice<- renderUI({
    reqtemp <- data()
    req(reqtemp)
    ## using selectizeInput with drag_drop and DT
    selectizeInput("view_vars", "Select variables to show:", choices  = colNamTot(),
                   selected = colNamSel(),
                   multiple = TRUE,
                   options = list(plugins = list('remove_button', 'drag_drop'), maxItems = nchans()))
  })
  
  
  
  output$ui_qual<- renderUI({
    reqtemp <- data()
    req(reqtemp)
    ## using selectizeInput with drag_drop and DT
    selectizeInput("quality", "Select Peptide Qvality Score:", choices  = colNamTot(),
                   selected = colNamSel(),
                   multiple = FALSE,
                   options = list(plugins = list('remove_button', 'drag_drop'), maxItems = 1))
  })
  
  output$remchan<- renderUI({
    reqtemp <- data()
    req(reqtemp)
    
    selectizeInput(
      'remchan', 'Choose Channels to remove', choices = finalNames(), multiple = TRUE,
      options = list(plugins = list('remove_button'))
    )
    
  })
  
  ## final names setting 
  standardNames <- reactive({
    if(input$datype == "intensity"){
      
      out<-   paste("rep",rep(1:input$reps,each = input$chans),"_",
                    rep(paste("C",c("ontrol",0:(input$chans - 2)), sep = ""),input$reps),sep = "")
      
      vals$labelNames = setNames(out, out)
      
    }else{
      out<-  paste("rep",rep(1:input$reps,each = input$chans),"_",
                   rep(paste("C",0:(input$chans - 1), sep = ""),input$reps),sep = "")
      vals$labelNames = setNames(out, out)
      
    }
    
    return(out)
  })
  
  output$setNames <- renderUI({
    
    list(
      h4("Change Names")
      , selectInput("nameToChange", "Standard name"
                    , names(vals$labelNames))
      , textInput("labelToAdd", "Change to")
      , actionButton("makeLabel", "Set label")
    )
    
  })
  
  observeEvent(input$makeLabel, {
    vals$labelNames[input$nameToChange] <- input$labelToAdd
  })
  
  finalNames<- reactive({
    

    vals$labelNames
      
    
    
    
  })
  
  
  
  ##############################################
  

  
  sigConc <- reactive({
    
    if(input$modtyp == "sigmoid"){
      sort(as.numeric(unlist(strsplit(input$concsig,","))))
    }else{
      return(NULL)
    }
    
  })
  
  
  output$concsig <- renderText({
    
    
    x <- as.numeric(unlist(strsplit(input$concsig,",")))
    if(input$datype== 'intensity'){
      n<- input$chans - 1
    }else{
      
      n <- input$chans*input$reps
    }
    
    if(length(x) < n){
 
        paste0("Not enough concentrations, please put ",n," concentrations")
      
    } else if(length(x) > n){

        paste0("Too many concentrations, please put ",n," concentrations")
        
      
    }else{
      x<- sort(x)
      print(c("Concentrations:",x)) 
    }
    
    
  })
  
  
  
  
  output$selectnames <- renderUI({
    selectizeInput("selectnames", "Change Names", choices = standardNames(),
                   options = list(maxOptions = input$reps * input$chans))
  })
  
  
  
  
  vals <- reactiveValues(
    labelNames = character()
  )
  
  #### LABEL UPLOAD
  
  
  #### data imported and finalNames is vector of colnames ####
  
  rReps <- reactive({
    
    input$reps
    
  })
  
  nchans <- reactive({
    input$reps * input$chans
    
  })
  
  ### probably dont need this
  ### datachange table
  
  
  ### selectize
  
  colNamTot<- reactive({
    reqtemp <- data()
    req(reqtemp)
    a2 <- data()
    colnames(a2)
  })
  
  colNamSel<- reactive({
    a1 <- data()
    a1 <- colnames(a1)
    sample(a1,nchans())
    
  })
  

  
 
  
  indexmatrix <- reactive({
    reqtemp <- data()
    req(reqtemp)
    if(input$reps <= 1){
      NULL
    } else {
      if(input$datype == "intensity"){
        channe <- channels()
        reps <- input$reps
        chans <- input$chans - 1
        ser <- finalNames()
        ser<- ser[-seq(1,by = input$chans,length.out = input$reps)]
        
      }else{
        ser <- finalNames()
        channe <- channels()
        reps <- input$reps
        chans <- input$chans
        
      }
      
      combinations <- t(combn(reps,2))
      combmat<- matrix(rep(as.vector(t(combinations)),chans),ncol = 2,byrow = TRUE)
      
      # create factor for repeat ...
      
      name.vec <- 1:(chans*reps)
      repfac <- rep(1:chans,times = reps)
      index <- rep(0:(reps-1),each = chans)
      combfac<- rep(1:(reps),each = chans)
      
      # total combinations  = chans * nrow()
      columnindex <- matrix(0 , ncol = 5 , nrow = chans * nrow(combinations) )
      
      colnames(columnindex) <- c("concentration","rep1","rep2","index1","index2")
      columnindex[,1] <- rep(1:chans ,each = nrow(combinations))
      
      columnindex[,2:3]<- combmat
      
      columnindex
      # create matrix which will be indexed by first 3 columns of column index
      
      index.mat <- matrix(name.vec,ncol = reps)
      
      for(i in 1:nrow(columnindex)){
        
        columnindex[i,4] <- index.mat[columnindex[i,1],columnindex[i,2]]
        columnindex[i,5] <- index.mat[columnindex[i,1],columnindex[i,3]]
        
      }
      

      
      
      dadt <- dataMerge()
      
      create.names <- rep("", nrow(columnindex))
      
      
      for( i in 1:nrow(columnindex)){
        create.names[i] <- paste(ser[columnindex[i,4]], "vs", ser[columnindex[i,5]])
      }
      
      final.mat <- data.frame(names = create.names,columnindex)
      final.mat
      final.mat[!is.na(rowSums(matrix(match(columnindex,channe),ncol=5))),]
    }
  })
  
  
  
  selRemCol <- reactive({
    reqtemp <- data()
    req(reqtemp)
    
    a1<- finalNames()
    a1<- as.character(a1)
    
    if(is.null(input$remchan)){
      0
    }else{
      -match(input$remchan,a1)
      
    }
    
  })
  
  channels <- reactive({
    reqtemp <- data()
    req(reqtemp)
    
    if(input$datype == 'intensity'){
      if(input$reps == 1){
        1:(input$chans*input$reps -1)
        
      }else{
        1:(input$chans*input$reps -2)
      }
    }else{
      
      a1 <- selRemCol()
      vec <- 1:nchans()
      if(a1 != 0){
        vec <- vec[selRemCol()]
      }
      vec
    }
  })
  
  
  
  ##### DATA HANDLING #####
  
  ## Intensities
  
  
  intData <- reactive({
    
    if(input$datype == "intensity"){
      
      
      tempdat <- data()
      
      
      # make matrix of descriptions and accessions to filter by common proteins
      
      accDesMat<- as.character(tempdat[as.character(tempdat[,input$accession]) == "",c(input$accession)])
      
      if(input$incpd == TRUE){
        data.merged <- tempdat[,c(input$view_vars, input$accession, input$sequence,input$quality,input$pdofpd)]
        
      }else{
        data.merged <- tempdat[,c(input$view_vars, input$accession, input$sequence,input$quality)]
      }
      
      
      channelnames <- standardNames()
      repindex <- rep(1:input$reps,each = input$chans)
      totfal <- rep(FALSE , (input$chans + 3 ))
      
      if(input$incpd == TRUE){
        
        newdf <- cbind(data.merged[,1:input$chans], data.merged[,c(input$accession, input$sequence,input$quality,input$pdofpd)])
        colnames(newdf) <- c(channelnames[repindex  == 1],"Accession","Sequence", "Quality",'Kd')
        newdf <- newdf[!is.na(rowSums(newdf[,1:input$chans])),]
        newdf <- newdf[newdf$Quality <= 0.05, ]
        newdf <- data.frame(newdf, outliers = rep(0,length(newdf[,1])), uniquePeps = rep(0,length(newdf[,1])), addedVals = rep(0,length(newdf[,1])),Kd = newdf$Kd)
        
        
      }else{
        newdf <- cbind(data.merged[,1:input$chans], data.merged[,c(input$accession, input$sequence,input$quality)])
        colnames(newdf) <- c(channelnames[repindex  == 1],"Accession","Sequence", "Quality")
        newdf <- newdf[!is.na(rowSums(newdf[,1:input$chans])),]
        newdf <- newdf[newdf$Quality <= 0.05, ]
      }  
      
      if(input$reps == 2){
        channelnames <- paste("rep",rep(1:input$reps,each = input$chans),"_",
                              rep(paste("C",c("ontrol",0:(input$chans - 2)), sep = ""),input$reps),sep = "")  
        newdf2 <- cbind(data.merged[,(input$chans+1):(2*input$chans) ],data.merged[,c(input$accession, input$sequence,input$quality)])
        colnames(newdf2) <- c(channelnames[repindex  == 2],"Accession","Sequence","Quality")
        newdf2 <- newdf2[!is.na(rowSums(newdf2[,1:input$chans])),]
        newdf2 <- newdf2[newdf2$Quality <= 0.05, ]
        
        
        
        newdf <- data.frame(newdf, outliers = rep(0,length(newdf[,1])), uniquePeps = rep(0,length(newdf[,1])), addedVals = rep(0,length(newdf[,1])))
        newdf2 <- data.frame(newdf2, outliers = rep(0,length(newdf2[,1])), uniquePeps = rep(0,length(newdf2[,1])), addedVals = rep(0,length(newdf2[,1])))
        
        
        common.proteins <- intersect(unique(newdf$Accession),unique(newdf2$Accession))
        
        newdf <- newdf[!is.na(match(newdf$Accession,common.proteins)),]
        newdf2 <- newdf2[!is.na(match(newdf2$Accession,common.proteins)),]
      }else{
        common.proteins <- unique(newdf$Accession)
      }   
      
      
      
      if(input$dorem == "no"){
      
        if(input$reps == 1 & input$modtyp == 'sigmoid'){
          if(input$incpd == TRUE){
            prot1 <- unique(newdf$Accession)
            sumkd <- rep(0,length(prot1))
            protdf<- newdf[1:length(prot1),colnames(newdf) != "Sequence"]
            for(i in 1:length(prot1)){
              ## grep for total intensity () includes non unique peps
              protdf[i,1:input$chans] <- apply(newdf[grep(prot1[i],newdf$Accession),1:input$chans],2,sum,na.rm = TRUE)
              protdf$Kd[i] <- sum(newdf$Kd[grep(prot1[i],newdf$Accession)],na.rm = TRUE)
              protdf$Accession[i] <- prot1[i]
              
              ## use == to get unique peptides per protein per repeat
              # protdf$uniquePeps[i] <- length(unique(newdf$Sequence[newdf$Accession == prot1[i]]))
              protdf$uniquePeps[i] <- length(unique(unique(newdf$Sequence[newdf$Accession == prot1[i]])))
              
            }
            
            protdf$uniquePeps[grep(";",protdf$Accession)] <- 0
            
            fc1 <- protdf[,1] / protdf[,2:input$chans]
            
            Kd <- protdf$Kd / protdf[,1]
            
            fcprotdf <- data.frame(log2(fc1),
                                   Accession = protdf$Accession,uniquePepr1 = protdf$uniquePeps,
                                   uniquePepr2 =protdf$uniquePeps, Kd = Kd)
            fcprotdf<- fcprotdf[!is.na(rowSums(fcprotdf[1:(input$chans - 1 )])),]
            fcprotdf
            
          }else{
          
            prot1 <- unique(newdf$Accession)
            protdf<- newdf[1:length(prot1),colnames(newdf) != "Sequence"]
            
            for(i in 1:length(prot1)){
              ## grep for total intensity () includes non unique peps
              protdf[i,1:input$chans] <- apply(newdf[grep(prot1[i],newdf$Accession),1:input$chans],2,sum,na.rm = TRUE)
              protdf$Accession[i] <- prot1[i]
              
              ## use == to get unique peptides per protein per repeat
              # protdf$uniquePeps[i] <- length(unique(newdf$Sequence[newdf$Accession == prot1[i]]))
              protdf$uniquePeps[i] <- length(unique(unique(newdf$Sequence[newdf$Accession == prot1[i]])))
              
            }
            
            protdf$uniquePeps[grep(";",protdf$Accession)] <- 0
            
            fc1 <- protdf[,1] / protdf[,2:input$chans]
            
            
            fcprotdf <- data.frame(log2(fc1),
                                   Accession = protdf$Accession,uniquePepr1 = protdf$uniquePeps,
                                   uniquePepr2 =protdf$uniquePeps)
            fcprotdf<- fcprotdf[!is.na(rowSums(fcprotdf[1:(input$chans - 1 )])),]
            fcprotdf
            
          }
          
        
        }else{
        prot1 <- unique(newdf$Accession)
        protdf<- newdf[1:length(prot1),colnames(newdf) != "Sequence"]
        
        
        for(i in 1:length(prot1)){
          ## grep for total intensity () includes non unique peps
          protdf[i,1:input$chans] <- apply(newdf[grep(prot1[i],newdf$Accession),1:input$chans],2,sum,na.rm = TRUE)
          
          protdf$Accession[i] <- prot1[i]
          
          ## use == to get unique peptides per protein per repeat
          # protdf$uniquePeps[i] <- length(unique(newdf$Sequence[newdf$Accession == prot1[i]]))
          protdf$uniquePeps[i] <- length(unique(c(unique(newdf$Sequence[newdf$Accession == prot1[i]]),unique(newdf2$Sequence[newdf2$Accession == prot1[i]]))))
          
        }
        
        protdf$uniquePeps[grep(";",protdf$Accession)] <- 0
        
        
          prot2 <- unique(newdf2$Accession)
          ### second data frame
          
          protdf2<- newdf2[1:length(prot2),colnames(newdf2) != "Sequence"]
          
          for(i in 1:length(prot2)){
          
            ## grep for total intensity () includes non unique peps
            protdf2[i,1:input$chans] <- apply(newdf2[grep(prot2[i],newdf2$Accession),1:input$chans],2,median,na.rm = TRUE)
            protdf2$Accession[i] <- prot2[i]
            
            ## use == to get unique peptides per protein per repeat
            protdf2$uniquePeps[i] <- length(unique(newdf2$Sequence[newdf2$Accession == prot2[i]]))
          }
          
          protdf2$uniquePeps[grep(";",protdf2$Accession)] <- 0
        
        
          com.prot <- intersect(protdf$Accession,protdf2$Accession)
        
          fc1 <- protdf[,1] / protdf[,2:input$chans]
          fc2 <-  protdf2[,1] / protdf2[,2:input$chans]
          
        
          fcprotdf <- data.frame(log2(fc1[match(com.prot,protdf$Accession),]),
                                 log2(fc2[match(com.prot,protdf2$Accession),]),
                                 Accession = com.prot,uniquePepr1 = protdf$uniquePeps[match(com.prot,protdf$Accession)],
                                 uniquePepr2 =protdf2$uniquePeps[match(com.prot,protdf2$Accession)])
          fcprotdf<- fcprotdf[!is.na(rowSums(fcprotdf[1:(input$reps*input$chans -2 )])),]
          fcprotdf
        }
        
        
      }else{
        
        totpepdf <- NULL
        totpepdf2 <- NULL
        
        for (z in 1:length(common.proteins)){
          temp <- newdf[newdf$Accession == common.proteins[z],]
          temp2 <- newdf2[newdf2$Accession == common.proteins[z],]
          # first step: check if all peptides are unique ...
          
          if(all.equal(grep(";",temp$Accession),integer(0)) == TRUE){
            
            uniPeptides1 <- length(unique(temp$Sequence))
            
          }  else {
            uniPeptides1 <- length(unique(temp$Sequence[-grep(";",temp$Accession)]))
          }
          
          if(all.equal(grep(";",temp2$Accession),integer(0)) == TRUE){
            
            uniPeptides2 <- length(unique(temp2$Sequence))
            
          }else{
            
            uniPeptides2 <- length(unique(temp2[-grep(";",temp2$Accession)]))
          }
          
          ## add unique peptide column
          
          temp$uniquePeps <- uniPeptides1
          temp2$uniquePeps <- uniPeptides2
          
          
          tempPep <- intersect(unique(temp$Sequence),unique(temp2$Sequence))
          
          
          if(all.equal(tempPep, character(0)) == TRUE){
            next
          }
          
          for(i in 1:length(tempPep)){
            
            if(sum(temp$Sequence == tempPep[i]) != sum(temp2$Sequence == tempPep[i])){
              
              dr <- peptide.match(temp[temp$Sequence == tempPep[i],],temp2[temp2$Sequence == tempPep[i],],input$chans)
              dr1 <- dr$dr1
              dr2 <- dr$dr2
              
            } else {
              
              dr1 <- temp[temp$Sequence == tempPep[i],]
              dr2 <- temp2[temp2$Sequence == tempPep[i],]
              
              
            }
            
            tempoindex <- rep(FALSE,nrow(dr1))
            for(j in 1:nrow(dr1)){
              percor <- cor.test(log2(as.numeric(dr1[j,1:input$chans])),log2(as.numeric(dr2[j,1:input$chans])))
              tempoindex[j] <- percor$estimate < 0.4
            }
            
            dr1[tempoindex,1:input$chans] <- NA
            dr2[tempoindex,1:input$chans] <- NA
            
            
            
            tempReplace <- dr1
            tempReplace2 <- dr2
            
            
            
            temp <- temp[temp$Sequence != tempPep[i],]
            temp <- rbind(temp,tempReplace)
            
            temp2 <- temp2[temp2$Sequence != tempPep[i],]
            temp2 <- rbind(temp2,tempReplace2)
            
          }
          
          temp<- temp[match(tempPep,temp$Sequence),]
          temp2<- temp2[match(tempPep,temp2$Sequence),]
          
          # newdf <- newdf[newdf$Accession != common.proteins[z],]
          # newdf<- rbind(newdf,temp)
          totpepdf <- rbind(totpepdf,temp)
          
          # newdf2 <- newdf2[newdf2$Accession != common.proteins[z],]
          # newdf2<- rbind(newdf2,temp2)
          totpepdf2 <- rbind(totpepdf2,temp2)
        }
        
        
        
        
        totpepdf<- totpepdf[!is.na(rowSums(totpepdf[,1:input$chans])),]
        totpepdf2<- totpepdf2[!is.na(rowSums(totpepdf2[,1:input$chans])),]
        
        totpepdf<- totpepdf[totpepdf$addedVals == 0, ]
        totpepdf<- totpepdf[totpepdf$addedVals == 0, ]
        
        
        totpepdf$uniquePeps[grep(';',totpepdf$Accession)] <- 0
        
        
        ### take sums
        
        pepframe<- data.frame(totpepdf[1:length(common.proteins),1:input$chans],totpepdf2[1:length(common.proteins),1:input$chans],
                              Accession = totpepdf$Accession[1:length(common.proteins)], uniquePeps = totpepdf$uniquePeps[1:length(common.proteins)])
        
        pepsum1 <- pepsum2 <- totpepdf[1:length(common.proteins),]
        pepsum1 <- pepsum1[,-match(c("Sequence","addedVals", "Quality","outliers"),colnames(pepsum1))]
        pepsum2 <- pepsum2[,-match(c("Sequence","addedVals", "Quality","outliers"),colnames(pepsum2))]
        
        pepsum1$pepNum <- pepsum1$pepNum <- rep(0,length(common.proteins))
        
        
        colnames(pepsum2)[1:input$chans]<- channelnames[(input$chans+1):(input$chans*input$reps)]
        
        for(i in 1:length(common.proteins)){
          pepsum1[i,1:input$chans] <- colSums(totpepdf[grep(common.proteins[i],totpepdf$Accession),1:input$chans])
          pepsum2[i,1:input$chans] <- colSums(totpepdf2[grep(common.proteins[i],totpepdf2$Accession),1:input$chans])
          
          pepsum1$pepNum[i] <- nrow(totpepdf[grep(common.proteins[i],totpepdf$Accession),1:input$chans])
          pepsum2$pepNum[i] <- nrow(totpepdf2[grep(common.proteins[i],totpepdf2$Accession),1:input$chans])
          
          
          pepsum1$Accession[i]<- pepsum2$Accession[i] <- common.proteins[i]
          pepsum1$uniquePeps[i]<- totpepdf$uniquePeps[as.logical(match(totpepdf$Accession,common.proteins[i],nomatch = FALSE))][1]
          
          
          
          
          
          #   pepframe[i, 1:input$chans] <- apply( totpepdf[totpepdf$Accession == common.proteins[i], 1:input$chans],2, sum)
          # pepframe[i, (input$chans+1):(input$chans*input$reps)] <- apply( totpepdf2[totpepdf2$Accession == common.proteins[i], 1:input$chans],2, sum)
          #
          # pepframe$Accession <- common.proteins[i]
          # pepframe$uniquePeps <- totpepdf$uniquePeps[totpepdf$Accession == common.proteins[i]][1]
        }
        
        indexpepsum <- ((rowSums(pepsum1[1:input$chans]) != 0) + ((rowSums(pepsum2[1:input$chans])) != 0)) == 2
        pepsum1<- pepsum1[indexpepsum, ]
        pepsum2<- pepsum2[indexpepsum, ]
        
        # indpepsum2 <- pepsum1$pepNum > 1 | pepsum2$pepNum > 1
        # pepsum1 <- pepsum1[indpepsum2,]
        # pepsum2 <- pepsum1[indpepsum2,]
        
        
        fc1 <- pepsum1[,1] / pepsum1[,2:input$chans]
        fc2 <- pepsum2[,1] / pepsum2[,2:input$chans]
        
        
        fcprotdf <- data.frame(log2(fc1),log2(fc2),pepsum1$Accession,uniquePepr1 = pepsum1$uniquePeps,
                               uniquePepr2 = pepsum2$uniquePeps,num1 =  pepsum1$pepNum, num2 = pepsum2$pepNum )
        
        
        fcprotdf
        
      }
      fcprotdf
      
    }else{
      return(NULL)
    }
    
    
  })
  
  ######################################################### isolate removed peptides 
  
  pepdwn <- reactive({
    
    if(input$intcalc == 0){
      return(NULL)
    }else{
      
      
      tempdat <- data()
      
      
      # make matrix of descriptions and accessions to filter by common proteins
      
      accDesMat<- as.character(tempdat[as.character(tempdat[,input$accession]) == "",c(input$accession)])
      
      data.merged <- tempdat[,c(input$view_vars, input$accession, input$sequence,input$quality)]
      
      
      channelnames <- standardNames()
      repindex <- rep(1:input$reps,each = input$chans)
      totfal <- rep(FALSE , (input$chans + 3 ))
      
      
      newdf <- cbind(data.merged[,1:input$chans], data.merged[,c(input$accession, input$sequence,input$quality)])
      colnames(newdf) <- c(channelnames[repindex  == 1],"Accession","Sequence", "Quality")
      newdf <- newdf[!is.na(rowSums(newdf[,1:input$chans])),]
      newdf <- newdf[newdf$Quality <= 0.05, ]
      
      
      
      newdf2 <- cbind(data.merged[,(input$chans+1):(2*input$chans) ],data.merged[,c(input$accession, input$sequence,input$quality)])
      colnames(newdf2) <- c(channelnames[repindex  == 2],"Accession","Sequence","Quality")
      newdf2 <- newdf2[!is.na(rowSums(newdf2[,1:input$chans])),]
      newdf2 <- newdf2[newdf2$Quality <= 0.05, ]
      
      
      
      newdf <- data.frame(newdf, outliers = rep(0,length(newdf[,1])), uniquePeps = rep(0,length(newdf[,1])), addedVals = rep(0,length(newdf[,1])))
      newdf2 <- data.frame(newdf2, outliers = rep(0,length(newdf2[,1])), uniquePeps = rep(0,length(newdf2[,1])), addedVals = rep(0,length(newdf2[,1])))
      
      
      common.proteins <- intersect(unique(newdf$Accession),unique(newdf2$Accession))
      
      newdf <- newdf[!is.na(match(newdf$Accession,common.proteins)),]
      newdf2 <- newdf2[!is.na(match(newdf2$Accession,common.proteins)),]
      
      
      
      totpepdf <- NULL
      totpepdf2 <- NULL
      
      for (z in 1:length(common.proteins)){
        temp <- newdf[newdf$Accession == common.proteins[z],]
        temp2 <- newdf2[newdf2$Accession == common.proteins[z],]
        # first step: check if all peptides are unique ...
        
        if(all.equal(grep(";",temp$Accession),integer(0)) == TRUE){
          
          uniPeptides1 <- length(unique(temp$Sequence))
          
        }  else {
          uniPeptides1 <- length(unique(temp$Sequence[-grep(";",temp$Accession)]))
        }
        
        if(all.equal(grep(";",temp2$Accession),integer(0)) == TRUE){
          
          uniPeptides2 <- length(unique(temp2$Sequence))
          
        }else{
          
          uniPeptides2 <- length(unique(temp2[-grep(";",temp2$Accession)]))
        }
        
        ## add unique peptide column
        
        temp$uniquePeps <- uniPeptides1
        temp2$uniquePeps <- uniPeptides2
        
        
        tempPep <- intersect(unique(temp$Sequence),unique(temp2$Sequence))
        
        
        if(all.equal(tempPep, character(0)) == TRUE){
          next
        }
        
        for(i in 1:length(tempPep)){
          
          if(sum(temp$Sequence == tempPep[i]) != sum(temp2$Sequence == tempPep[i])){
            
            dr <- peptide.match(temp[temp$Sequence == tempPep[i],],temp2[temp2$Sequence == tempPep[i],],input$chans)
            dr1 <- dr$dr1
            dr2 <- dr$dr2
            
          } else {
            
            dr1 <- temp[temp$Sequence == tempPep[i],]
            dr2 <- temp2[temp2$Sequence == tempPep[i],]
            
            
          }
          
          tempoindex <- rep(FALSE,nrow(dr1))
          for(j in 1:nrow(dr1)){
            percor <- cor.test(log2(as.numeric(dr1[j,1:input$chans])),log2(as.numeric(dr2[j,1:input$chans])))
            tempoindex[j] <- percor$estimate < 0.4
          }
          
          dr1[tempoindex,1:input$chans] <- NA
          dr2[tempoindex,1:input$chans] <- NA
          
          
          
          tempReplace <- dr1
          tempReplace2 <- dr2
          
          
          
          temp <- temp[temp$Sequence != tempPep[i],]
          temp <- rbind(temp,tempReplace)
          
          temp2 <- temp2[temp2$Sequence != tempPep[i],]
          temp2 <- rbind(temp2,tempReplace2)
          
        }
        
        temp<- temp[match(tempPep,temp$Sequence),]
        temp2<- temp2[match(tempPep,temp2$Sequence),]
        
        # newdf <- newdf[newdf$Accession != common.proteins[z],]
        # newdf<- rbind(newdf,temp)
        totpepdf <- rbind(totpepdf,temp)
        
        # newdf2 <- newdf2[newdf2$Accession != common.proteins[z],]
        # newdf2<- rbind(newdf2,temp2)
        totpepdf2 <- rbind(totpepdf2,temp2)
      }
      
      
      
      
      totpepdf<- totpepdf[!is.na(rowSums(totpepdf[,1:input$chans])),]
      totpepdf2<- totpepdf2[!is.na(rowSums(totpepdf2[,1:input$chans])),]
      
      totpepdf<- totpepdf[totpepdf$addedVals == 0, ]
      totpepdf<- totpepdf[totpepdf$addedVals == 0, ]
      
      
      totpepdf$uniquePeps[grep(';',totpepdf$Accession)] <- 0
      
      
      ### take sums
      
      pepframe<- data.frame(totpepdf[1:length(common.proteins),1:input$chans],totpepdf2[1:length(common.proteins),1:input$chans],
                            Accession = totpepdf$Accession[1:length(common.proteins)], uniquePeps = totpepdf$uniquePeps[1:length(common.proteins)])
      
      pepsum1 <- pepsum2 <- totpepdf[1:length(common.proteins),]
      pepsum1 <- pepsum1[,-match(c("Sequence","addedVals", "Quality","outliers"),colnames(pepsum1))]
      pepsum2 <- pepsum2[,-match(c("Sequence","addedVals", "Quality","outliers"),colnames(pepsum2))]
      
      pepsum1$pepNum <- pepsum1$pepNum <- rep(0,length(common.proteins))
      
      
      colnames(pepsum2)[1:input$chans]<- channelnames[(input$chans+1):(input$chans*input$reps)]
      
      for(i in 1:length(common.proteins)){
        pepsum1[i,1:input$chans] <- colSums(totpepdf[grep(common.proteins[i],totpepdf$Accession),1:input$chans])
        pepsum2[i,1:input$chans] <- colSums(totpepdf2[grep(common.proteins[i],totpepdf2$Accession),1:input$chans])
        
        pepsum1$pepNum[i] <- nrow(totpepdf[grep(common.proteins[i],totpepdf$Accession),1:input$chans])
        pepsum2$pepNum[i] <- nrow(totpepdf2[grep(common.proteins[i],totpepdf2$Accession),1:input$chans])
        
        
        pepsum1$Accession[i]<- pepsum2$Accession[i] <- common.proteins[i]
        pepsum1$uniquePeps[i]<- totpepdf$uniquePeps[as.logical(match(totpepdf$Accession,common.proteins[i],nomatch = FALSE))][1]
        
      }
      
      indexpepsum <- ((rowSums(pepsum1[1:input$chans]) != 0) + ((rowSums(pepsum2[1:input$chans])) != 0)) == 2
      pepsum1<- pepsum1[indexpepsum, ]
      pepsum2<- pepsum2[indexpepsum, ]
      
      end <- -row.match(newdf[,1:(input$chans+2)],totpepdf[,1:(input$chans+2)])
      end2 <- -row.match(newdf2[,1:(input$chans+2)],totpepdf2[,1:(input$chans+2)])
      
      rempep <- newdf[end,]
      rempep2 <- newdf2[end2,]
      colnames(rempep)<- paste("r1_",colnames(rempep),sep = '')
      colnames(rempep2)<- paste("r2_",colnames(rempep2),sep = '')
      
      
      
      if(nrow(rempep) > nrow(rempep2)){
        
        finalrempep <- rempep
        finalrempep[] <- NA
        finalrempep[1:nrow(rempep2),]<- rempep2
        finalrempep<- cbind(rempep,finalrempep)
      }else{
        finalrempep <- rempep2
        finalrempep[] <- NA
        finalrempep[1:nrow(rempep),] <- rempep
        finalrempep <- cbind(finalrempep,rempep2) 
      }
      colnames(finalrempep)<- c(paste("r1_",colnames(rempep),sep = ''),paste("r2_",colnames(rempep2),sep = ''))
      finalrempep
      
    }
    
  })
  
  ## log fold changes
  
  
  dataMerge <- reactive({
    # testingkd <- rKd()
    # print(head(testingkd))
    inten <- intData()
    
    if(input$datype != "intensity"){
      
      if(input$modtyp == "sigmoid"){
        
        data_orig2 <- data()
        
        ## if user has specified accession & Description 
        if(input$toacc == FALSE){
          pattern<-"GN=(\\S+)"
          g_fromD1<-str_extract(data_orig2$Description,pattern)
          gID_D1a<-str_split_fixed(g_fromD1,"GN=",n=2)
          gID_D1a<-as.vector(gID_D1a[,2])
          gID_D1<-as.matrix(replace(gID_D1a,gID_D1a=="","NA"))
          #Addition of the gene ID column
          data_orig2["geneID"]<- (gID_D1)
        }else{
          if(input$genefile == FALSE){
            tempacc <- data_orig2[,input$accession]
            data_orig2 <- data_orig2[,(-input$accession)]
            data_orig2$Accession <- tempacc
            
            
            if(input$organism == 'H.sapiens'){
              query = '<query model="genomic" view="Protein.primaryAccession Protein.uniprotAccession Protein.genes.symbol" sortOrder="Protein.primaryAccession ASC" >
              <constraint path="Protein.organism.shortName" op="=" value="H. sapiens" code="A" />
              </query>'
              
              ret=POST('http://www.humanmine.org/humanmine/service/query/results',
                       body=list(query=query, format='json'),
                       encode='form')
              
            }else if(input$organism == 'D. melanogaster'){
              query = '<query model="genomic" view="Protein.primaryAccession Protein.uniprotAccession Protein.genes.symbol" sortOrder="Protein.primaryAccession ASC" >
              <constraint path="Protein.organism.shortName" op="=" value="D. melanogaster" code="A" />
              </query>'
              
              ret=POST('http://www.flymine.org/flymine/service/query/results',
                       body=list(query=query, format='json'),
                       encode='form')
            } else if (input$organism == 'R.norvegicus'){
              
              query = '<query model="genomic" view="Protein.primaryAccession Protein.uniprotAccession Protein.genes.symbol" sortOrder="Protein.primaryAccession ASC" >
              <constraint path="Protein.organism.shortName" op="=" value="R. norvegicus" code="A" />
              </query>'
              
              ret = POST('http://www.ratmine.org/ratmine/service/query/results',
                         body=list(query=query, format='json'),
                         encode='form')
            } else if(input$organism == 'C. elegans'){
              query = '<query model="genomic" view="Protein.primaryAccession Protein.uniprotAccession Protein.genes.symbol" sortOrder="Protein.primaryAccession ASC" >
          <constraint path="Protein.organism.shortName" op="=" value="C. elegans" code="A" />
          </query>'
              
              et=POST('http://www.humanmine.org/humanmine/service/query/results',
                      body=list(query=query, format='json'),
                      encode='form')
            }
            
            
            response <- jsonlite::fromJSON(httr::content(ret,as='text'))
            
            data.prots <- response$results
            
            GeneID <- data.prots[match(data_orig2$Accession,data.prots[,1]),3]
            GeneID <- make.names(GeneID,unique = TRUE)
            
            data_orig2$geneID <- GeneID
            
            # uniGene <- uniprotGene()
            # uniGene$Gene.names <- gsub(' .*','',uniGene$Gene.names)
            # uniGene$Gene.names[uniGene$Gene.names == '']<- NA
            # GeneID <- uniGene$Gene[match(data_orig2$Accession,uniGene$Entry)]
            # GeneID <- make.names(GeneID,unique = TRUE)
            # 
            # data_orig2$geneID <- GeneID
          }else{
            
            tempacc <- data_orig2[,input$accession]
            data_orig2 <- data_orig2[,(-input$accession)]
            data_orig2$Accession <- tempacc
            
            uniGene <- uploadGene()
            GeneID <- uniGene$Gene[match(data_orig2$Accession,uniGene$Entry)]
            GeneID <- make.names(GeneID,unique = TRUE)
            
            data_orig2$geneID <- GeneID
            
          }
        }
        
        if(input$incpd == TRUE){
          
          if(input$toacc == FALSE){
            data.merged<-data.frame(data_orig2[,input$view_vars],
                                    
                                    data_orig2$Accession, data_orig2$geneID,
                                    data_orig2$X..Unique.Peptides,data_orig2[,input$pdofpd])
          }else{
            data.merged<-data.frame(data_orig2[,input$view_vars],
                                    
                                    data_orig2$Accession, data_orig2$geneID,
                                    data_orig2[,input$unipeps],data_orig2[,input$pdofpd])
          }
          final.Names <- finalNames()
          colnames(data.merged)<- c(final.Names,
                                    "Accession", "GeneID","UniquePeps",'Kd')
          
          
          tmp<-data.merged[,1:input$chans]
          tf <- as.data.frame(lapply(tmp, function(x) (is.na(x))))
          countsNAs<-as.data.frame(apply(tf,1,function(x)table(x)["TRUE"]))
          n_of_miss<-as.data.frame(as.numeric(str_replace_all(as.list(countsNAs[,1]),"NA",0)))
          data.merged <- data.frame(data.merged,n_of_miss)
          colnames(data.merged)<-c( final.Names,
                                    "Accession", "GeneID","UniquePeps",'Kd', "MissingVal")
          #
          missing_val <- 0
          #
          data.merged<- as.data.frame(data.merged[data.merged$MissingVal <= missing_val, ])
          #
          # #Specify the number of missing points.For zero missing point  is ==0
          #
          # #filiter for 2 unique peptides
          data.merged <- data.merged[data.merged$UniquePeps > 1, ]
          #
          #
          data.merged <- data.frame((normalize.loess(2^(data.merged[,channels()]))),
                                    Accession=data.merged$Accession,
                                    GeneID=data.merged$GeneID,
                                    UniquePeps=data.merged$UniquePeps,
                                    Kd = 1 / data.merged$Kd,
                                    MissingVal=data.merged$MissingVal)
          
                data.merged
        
        } else {
          if(input$toacc == FALSE){
            data.merged<-data.frame(data_orig2[,input$view_vars],
                                    
                                    data_orig2$Accession, data_orig2$geneID,
                                    data_orig2$X..Unique.Peptides)
            
          }else{
            data.merged<-data.frame(data_orig2[,input$view_vars],
                                    
                                    data_orig2$Accession, data_orig2$geneID,
                                    data_orig2[,input$unipeps])
          }
          final.Names <- finalNames()
          
          colnames(data.merged)<- c(final.Names,
                                    "Accession", "GeneID","UniquePeps")
          
          
          tmp<-data.merged[,1:input$chans]
          tf <- as.data.frame(lapply(tmp, function(x) (is.na(x))))
          countsNAs<-as.data.frame(apply(tf,1,function(x)table(x)["TRUE"]))
          n_of_miss<-as.data.frame(as.numeric(str_replace_all(as.list(countsNAs[,1]),"NA",0)))
          data.merged <- data.frame(data.merged,n_of_miss)
          colnames(data.merged)<-c( final.Names,
                                    "Accession", "GeneID","UniquePeps","MissingVal")
          #
          missing_val <- 0
          #
          data.merged<- as.data.frame(data.merged[data.merged$MissingVal <= missing_val, ])
          #
          # #Specify the number of missing points.For zero missing point  is ==0
          #
          # #filiter for 2 unique peptides
          data.merged <- data.merged[data.merged$UniquePeps > 1, ]
          #
          #
          data.merged <- data.frame((normalize.loess(2^(data.merged[,channels()]))),
                                    Accession=data.merged$Accession,
                                    GeneID=data.merged$GeneID,
                                    UniquePeps=data.merged$UniquePeps,
                                    MissingVal=data.merged$MissingVal)
          data.merged
        }
        
      } else {
        tempdat <- data()
        temp<- channels()
        nvec <- length(temp)
        
        data.merged <- tempdat[,input$view_vars]
        colnames(data.merged) <- finalNames()
        
        data.names<- c(input$view_vars,"Accession", "GeneID","UniquePeps","MissingVal")
        
        if(input$datype == "FC"){
          data.merged[,temp] <- log2(data.merged[,temp])
        }
        #tidy below !!
        
        if(input$toacc == FALSE){
          
        
          pattern<-"GN=(\\S+)"
          g_fromD1<-str_extract(tempdat$Description,pattern)
          gID_D1a<-str_split_fixed(g_fromD1,"GN=",n=2)
          gID_D1a<-as.vector(gID_D1a[,2])
          gID_D1<-as.matrix(replace(gID_D1a,gID_D1a=="","NA"))
          
          Accession <- tempdat$Accession
          UniquePeps <- MissingVal<- tempdat[,grep("Unique",colnames(tempdat))]
          data.merged <- cbind(data.merged,Accession, GeneID = gID_D1,UniquePeps,MissingVal)
          
        } else {
          
          tempdat$Accession <- tempdat[,input$accession]
          if(input$genefile == FALSE){
            
            # uniGene <- uniprotGene()
            # uniGene$Gene.names <- gsub(' .*','',uniGene$Gene.names)
            # print(head(uniGene))
            tempdat$Accession <- tempdat[,input$accession]
            if(input$genefile == FALSE){
              
              # uniGene <- uniprotGene()
              # uniGene$Gene.names <- gsub(' .*','',uniGene$Gene.names)
              
              if(input$organism == 'H.sapiens'){
                query = '<query model="genomic" view="Protein.primaryAccession Protein.uniprotAccession Protein.genes.symbol" sortOrder="Protein.primaryAccession ASC" >
                <constraint path="Protein.organism.shortName" op="=" value="H. sapiens" code="A" />
                </query>'
                
                ret=POST('http://www.humanmine.org/humanmine/service/query/results',
                         body=list(query=query, format='json'),
                         encode='form')
                
              }else if(input$organism == 'D.melanogaster'){
                query = '<query model="genomic" view="Protein.primaryAccession Protein.uniprotAccession Protein.genes.symbol" sortOrder="Protein.primaryAccession ASC" >
                <constraint path="Protein.organism.shortName" op="=" value="D. melanogaster" code="A" />
                </query>'
                
                ret=POST('http://www.flymine.org/flymine/service/query/results',
                         body=list(query=query, format='json'),
                         encode='form')
              }else if (input$organism == 'R.norvegicus'){
                
                query = '<query model="genomic" view="Protein.primaryAccession Protein.uniprotAccession Protein.genes.symbol" sortOrder="Protein.primaryAccession ASC" >
                <constraint path="Protein.organism.shortName" op="=" value="R. norvegicus" code="A" />
                </query>'
                
                ret = POST('http://www.ratmine.org/ratmine/service/query/results',
                           body=list(query=query, format='json'),
                           encode='form')
              }else if(input$organism == 'C. elegans'){
                query = '<query model="genomic" view="Protein.primaryAccession Protein.uniprotAccession Protein.genes.symbol" sortOrder="Protein.primaryAccession ASC" >
          <constraint path="Protein.organism.shortName" op="=" value="C. elegans" code="A" />
                </query>'
                
                et=POST('http://www.humanmine.org/humanmine/service/query/results',
                        body=list(query=query, format='json'),
                        encode='form')
              }
              
              
              response = jsonlite::fromJSON(httr::content(ret,as='text'))
              
              data.prots <- response$results
              
              GeneID <- data.prots[match(tempdat$Accession,data.prots[,1]),3]
              GeneID <- make.names(GeneID,unique = TRUE)
          } else{
            uniGene <- uploadGene()
            colnames(uniGene)[1:2] <- c('Entry','Gene.names')
            uniGene$Gene.names[uniGene$Gene.names == '']<- NA
            GeneID <- uniGene$Gene[match(tempdat$Accession,uniGene$Entry)]
            GeneID <- make.names(GeneID,unique = TRUE)
          }
          
          Accession <- tempdat$Accession
          UniquePeps <- MissingVal<- tempdat[,input$unipeps]
          
          # data_orig2$geneID <- GeneID
          data.merged <- cbind(data.merged,Accession, GeneID = GeneID,UniquePeps,MissingVal)
          
        }
        
        
        
        missing<- rowSums(is.na(data.merged[,channels()]))
        
        missing_val <- 0
        data.merged$MissingVal <- missing
        
        ## subset by missing
        data.merged <- data.merged[data.merged$MissingVal <= missing_val, ]
        data.merged<- as.data.frame(data.merged[data.merged$MissingVal <= missing_val, ])
        
        data.merged <- data.merged[data.merged$UniquePeps > 1, ]
        #
        #
        data.merged <- data.frame(log2(normalize.loess(2^(data.merged[,channels()]))),
                                  Accession=data.merged$Accession,
                                  GeneID=data.merged$GeneID,
                                  UniquePeps=data.merged$UniquePeps,
                                  MissingVal=data.merged$MissingVal)
        data.merged
        ## filter for 2 unique peptide
      }
      
    } } else {
      
      ### intensities to protein done.. here no description so we need toi use accession and intermine
      
      tempdat <- inten
      
      if(input$genefile == FALSE){
        
        # uniGene <- uniprotGene()
        # uniGene$Gene.names <- gsub(' .*','',uniGene$Gene.names)
        if(input$organism == 'H.sapiens'){
          query = '<query model="genomic" view="Protein.primaryAccession Protein.uniprotAccession Protein.genes.symbol" sortOrder="Protein.primaryAccession ASC" >
          <constraint path="Protein.organism.shortName" op="=" value="H. sapiens" code="A" />
          </query>'
          
          ret=POST('http://www.humanmine.org/humanmine/service/query/results',
                   body=list(query=query, format='json'),
                   encode='form')
          
        }else if(input$organism == 'C. elegans'){
          query = '<query model="genomic" view="Protein.primaryAccession Protein.uniprotAccession Protein.genes.symbol" sortOrder="Protein.primaryAccession ASC" >
          <constraint path="Protein.organism.shortName" op="=" value="C. elegans" code="A" />
          </query>'
          
          et=POST('http://www.humanmine.org/humanmine/service/query/results',
                  body=list(query=query, format='json'),
                  encode='form')
        }else if(input$organism == 'D. melanogaster'){
          query = '<query model="genomic" view="Protein.primaryAccession Protein.uniprotAccession Protein.genes.symbol" sortOrder="Protein.primaryAccession ASC" >
          <constraint path="Protein.organism.shortName" op="=" value="D. melanogaster" code="A" />
          </query>'
          
          ret = POST('http://www.flymine.org/flymine/service/query/results',
                     body=list(query=query, format='json'),
                     encode='form')
        } else if (input$organism == 'R.norvegicus'){
          
          query = '<query model="genomic" view="Protein.primaryAccession Protein.uniprotAccession Protein.genes.symbol" sortOrder="Protein.primaryAccession ASC" >
          <constraint path="Protein.organism.shortName" op="=" value="R. norvegicus" code="A" />
          </query>'
          
          ret = POST('http://www.ratmine.org/ratmine/service/query/results',
                     body=list(query=query, format='json'),
                     encode='form')
        }
        
        
        response = jsonlite::fromJSON(httr::content(ret,as='text'))
        
        data.prots <- response$results
        
        GeneID <- data.prots[match(tempdat$pepsum1.Accession,data.prots[,1]),3]
        GeneID <- make.names(GeneID,unique = TRUE)
      }else{
        uniGene <- uploadGene()
        colnames(uniGene) <- c('Entry','Gene.names')
      }
      
      if(input$dorem =='yes'){
        
        Accession <- tempdat$pepsum1.Accession
      }else{
        Accession <- tempdat$Accession
      }
      uniGene$Gene.names[uniGene$Gene.names == '']<- NA
      GeneID <- uniGene$Gene[match(Accession,uniGene$Entry)]
      GeneID <- make.names(GeneID,unique = TRUE)
      
      tempdat$GeneID <- GeneID
      
      UniquePeps <- tempdat$uniquePepr1
      
      if(input$reps == 1){
        
        data.merged <- tempdat
        final.Names <- paste0('rep1_C',0:(input$chans - 2))
        if(input$incpd == TRUE){
          
          data.merged <- data.frame(data.merged[,1:(input$chans -1 )], Accession = data.merged$Accession,
                                    GeneID = data.merged$GeneID, UniquePeps = data.merged$uniquePepr1, Kd = data.merged$Kd)
          colnames(data.merged)<- c(final.Names,
                                    "Accession", "GeneID","UniquePeps", 'Kd')
          
          tmp<-data.merged[,1:(input$chans -1 )]
          tf <- as.data.frame(lapply(tmp, function(x) (is.na(x))))
          countsNAs<-as.data.frame(apply(tf,1,function(x)table(x)["TRUE"]))
          n_of_miss<-as.data.frame(as.numeric(str_replace_all(as.list(countsNAs[,1]),"NA",0)))
          data.merged <- data.frame(data.merged,n_of_miss)
          colnames(data.merged)<-c( final.Names,
                                    "Accession", "GeneID","UniquePeps",'Kd',"MissingVal")
          #
          missing_val <- 0
          #
          data.merged<- as.data.frame(data.merged[data.merged$MissingVal <= missing_val, ])
          #
          # #Specify the number of missing points.For zero missing point  is ==0
          #
          # #filiter for 2 unique peptides
          data.merged <- data.merged[data.merged$UniquePeps > 1, ]
          #
          #
          data.merged <- data.frame((normalize.loess(2^(data.merged[,1:(input$chans -1)]))),
                                    Accession=data.merged$Accession,
                                    GeneID=data.merged$GeneID,
                                    UniquePeps=data.merged$UniquePeps,
                                    Kd = data.merged$Kd
          )
          
          
          data.merged
          
        } else {
          
          data.merged <- data.frame(data.merged[,1:(input$chans -1 )], Accession = data.merged$Accession,
                                    GeneID = data.merged$GeneID, UniquePeps = data.merged$uniquePepr1)
          colnames(data.merged)<- c(final.Names,
                                    "Accession", "GeneID","UniquePeps")
          
          
          tmp<-data.merged[,1:(input$chans -1 )]
          tf <- as.data.frame(lapply(tmp, function(x) (is.na(x))))
          countsNAs<-as.data.frame(apply(tf,1,function(x)table(x)["TRUE"]))
          n_of_miss<-as.data.frame(as.numeric(str_replace_all(as.list(countsNAs[,1]),"NA",0)))
          data.merged <- data.frame(data.merged,n_of_miss)
          colnames(data.merged)<-c( final.Names,
                                    "Accession", "GeneID","UniquePeps","MissingVal")
          #
          missing_val <- 0
          #
          data.merged<- as.data.frame(data.merged[data.merged$MissingVal <= missing_val, ])
          #
          # #Specify the number of missing points.For zero missing point  is ==0
          #
          # #filiter for 2 unique peptides
          data.merged <- data.merged[data.merged$UniquePeps > 1, ]
          #
          #
          data.merged <- data.frame((normalize.loess(2^(data.merged[,1:(input$chans -1)]))),
                                    Accession=data.merged$Accession,
                                    GeneID=data.merged$GeneID,
                                    UniquePeps=data.merged$UniquePeps)
          
          
          data.merged
        }
        
      } else{
        
        data.merged <- tempdat[,1:(input$chans*input$reps - 2)]
        data.merged <- tempdat[,1:(input$chans*input$reps - 2)]
        if( input$dorem == 'no'){
          Accession <- tempdat$Accession
          MissingVal <- UniquePeps
          data.merged <- cbind(data.merged,Accession, GeneID = GeneID,UniquePeps,MissingVal)
          
          data.merged <- data.merged[!is.na(rowSums(data.merged[,1:(input$chans*input$reps -2)])),]
          data.merged <- data.merged[data.merged$UniquePeps > 1, ]
          
          data.merged <- data.frame((log2(normalize.loess(2^(data.merged[,1:(input$chans*input$reps - 2)])))),
                                    Accession=data.merged$Accession,
                                    GeneID=data.merged$GeneID,
                                    UniquePeps=data.merged$UniquePeps)
          
          data.merged <- data.merged[!is.na(rowSums(data.merged[,1:(input$chans*input$reps -2)])),]
          data.merged
        }else{
          
          Accession <- tempdat$pepsum1.Accession
          data.merged <- cbind(data.merged,Accession, GeneID = GeneID,UniquePeps, num1 = tempdat$num1, num2 = tempdat$num2)
          
          
          data.merged <- data.merged[!is.na(rowSums(data.merged[,1:(input$chans*input$reps -2)])),]
          data.merged <- data.merged[data.merged$UniquePeps > 1, ]
          
          data.merged <- data.frame((log2(normalize.loess(2^(data.merged[,1:(input$chans*input$reps - 2)])))),
                                    Accession=data.merged$Accession,
                                    GeneID=data.merged$GeneID,
                                    UniquePeps=data.merged$UniquePeps, num1 = data.merged$num1, num2 = data.merged$num2)
          
          data.merged <- data.merged[data.merged$num1 > 1 | data.merged$num2 > 1 , ]
          data.merged <- data.merged[!is.na(rowSums(data.merged[,1:(input$chans*input$reps -2)])),]
          data.merged
        }
        
              }
      
     
    }
    

    data.merged
  })
  
  sigFinalNames <- reactive({
    paste0("C_",1:input$chans)
  })
  
  sigPredNames <- reactive({
    paste0("predX",1:input$chans)
  })
  
  dataMerge2 <- reactive({
    nvec <- channels()
    nvec <- length(nvec)
    data.merged <- dataMerge()
    
    if(input$modtyp == "sigmoid"){
      if(input$datype == 'intensity'){
        
        
        conc <- sigConc()
        if(input$incpd == TRUE){

          final.Names <- paste0('rep1_C',0:(input$chans - 2 ))
          pred.names <- paste0('predX',1:(input$chans -1))
          colnames(data.merged[,1:(input$chans -1)]) <- final.Names
          data_merged_positives<- data.merged
          # na.omit(data.merged[data.merged[,1:input$chans] >= 1,])
          
          data_merged_positives2<-( (1/data_merged_positives[,1:(input$chans -1 )]))*100
          Reps_FC<-data.frame(data_merged_positives2 ,
                              Accession =  data_merged_positives$Accession,GeneID = data_merged_positives$GeneID,
                              UniquePeps = data_merged_positives$UniquePeps, depletionConstant = data_merged_positives$Kd
          )
          
          
          ryegrass.m1<- vector(mode = "list",length = nrow(Reps_FC))
          pvals<-list()
          stderr<-list()
          model_pred<-list()
          coeff_predicted<-list()
          for(i in 1:nrow(Reps_FC)){
            #print(i)
            
            #nrow(full_df_2)
            #maxIt and relTol to be user defined
            ryegrass.m1[[i]]<-try(drm(as.numeric(Reps_FC[i,1:(input$chans - 1 )]) ~ as.numeric(conc),
                                      na.action = na.omit,
                                      control = drmc(constr = FALSE, errorm = FALSE, noMessage = TRUE, maxIt = 1000, relTol = 1e-06),
                                      fct = LL.4(fixed=c(NA, NA, NA, NA), #see note @top this file
                                                 names = c("Slope", "Lower Limit", "Upper Limit", "RB50"))),silent = T)
            
          }
          
          failed_sigm=0
          for(i in 1:length(ryegrass.m1)){
            #print(i)
            #checking_val if FALSE  the model has failed to calculate the pval
            checking_val<-try(is.numeric(coefficients(ryegrass.m1[[i]])[["Slope:(Intercept)"]]),silent = T)
            
            if(checking_val=="TRUE"){
              #print(checking_val)
              pvals[[i]]<-t(as.data.frame(summary(ryegrass.m1[[i]])$coefficients[13:16]))
              colnames(pvals[[i]])<-c("SlopePval", "Lower_LimitPval", "Upper_LimitPval", "RB50Pval")
              coeff_predicted[[i]]<-t(data.frame(coefficients(ryegrass.m1[[i]])))
              colnames(coeff_predicted[[i]])<-c("SlopeCoef", "Lower_LimitCoef", "Upper_LimitCoef","RB50Coef")
              stderr[[i]]<-t(as.data.frame(summary(ryegrass.m1[[i]])$coefficients[5:8]))
              colnames(stderr[[i]])<-c("SlopeErr", "Lower_LimitErr", "Upper_LimitErr","RB50Err")
              
              model_pred[[i]]<-predict(ryegrass.m1[[i]])
            }else{
              failed_sigm=failed_sigm+1
              
              fit <- lm(as.numeric(Reps_FC[i,1:(input$chans -1)]) ~ poly(log10(conc),2 ))
              #extract the pval
              pval<-as.numeric(summary(fit)$coefficients[,4] )
              pvals[[i]]<- t(as.data.frame(c(pval,"lm-fit:intercept.slope.quadratic") ))
              colnames(pvals[[i]])<-c("SlopePval", "Lower_LimitPval", "Upper_LimitPval","RB50Pval")
              stderr[[i]]<- data.frame(NA,NA,NA,NA)
              colnames(stderr[[i]])<-c("SlopeErr", "Lower_LimitErr", "Upper_LimitErr","RB50Err")
              coeff_predicted[[i]]<-data.frame(NA,NA,NA,NA)
              colnames(coeff_predicted[[i]])<-c("SlopeCoef", "Lower_LimitCoef", "Upper_LimitCoef","RB50Coef")
              model_pred[[i]]<- as.numeric(fitted(fit))
            } #just adding NAs for the times the model failed
          }
          
          
          modelsReps<-data.frame(
            do.call(rbind.data.frame,lapply(model_pred,function(x) as.numeric(x))),
            Reps_FC$GeneID,
            do.call(rbind.data.frame,lapply(pvals,function(x) x) ),
            do.call(rbind.data.frame,lapply(coeff_predicted,function(x) x)) ,
            do.call(rbind.data.frame,lapply(stderr,function(x) x) )
          )
          colnames(modelsReps)<-c(pred.names,"GeneID",
                                  "SlopePval", "Lower_LimitPval","Upper_LimitPval", "RB50Pval",
                                  "SlopeCoef", "Lower_LimitCoef","Upper_LimitCoef", "RB50Coef",
                                  "SlopeErr", "Lower_LimitErr","Upper_LimitErr", "RB50Err"
          )
          
          data_merged_2 <-merge.data.frame(modelsReps,Reps_FC,by = 'GeneID')
          
          
          data_merged_2<-data.frame(data_merged_2,"Top_minus_min"=data_merged_2$predX1-data_merged_2[,paste("predX",(input$chans-1),sep = "")])
          
          crap <- crapome()
          # 
          tempcrap <- crap$Gene
          tempcrap<- crap$percentagePresent[match(data_merged_2$GeneID,crap$Gene)]
          data_merged_2 <- data.frame(data_merged_2, CRAPomePercent = tempcrap)
          
          proteome<-as.vector(toupper(data_merged_2$GeneID))
          mykinases <- intersect(proteome,kinome)
          
          if(!is.null(input$venninp)){
           myvenninp <- uploadVenn()
           uploadsInt <- intersect(proteome,myvenninp)
           
           myvenvec <- match(proteome, uploadsInt, nomatch = NA)
           
           ## now combine with my kin vec 
           
           mykinvec <- match(proteome,mykinases, nomatch = NA)
           mykinvec[!is.na(mykinvec) & !is.na(myvenvec)] <- 'KINASE + UPLOAD' 
           mykinvec[!is.na(mykinvec) & is.na(myvenvec)] <- 'KINASE'
           mykinvec[is.na(mykinvec) & !is.na(myvenvec)] <- 'UPLOAD'
           print(mykinvec)
           
          } else {

            mykinvec <- match(proteome,mykinases, nomatch = NA)
            mykinvec[!is.na(mykinvec)]<- 'KINASE'
          
          }
          
          data_merged_2<- data.frame(data_merged_2, correctedRB50 = (data_merged_2$RB50Coef*data_merged_2$depletionConstant),Kinase = mykinvec)
          data_merged_2
          
        } else {
          
          final.Names <- paste0('rep1_C',0:(input$chans - 2 ))
          pred.names <- paste0('predX',1:(input$chans -1))
          colnames(data.merged[,1:(input$chans -1)]) <- final.Names
          data_merged_positives<- data.merged
          # na.omit(data.merged[data.merged[,1:input$chans] >= 1,])
          
          data_merged_positives2<-( (1/data_merged_positives[,1:(input$chans -1 )]))*100
          Reps_FC<-data.frame(data_merged_positives2 ,
                              Accession =  data_merged_positives$Accession,GeneID = data_merged_positives$GeneID,
                              UniquePeps = data_merged_positives$UniquePeps
          )
          
          
          ryegrass.m1<- vector(mode = "list",length = nrow(Reps_FC))
          pvals<-list()
          stderr<-list()
          model_pred<-list()
          coeff_predicted<-list()
          for(i in 1:nrow(Reps_FC)){
            #print(i)
            
            #nrow(full_df_2)
            #maxIt and relTol to be user defined
            ryegrass.m1[[i]]<-try(drm(as.numeric(Reps_FC[i,1:(input$chans - 1 )]) ~ as.numeric(conc),
                                      na.action = na.omit,
                                      control = drmc(constr = FALSE, errorm = FALSE, noMessage = TRUE, maxIt = 1000, relTol = 1e-06),
                                      fct = LL.4(fixed=c(NA, NA, NA, NA), #see note @top this file
                                                 names = c("Slope", "Lower Limit", "Upper Limit", "RB50"))),silent = T)
            
          }
          
          failed_sigm=0
          for(i in 1:length(ryegrass.m1)){
            #print(i)
            #checking_val if FALSE  the model has failed to calculate the pval
            checking_val<-try(is.numeric(coefficients(ryegrass.m1[[i]])[["Slope:(Intercept)"]]),silent = T)
            
            if(checking_val=="TRUE"){
              #print(checking_val)
              pvals[[i]]<-t(as.data.frame(summary(ryegrass.m1[[i]])$coefficients[13:16]))
              colnames(pvals[[i]])<-c("SlopePval", "Lower_LimitPval", "Upper_LimitPval", "RB50Pval")
              coeff_predicted[[i]]<-t(data.frame(coefficients(ryegrass.m1[[i]])))
              colnames(coeff_predicted[[i]])<-c("SlopeCoef", "Lower_LimitCoef", "Upper_LimitCoef","RB50Coef")
              stderr[[i]]<-t(as.data.frame(summary(ryegrass.m1[[i]])$coefficients[5:8]))
              colnames(stderr[[i]])<-c("SlopeErr", "Lower_LimitErr", "Upper_LimitErr","RB50Err")
              
              model_pred[[i]]<-predict(ryegrass.m1[[i]])
            }else{
              failed_sigm=failed_sigm+1
              
              fit <- lm(as.numeric(Reps_FC[i,1:(input$chans -1)]) ~ poly(log10(conc),2 ))
              #extract the pval
              pval<-as.numeric(summary(fit)$coefficients[,4] )
              pvals[[i]]<- t(as.data.frame(c(pval,"lm-fit:intercept.slope.quadratic") ))
              colnames(pvals[[i]])<-c("SlopePval", "Lower_LimitPval", "Upper_LimitPval","RB50Pval")
              stderr[[i]]<- data.frame(NA,NA,NA,NA)
              colnames(stderr[[i]])<-c("SlopeErr", "Lower_LimitErr", "Upper_LimitErr","RB50Err")
              coeff_predicted[[i]]<-data.frame(NA,NA,NA,NA)
              colnames(coeff_predicted[[i]])<-c("SlopeCoef", "Lower_LimitCoef", "Upper_LimitCoef","RB50Coef")
              model_pred[[i]]<- as.numeric(fitted(fit))
            } #just adding NAs for the times the model failed
          }
          
          
          modelsReps<-data.frame(
            do.call(rbind.data.frame,lapply(model_pred,function(x) as.numeric(x))),
            Reps_FC$GeneID,
            do.call(rbind.data.frame,lapply(pvals,function(x) x) ),
            do.call(rbind.data.frame,lapply(coeff_predicted,function(x) x)) ,
            do.call(rbind.data.frame,lapply(stderr,function(x) x) )
          )
          colnames(modelsReps)<-c(pred.names,"GeneID",
                                  "SlopePval", "Lower_LimitPval","Upper_LimitPval", "RB50Pval",
                                  "SlopeCoef", "Lower_LimitCoef","Upper_LimitCoef", "RB50Coef",
                                  "SlopeErr", "Lower_LimitErr","Upper_LimitErr", "RB50Err"
          )
          
          data_merged_2 <-merge.data.frame(modelsReps,Reps_FC,by = 'GeneID')
          
          
          data_merged_2<-data.frame(data_merged_2,"Top_minus_min"=data_merged_2$predX1-data_merged_2[,paste("predX",(input$chans-1),sep = "")])
          
          ## crapome integration 
          
          crap <- crapome()
          # 
          tempcrap <- crap$Gene
          tempcrap<- crap$percentagePresent[match(data_merged_2$GeneID,crap$Gene)]
          data_merged_2 <- data.frame(data_merged_2, CRAPomePercent = tempcrap)
          
          proteome<-as.vector(toupper(data_merged_2$GeneID))
          mykinases <- intersect(proteome,kinome)
          
          mykinvec <- match(proteome,mykinases, nomatch = NA)
          mykinvec[!is.na(mykinvec)]<- 'KINASE'
          
          
          data_merged_2<- data.frame(data_merged_2,Kinase = mykinvec)
          data_merged_2
   
        }
        
        
      } else{

        conc <- sigConc()
        final.Names <- sigFinalNames()
        pred.names <- sigPredNames()
        colnames(data.merged[,1:input$chans]) <- final.Names
        data_merged_positives<- data.merged
        # na.omit(data.merged[data.merged[,1:input$chans] >= 1,])
        
        data_merged_positives2<-( (1/data_merged_positives[,1:input$chans]))*100
       if(input$incpd == TRUE){
         Reps_FC<-data.frame(data_merged_positives2 ,
                             Accession =  data_merged_positives$Accession,GeneID = data_merged_positives$GeneID,
                             UniquePeps = data_merged_positives$UniquePeps, depletionConstant = data_merged_positives$Kd, MissingVal = data_merged_positives$MissingVal
         )
       }else{
         Reps_FC<-data.frame(data_merged_positives2 ,
                             Accession =  data_merged_positives$Accession,GeneID = data_merged_positives$GeneID,
                             UniquePeps = data_merged_positives$UniquePeps, MissingVal = data_merged_positives$MissingVal
         )
       }
        
        ryegrass.m1<- vector(mode = "list",length = nrow(Reps_FC))
        pvals<-list()
        stderr<-list()
        model_pred<-list()
        coeff_predicted<-list()
        for(i in 1:nrow(Reps_FC)){
          #print(i)
          
          #nrow(full_df_2)
          #maxIt and relTol to be user defined
          ryegrass.m1[[i]]<-try(drm(as.numeric(Reps_FC[i,1:input$chans]) ~ as.numeric(conc),
                                    na.action = na.omit,
                                    control = drmc(constr = FALSE, errorm = FALSE, noMessage = TRUE, maxIt = 1000, relTol = 1e-06),
                                    fct = LL.4(fixed=c(NA, NA, NA, NA), #see note @top this file
                                               names = c("Slope", "Lower Limit", "Upper Limit", "RB50"))),silent = T)
          
        }
        
        failed_sigm=0
        for(i in 1:length(ryegrass.m1)){
          #print(i)
          #checking_val if FALSE  the model has failed to calculate the pval
          checking_val<-try(is.numeric(coefficients(ryegrass.m1[[i]])[["Slope:(Intercept)"]]),silent = T)
          
          if(checking_val=="TRUE"){
            #print(checking_val)
            pvals[[i]]<-t(as.data.frame(summary(ryegrass.m1[[i]])$coefficients[13:16]))
            colnames(pvals[[i]])<-c("SlopePval", "Lower_LimitPval", "Upper_LimitPval", "RB50Pval")
            coeff_predicted[[i]]<-t(data.frame(coefficients(ryegrass.m1[[i]])))
            colnames(coeff_predicted[[i]])<-c("SlopeCoef", "Lower_LimitCoef", "Upper_LimitCoef","RB50Coef")
            stderr[[i]]<-t(as.data.frame(summary(ryegrass.m1[[i]])$coefficients[5:8]))
            colnames(stderr[[i]])<-c("SlopeErr", "Lower_LimitErr", "Upper_LimitErr","RB50Err")
            
            model_pred[[i]]<-predict(ryegrass.m1[[i]])
          }else{
            failed_sigm=failed_sigm+1
            
            fit <- lm(as.numeric(Reps_FC[i,1:input$chans]) ~ poly(log10(conc),2 ))
            #extract the pval
            pval<-as.numeric(summary(fit)$coefficients[,4] )
            pvals[[i]]<- t(as.data.frame(c(pval,"lm-fit:intercept.slope.quadratic") ))
            colnames(pvals[[i]])<-c("SlopePval", "Lower_LimitPval", "Upper_LimitPval","RB50Pval")
            stderr[[i]]<- data.frame(NA,NA,NA,NA)
            colnames(stderr[[i]])<-c("SlopeErr", "Lower_LimitErr", "Upper_LimitErr","RB50Err")
            coeff_predicted[[i]]<-data.frame(NA,NA,NA,NA)
            colnames(coeff_predicted[[i]])<-c("SlopeCoef", "Lower_LimitCoef", "Upper_LimitCoef","RB50Coef")
            model_pred[[i]]<- as.numeric(fitted(fit))
          } #just adding NAs for the times the model failed
        }
        
        
        modelsReps<-data.frame(
          do.call(rbind.data.frame,lapply(model_pred,function(x) as.numeric(x))),
          Reps_FC$GeneID,
          do.call(rbind.data.frame,lapply(pvals,function(x) x) ),
          do.call(rbind.data.frame,lapply(coeff_predicted,function(x) x)) ,
          do.call(rbind.data.frame,lapply(stderr,function(x) x) )
        )
        colnames(modelsReps)<-c(pred.names,"GeneID",
                                "SlopePval", "Lower_LimitPval","Upper_LimitPval", "RB50Pval",
                                "SlopeCoef", "Lower_LimitCoef","Upper_LimitCoef", "RB50Coef",
                                "SlopeErr", "Lower_LimitErr","Upper_LimitErr", "RB50Err"
        )
        data_merged_2<-merge.data.frame(modelsReps,Reps_FC,by="GeneID")
        
          data_merged_2<-data.frame(data_merged_2,"Top_minus_min"=data_merged_2$predX1-data_merged_2[,paste("predX",input$chans,sep = "")])
    
        
        ## crapome integration 
        
    
        crap <- crapome()
        # 
        tempcrap <- crap$Gene
        tempcrap<- crap$percentagePresent[match(data_merged_2$GeneID,crap$Gene)]
        data_merged_2 <- data.frame(data_merged_2, CRAPomePercent = tempcrap)
        
        proteome<-as.vector(toupper(data_merged_2$GeneID))
        mykinases <- intersect(proteome,kinome)
        
        mykinvec <- match(proteome,mykinases, nomatch = NA)
        mykinvec[!is.na(mykinvec)]<- 'KINASE'
        
        if(input$incpd == TRUE){
          
          data_merged_2<- data.frame(data_merged_2, correctedRB50 = (data_merged_2$RB50Coef*data_merged_2$depletionConstant),Kinase = mykinvec)
        }else{
          data_merged_2<- data.frame(data_merged_2,Kinase = mykinvec)
        }  
        data_merged_2

      }
      
    } else{
      if(input$datype == "intensity"){
        conc<- rep(0:(input$chans - 2), times = input$reps)
        
      }else{
        
        conc<- rep(0:(input$chans - 1), times = rReps())[channels()]
        
      }
      
      design<-model.matrix(~poly(conc,2))
      colnames(design)<-c("Intercept","Slope","Quadratic")
      
      # reactive start
      
      fit <- lmFit(data.merged[,1:length(conc)], method = "ls" , design = design )
      fit <- eBayes(fit)
      
      res <- topTable(fit, coef = "Slope", number = nrow(data.merged), adjust="BH") #pval for the slope
      res2 <- topTable(fit, coef = 1, number = nrow(data.merged), adjust="BH")#pval for the intercept
      res3 <- topTable(fit, coef = "Quadratic", number = nrow(data.merged), adjust="BH") #pval for the quadratic term ()
      
      
      #add the pvalues to the dataframe
      
      tmp_1<- cbind(data.merged[rownames(res),],res)
      
      tmp_2<- cbind(data.merged[rownames(res2),],res2)
      
      tmp_3<- cbind(data.merged[rownames(res3),],res3)
      
      ####
      
      
      
      tobeselected <- merge.data.frame(tmp_1,tmp_2,by="Accession")
      tobeselected <- merge.data.frame(tobeselected,tmp_3,by="Accession")
      print(colnames(tobeselected))
      
      selectnames <- c(paste0(colnames(data.merged)[1:length(conc)],".x"),
                       "logFC.x", "AveExpr.x", "P.Value", "adj.P.Val", "P.Value.x",
                       "adj.P.Val.x", "P.Value.y", "adj.P.Val.y","Accession","GeneID.x","UniquePeps")
      
      
      data.merged<- tobeselected[,match(selectnames, colnames(tobeselected))]
      print(colnames(data.merged))
      proteome<-as.vector(toupper(tobeselected$GeneID.x))
      mykinases <- intersect(proteome,kinome)
      
      mykinvec <- match(proteome,mykinases, nomatch = NA)
      mykinvec[!is.na(mykinvec)]<- 'KINASE'
      
      crap <- crapome()
      # 
      tempcrap <- crap$Gene
      tempcrap<- crap$percentagePresent[match(toupper(data.merged$GeneID.x),crap$Gene)]
      data.merged <- data.frame(data.merged, CRAPomePercent = tempcrap,Kinase = mykinvec)  
      print(colnames(data.merged))

      nam <- finalNames()
      if(input$datype == 'intensity'){
        nam <- nam[-seq(1,input$reps*input$chans,by = input$chans)]
      }
      colnames(data.merged)[1:length(nam)] <- nam
     
      data.merged
    }
    
  })
  
  
  rSu<- reactive({
    if(input$reps == 1 ){
      
      
      data.merged<- dataMerge()
      data.merged
      
    }else{
      data.merged <- dataMerge2()
      data.merged
      
    }
  })
  
  
  ###################
  ## output
  ###################
  
  output$venn <-  renderPlot({
    req(data())
    if(!is.null(input$venninp)){
      
      upload <- uploadVenn()
      upload<- as.vector(toupper(upload))
      print(upload)
      data.merged <- dataMerge()
      
      proteome<-as.vector(toupper(data.merged$GeneID))
      universe <- unique(c(kinome,proteome,upload))
      
      count<- matrix(0, ncol = 3 , nrow = length(universe))
      colnames(count) <- c("kinome",'proteome','upload')
      
      for(i in 1:length(universe)){
        count[i,1]<- universe[i] %in% kinome
        count[i,2]<- universe[i] %in% proteome
        count[i,3]<- universe[i] %in% upload
      }
      
      vennDiagram(vennCounts(count), circle.col = c("blue","red","green"), cex = 1,lwd = 2)
      
      
    }else{
      
      data.merged <- dataMerge()
      
      proteome<-as.vector(toupper(data.merged$GeneID))
      universe <- unique(c(kinome,proteome))
      # Generate a matrix, with the sets in columns and possible letters on rows
      Counts <- matrix(0, nrow=length(universe), ncol=2)
      # Populate the said matrix
      for (i in 1:length(universe)) {
        Counts[i,1] <- universe[i] %in% kinome
        Counts[i,2] <- universe[i] %in% proteome
        #Counts[i,3] <- universe[i] %in% Metacore
        #Counts[i,3] <- universe[i] %in% EXOCYTOSIS
        
      }
      
      colnames(Counts) <- c("Kinome","Proteome")
      cols<-c("Red", "Blue")
      
      #### VENN
      
      vennDiagram(vennCounts(Counts), circle.col=cols,
                  cex=1, #title size
                  lwd=2 #circle line size
      )
    }
    
  })
  
  output$bar<- renderPlot({
    vec<- channels()
    vec<- length(vec)
    palette.bar <- rep(terrain.colors(input$reps),each = ifelse(input$datype != 'intensity',input$chans,input$chans - 1 ) )
    data.merged<- dataMerge()
    boxplot(data.merged[,1:vec], col=palette.bar,las=2, cex.axis=1, main=c("Box Plots"),
            ylab=c("Log2(ratios)"))
    legend("topright",legend = paste("rep",1:input$reps,sep = ' '), fill = terrain.colors(input$reps))
    
  })
  
  
  output$plot2<-renderPlot({
    
    vec <- channels()
    pal <- rainbow(length(vec))
    leg.nam <- finalNames()
    data.merged <- dataMerge()
    missing_val <- 0
    
    par(mfrow = c(1,2))   
    plot(x=rank(data.merged[,1]),y=data.merged[,1], cex.axis=1,
         main=c(paste("N. of Missing val. ",missing_val ,
                      " \n Change Distribution after LOESS",sep="")),
         col=pal[1], ylab=c("Log2(ratios)"), xlab = c("Proteins"))
    if(length(vec) > 1){
      for(i in 2:length(vec)){
        points(x=rank(data.merged[,i]),y=data.merged[,i], cex.axis=1, main=c(""),col=pal[i])
      }
    }
    legend("topleft", legend = leg.nam, fill = pal)
    
    
    plot(density(x=data.merged[,1]),col= pal[1], main=c(" Density after LOESS"), ylim = c(0,5 ),cex.axis = 1)
    legend("topleft", legend = leg.nam, fill = pal)
    
    if(length(vec) > 1){
      for(i in 2:length(vec)){
        lines(density(x=data.merged[,i]), col=pal[i], main = c(""))
      }
    }
  })
  
  
  output$plot3<- renderPlot({
    
    vec <- channels()
    
    
    data.merged<- dataMerge()
    # print(colnames(data.merged))
    meanSdPlot(as.matrix(data.merged[,1:length(vec)]))
  })
  
  output$plot4<- renderPlot({
    data.merged<- dataMerge2()
    m0 <- ggplot(data.merged, aes(x=data.merged$P.Value))
    m0<-m0 + geom_histogram(aes(fill = ..count..),binwidth = 0.01) +
      scale_fill_gradient("Count", low = "green", high = "red")+
      xlab("P.val slope")
    
    
    m1 <- ggplot(data.merged, aes(x=data.merged$P.Value.x))
    m1<- m1 + geom_histogram(aes(fill = ..count..),binwidth = 0.01) +
      scale_fill_gradient("Count", low = "green", high = "red")+
      xlab("Pval intercept")
    
    m2 <- ggplot(data.merged, aes(x=data.merged$P.Value.y))
    m2 <- m2 + geom_histogram(aes(fill = ..count..),binwidth = 0.01) +
      scale_fill_gradient("Count", low = "green", high = "red")+
      xlab("Pval quadratic")
    
    grid.arrange(m0,m1,m2)
    
  })
  
  output$volcanoint<- renderPlot({
    
    res<- dataMerge2()
    # avgthr=0.2 #sign threshold for the averege fold change 0.3(log2)  is 1.3 FC
    
    par(mar=c(5,5,5,10), xpd=TRUE)
    # Make a basic volcano plot
    with(res, plot(res$AveExpr.x, -log10(res$P.Value.x), pch=20, main="Volcano plot (Intercept pval )",xlab=c("Log2_AvgFC"),ylab=c("-Log10(Pval)"), xlim=c(-abs(max(res$AveExpr.x)+1),abs(max(res$AveExpr.x)+1))))
    
    # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
    s=subset(res, P.Value.x<input$pvalsli )
    with(s, points(s$AveExpr.x, -log10(s$P.Value.x), pch=20, col="red"))
    
    s=subset(res, abs(res$AveExpr.x)>input$avthrssli)
    with(s, points(s$AveExpr.x, -log10(s$P.Value.x), pch=20, col="orange"))
    
    s=subset(res, P.Value.x<input$pvalsli & abs(res$AveExpr.x)>input$avthrssli)
    with(s, points(s$AveExpr.x, -log10(s$P.Value.x), pch=20, col="green"))
    
    # Label points with the textxy function from the calibrate plot
    s=subset(res, P.Value.x<input$pvalsli & abs(res$AveExpr.x)>input$avthrssli)
    with(s, textxy( s$AveExpr.x, -log10(s$P.Value.x),  labs=s$GeneID.x, cex=1)
    )
    legend("bottomleft", title="Legend",cex = 0.8,
           c("Not significant",
             paste("P.Value",input$pvalsli, sep = " "), #red
             paste("AvgFC >", input$avthrssli, sep=""), #orange,
             paste("P.Value", input$pvalsli,"& AvgFC >", input$avthrssli, sep="") #green
           ),
           col=c("black","red","orange","green"),
           horiz=F, pch=c(19))
    
    
    
  })
  
  output$volcanoquad<- renderPlot({
    res<- dataMerge2()
    # avgthr=0.2 #sign threshold for the averege fold change 0.3(log2)  is 1.3 FC
    
    par(mar=c(5,5,5,10), xpd=TRUE)
    # Make a basic volcano plot
    with(res, plot(res$AveExpr.x, -log10(res$P.Value.y), pch=20, main="Volcano plot (Quadratic pval )",xlab=c("Log2_AvgFC"),ylab=c("-Log10(Pval)"), xlim=c(-abs(max(res$AveExpr.x)+1),abs(max(res$AveExpr.x)+1))))
    
    # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
    s=subset(res, P.Value.y<input$pvalsli )
    with(s, points(s$AveExpr.x, -log10(s$P.Value.y), pch=20, col="red"))
    
    s=subset(res, abs(res$AveExpr.x)>input$avthrssli)
    with(s, points(s$AveExpr.x, -log10(s$P.Value.y), pch=20, col="orange"))
    
    s=subset(res, P.Value.y<input$pvalsli & abs(res$AveExpr.x)>input$avthrssli)
    with(s, points(s$AveExpr.x, -log10(s$P.Value.y), pch=20, col="green"))
    
    # Label points with the textxy function from the calibrate plot
    s=subset(res, P.Value.y<input$pvalsli & abs(res$AveExpr.x)>input$avthrssli)
    with(s, textxy( s$AveExpr.x, -log10(s$P.Value.y),  labs=s$GeneID.x, cex=.9)
    )
    legend("bottomleft", title="Legend",cex = 0.7,
           c("Not significant",
             paste("P.Value",input$pvalsli, sep = " "), #red
             paste("AvgFC >", input$avthrssli, sep=""), #orange,
             paste("P.Value", input$pvalsli,"& AvgFC >", input$avthrssli, sep="") #green
           ),
           col=c("black","red","orange","green"),
           horiz=F, pch=c(19))
    
    
    
  })
  
  
  output$repvsrep1 <- renderUI({
    # a <- finalNames()
    # a <- paste0(a,".x")
    # b <- channels()
    a<- indexmatrix()
    selectizeInput(inputId = "repvsrep1",label = "Select Condition", choices = as.character(a$names), multiple = FALSE)
  })
  
  output$repvsrep2 <- renderUI({
    a <- finalNames()
    a <- paste0(a,".x")
    
    b <- channels()
    c <- a[b]
    if(is.null(input$repvsrep1)){
      e <- c
    } else {
      d <- match(input$repvsrep1,c)
      e <- c[-d]
    }
    
    selectInput(inputId = "repvsrep2",label = "Select Condition", choices = e, multiple = FALSE)
  })
  
  output$repvsrep <- renderPlot({
    if(input$reps > 1){
      
      data.merged <- dataMerge2()
      index <- indexmatrix()
      val = max(c(max(data.merged[,index[index[,1] == input$repvsrep1,5]],na.rm = TRUE),max(data.merged[,index[index[,1] == input$repvsrep1,6]],na.rm = TRUE)))
      
      
      plot(x=data.merged[,index[index[,1] == input$repvsrep1,5]],y=data.merged[,index[index[,1] == input$repvsrep1,6]], xlim=c(0,val +0.2), col="green" , ylim=c(0,val+0.2),
           cex.axis=1.2, main=c("LogFC-LogFC Plots"),
           xlab=paste("rep",index[index[,1] == input$repvsrep1,3]),  ylab=paste("rep",index[index[,1] == input$repvsrep1,4]))
      text(data.merged[,index[index[,1] == input$repvsrep1,5]], data.merged[,index[index[,1] == input$repvsrep1,6]], labels=data.merged$GeneID, cex= 1,pos=4)
      lines(x = c(0,val), y = c(0,val),col="red")
    } else {
      plot.new()
      title(main = 'Plot not available: only 1 replicate')
      
    }
    
  })
  
  output$plot5<- renderPlot({
    req(data())
    if(input$modtyp == 'sigmoid'){
      su <- dataMerge()
    }else{
      
      su<-  rSu()
    }
    
    index<- channels()
    index<- length(index)
    cmat <-cor(su[,1:index],use="pairwise", method = "pearson")
    corrgram(cmat, order=TRUE, lower.panel=panel.shadeNtext,
             upper.panel=panel.pie, text.panel=panel.txt,
             main="Corrgram Plots" )
    
    
  })
  
  
  ### MeanDiff
  
  output$meandiff1 <- renderUI({

    a <- indexmatrix()
    
    selectizeInput(inputId = "meandiff1",label = "Select Condition", choices = as.character(a$names), multiple = FALSE)
  })
  
  output$meandiff2 <- renderUI({
    a <- finalNames()
    a <- paste0(a,".x")
    
    b <- channels()
    c <- a[b]
    if(is.null(input$meandiff1)){
      e <- c
    } else {
      d <- match(input$meandiff1,c)
      e <- c[-d]
    }
    
    selectInput(inputId = "meandiff2",label = "Select Condition", choices = e, multiple = FALSE)
  })
  
  output$plot6<- renderPlot({

    req(data())
    if(input$modtyp == 'sigmoid'){
      plot.new()
      title(main = 'Plot not available: only 1 replicate')
    }else{
      # data.merged <- dataMerge()
      su <- rSu()
      index <- indexmatrix()
      minxy <- min(c(min(su[,index[index[,1] == input$meandiff1,5]]),min(su[,index[index[,1] == input$meandiff1,6]])))
      maxxy <- max(c(max(su[,index[index[,1] == input$meandiff1,5]]),max(su[,index[index[,1] == input$meandiff1,6]])))
      
      tmd(
        xyplot(su[,index[index[,1] == input$meandiff1,5]] ~ su[,index[index[,1] == input$meandiff1,6]]), main=input$meandiff1,xlim = c(minxy,maxxy),
        ylim = c(minxy,maxxy),
        panel=function(x, y, ...) {
          panel.xyplot(x, y, ...);
          ltext(x=x, y=y, labels=su$GeneID, pch=c(13,3,16), cex=0.5, lwd=2, pos=1, offset=1, pch = 19)
        }
      )
    }
    

  })
  
  
  output$plot7<- renderPlot({
    req(data())
    if( input$datype == 'intensity'){
      
      nchan<- input$chans - 1
      
    }else{
      
      nchan<- input$chans
    }
    
    su <- rSu()
    reps <- input$reps
    index <- channels()
    print(.libPaths())
    print(all.equal(1:nchan,index))
    print(input$reps)
    pca <- prcomp(su[,1:length(index)], scale=F)
    
    DTA<-data.frame( as.numeric(t(su[,1:length(index)])%*%pca$x[,1]),
                     as.numeric(t(su[,1:length(index)])%*%pca$x[,2]))
    
    print(DTA)
    
    p<-ggplot(DTA, aes(x=DTA$as.numeric.t.su...1.length.index........pca.x...1..,
                       y=DTA$as.numeric.t.su...1.length.index........pca.x...2..))
    

    shapeval <- c(15:18,7:12)
    p <- p + geom_point(aes(colour = factor(rep(1:print(input$reps),each = (nchan)),labels = paste("Rep",1:print(input$reps)))[index],
                            shape = factor(rep(1:nchan,print(input$reps)),labels = c(paste("C",0:(nchan-1),sep = "")))[index] ), size = 5 ) + scale_shape_manual(values=shapeval[1:nchan]) + labs(x = "PC1", y = "PC2", title="PCA") + labs(color = "Replicates", shape="Concentration")
    
    print(p)
    
  })
  
  
  
  
  output$plot8<- renderPlot({
    if( input$modtyp == 'sigmoid'){
      plot.new()
      title(main = 'Plot not available: Sigmoidal fit applied')
      
    }else{
      
      res<- dataMerge2()
      # avgthr=0.2 #sign threshold for the averege fold change 0.3(log2)  is 1.3 FC
      
      par(mar=c(5,5,5,10), xpd=TRUE)
      # Make a basic volcano plot
      with(res, plot(res$AveExpr.x, -log10(res$P.Value), pch=20, main="Volcano plot (slope pval )",
                     xlab=c("Log2_AvgFC"),ylab=c("-Log10(Pval)"),
                     xlim=c(-abs(max(res$AveExpr.x)+1),abs(max(res$AveExpr.x)+1))))
      
      # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
      s=subset(res, P.Value< input$pvalsli )
      with(s, points(s$AveExpr.x, -log10(s$P.Value), pch=20, col="red"))
      
      s=subset(res, abs(res$AveExpr.x)>input$avthrssli)
      with(s, points(s$AveExpr.x, -log10(s$P.Value), pch=20, col="orange"))
      
      s=subset(res, P.Value< input$pvalsli & abs(res$AveExpr.x)>input$avthrssli)
      with(s, points(s$AveExpr.x, -log10(s$P.Value), pch=20, col="green"))
      
      # Label points with the textxy function from the calibrate plot
      s=subset(res, P.Value< input$pvalsli & abs(res$AveExpr.x)>input$avthrssli)
      with(s, textxy( s$AveExpr.x, -log10(s$P.Value),  labs=s$GeneID, cex=1)
      )
      legend("bottomleft", title="Legend",cex = 0.7,
             c("Not significant",
               paste("P.Value",input$pvalsli, sep = " "), #red
               paste("AvgFC >", input$avthrssli, sep=""), #orange,
               paste("P.Value", input$pvalsli,"& AvgFC >", input$avthrssli, sep="") #green
             ),
             col=c("black","red","orange","green"),
             horiz=F, pch=c(19))
    }
  })
  
  output$plot9<- renderD3heatmap({
    req(data())
    if(input$modtyp == 'sigmoid'){
      su <- dataMerge()
      print(head(su))
      if(input$datype =='intensity'){
        
        su1 <- su[,1:(input$chans - 1 )]
      }else{
        su1 <- su[,1:(input$chans )]
      }
    } else {
      
      su <- rSu()
      vec <- channels()
      vec <- length(vec)
      
      su1 <- su[,1:(vec)]
      vec.nam <- finalNames()
      vec.nam <- vec.nam[channels()]
      colnames(su)[match('GeneID.x',colnames(su))] <- 'GeneID'
    }

   
    
    d3heatmap(su1, Colv = FALSE,labRow = as.character(make.names(su$GeneID,unique = TRUE)), dendrogram = 'row' )
    

  })
  
  output$test<- DT::renderDataTable({
    
    fin<- standardNames()
    
    if(length(input$view_vars) == length(fin)){
      
      DT::datatable(data.frame(Original = input$view_vars,Standard = standardNames(), Final = finalNames()))
    }else{
      DT::datatable(data.frame(Error = "Incorrect dimensions", Comments=  "ensure the nunmber of names selected is the same as the number of channels and repeats inputted"))
    }
    
    
  })
  
  ############## Info Boxes
  
  ### P value
  pvalQc <- reactive({
    
    data.merged <- dataMerge2()
    PVal <-  c(sum(data.merged$P.Value <= 0.05),sum(data.merged$P.Value.x <= 0.05),
               sum(data.merged$P.Value.y <= 0.05))
    
    Names = c("Slope", "Intercept", "Quadratic")
    
    
    
    data.frame(Names = Names, PVal = PVal)
    
  })
  
  output$infopvalslo <- renderInfoBox({
    temp <- pvalQc()
    infoBox(
      title = "Slope",paste0(temp[1,2], " p values < 0.05"),color = ifelse(temp[1,2] > 0, "green","orange" ), icon = icon(ifelse(temp[1,2] > 0, "check","warning" ))
    )
  })
  
  output$infopvalint <- renderInfoBox({
    temp <- pvalQc()
    infoBox(
      title = "Intercept",paste0(temp[2,2], " p values < 0.05"),color = ifelse(temp[2,2] > 0, "green","orange" ), icon = icon(ifelse(temp[2,2] > 0, "check","warning" ))
    )
  })
  
  output$infopvalquad <- renderInfoBox({
    temp <- pvalQc()
    infoBox(
      title = "Quadratic",paste0(temp[3,2], " p values < 0.05"),color = ifelse(temp[3,2] > 0, "green","orange" ), icon = icon(ifelse(temp[2,2] > 0, "check","warning" ))
    )
  })
  
  output$corrinfo <- renderInfoBox({
    req(data())
    if(input$modtyp == 'sigmoid'){
      su <- dataMerge()
    }else{
      su <- rSu()
      
    }
    cmat <-cor(su[,channels()],use="pairwise", method = "pearson")
    cmat[upper.tri(cmat)]<- 0
    cmat<- melt(cmat)
    
    
    infoBox(
      title = NULL,value = ifelse(nrow(cmat[cmat$value < 0, ]) == 0, "No Anti-Correlation between Channels","Anti Correlation between some Channels" ),
      color = ifelse(nrow(cmat[cmat$value < 0, ]) == 0, "green","orange" ),
      icon = icon(ifelse(nrow(cmat[cmat$value < 0, ]) == 0, "check","warning" ))
    )
    
  })
  
  output$siginfodt <- renderInfoBox({
    
    if(input$modtyp == 'lin'){
      return(NULL)
    } else {
      
      conc<- sigConc()
      if(input$datype == 'intensity'){
        
        
        top<-15 #max prot to plot
        
        data_merged_2 <- dataMerge2()
        
        # RB50<-na.omit(data_merged_2[data_merged_2$RB50Err < as.numeric(summary(data_merged_2$RB50Err)[5]) & data_merged_2$RB50Pval<0.05 & data_merged_2$predX1-data_merged_2$predX9 >0 & data_merged_2$predX1 <= 100,])
        RB50 <- data.frame(na.omit(data_merged_2[data_merged_2$RB50Err < as.numeric(summary(data_merged_2$RB50Err)[5]) & data_merged_2$RB50Pval < 0.05
                                                 & data_merged_2$predX1-data_merged_2[,paste("predX",(input$chans - 1),sep = "")] >0 & data_merged_2$predX1 <= 100,]))
        
        
        RB50_ordered<- na.omit(RB50[order(RB50$RB50Pval, decreasing = F),][1:top,])
        
        
      } else{
        pred.names <- sigPredNames()
        final.Names <- finalNames()
        top<-15 #max prot to plot
        
        data_merged_2 <- dataMerge2()
        
        # RB50<-na.omit(data_merged_2[data_merged_2$RB50Err < as.numeric(summary(data_merged_2$RB50Err)[5]) & data_merged_2$RB50Pval<0.05 & data_merged_2$predX1-data_merged_2$predX9 >0 & data_merged_2$predX1 <= 100,])
        RB50 <- data.frame(na.omit(data_merged_2[data_merged_2$RB50Err < as.numeric(summary(data_merged_2$RB50Err)[5]) & data_merged_2$RB50Pval < 0.05
                                                 & data_merged_2$predX1-data_merged_2[,paste("predX",input$chans,sep = "")] >0 & data_merged_2$predX1 <= 100,]))
        
        
        RB50_ordered<- na.omit(RB50[order(RB50$RB50Pval, decreasing = F),][1:top,])
        
      }
            
      infoBox(
        title = "RB50",paste0(nrow(RB50_ordered), " Significant RB50"),color = ifelse(nrow(RB50_ordered) > 0, "green","orange" ), icon = icon(ifelse(nrow(RB50_ordered) > 0, "check","warning" ))
      )
      
      
    }
    
  })
  
  output$siginfoslop <- renderInfoBox({
    
    if(input$modtyp == 'lin'){
      return(NULL)
    } else {
      
      conc<- sigConc()
      if(input$datype == 'intensity'){
        
        top<-15 #max prot to plot
        
        conc<- sigConc()
        pred.names <- paste0('predX',1:(input$chans -1))
        final.Names <- paste0('rep1_C',0:(input$chans - 2))
        
        data_merged_2 <- dataMerge2()
        
        #Here make the subselections for using the ggplot functions SLOPE
        slope<-na.omit(data_merged_2[data_merged_2$SlopePval<0.05 ,])
        slope_ordered<-na.omit(slope[order(slope$SlopePval, decreasing = F),][1:top,])
        
        
      } else{
        
        top<-15 #max prot to plot
        
        conc<- sigConc()
        pred.names <- sigPredNames()
        final.Names <- finalNames()
        
        data_merged_2 <- dataMerge2()
        
        #Here make the subselections for using the ggplot functions SLOPE
        slope<-na.omit(data_merged_2[data_merged_2$SlopePval<0.05 ,])
        slope_ordered<-na.omit(slope[order(slope$SlopePval, decreasing = F),][1:top,])
        
      }
      
      infoBox(
        title = "Slope Coefficient",paste0(nrow(slope_ordered), " Significant Slope"),color = ifelse(nrow(slope_ordered) > 0, "green","orange" ), icon = icon(ifelse(nrow(slope_ordered) > 0, "check","warning" ))
      )
      
      
    }
    
  })
  
  output$siginfodiff <- renderInfoBox({
    
    if(input$modtyp == 'lin'){
      return(NULL)
    } else {
      
      conc<- sigConc()
      if(input$datype == 'intensity'){
        

        data_merged_2 <- dataMerge2()
        pred.names <- paste0('predX',1:(input$chans -1))
        final.Names <- paste0('rep1_C',0:(input$chans - 2))
        
        topperc<-30 #difference in % between top and bottom
        # data_merged_2 <- dataMerge2()
        diffinter<- data_merged_2[(data_merged_2$predX1 -data_merged_2[,paste("predX",(input$chans-1),sep = "")]) > topperc & data_merged_2$predX1 <= 100, ]
        
        
        
      } else{
        conc<- sigConc()
        
        pred.names <- sigPredNames()
        final.Names <- finalNames()
        
        topperc<-30 #difference in % between top and bottom
        data_merged_2 <- dataMerge2()
        diffinter<- data_merged_2[(data_merged_2$predX1 -data_merged_2[,paste("predX",input$chans,sep = "")]) > topperc & data_merged_2$predX1 <= 100, ]
        
        
      }
      
      infoBox(
        title = "Top - Bottom Difference",paste0(nrow(diffinter), " Significant Difference"),color = ifelse(nrow(diffinter) > 0, "green","orange" ), icon = icon(ifelse(nrow(diffinter) > 0, "check","warning" ))
      )
      
      
    }
    
  })
  ############################
  # SIGMOIDAL PLOTS
  ############################ 
  
 
  output$DiffTopBottom <- renderPlot({
    
    if(input$modtyp != 'sigmoid'){
      plot.new()
      legend('topleft', c("Linear fit applied, no sigmoidal plots available"),bty = 'n')
    }else{
      if(input$datype == 'intensity'){
        data_merged_2 <- dataMerge2()
        conc<- sigConc()
        pred.names <- paste0('predX',1:(input$chans -1))
        final.Names <- paste0('rep1_C',0:(input$chans - 2))
        
        topperc<-30 #difference in % between top and bottom
        # data_merged_2 <- dataMerge2()
        diffinter<- data_merged_2[(data_merged_2$predX1 -data_merged_2[,paste("predX",(input$chans-1),sep = "")]) > topperc & data_merged_2$predX1 <= 100, ]
        
        
        if(nrow(diffinter)>0){
          Diff_Top_bottom_pred<-shape_for_ggplot_pred(diffinter,log2(conc),pred.names)
          Diff_Top_bottom_perc<-shape_for_ggplot_perc(diffinter,log2(conc),final.Names)
          what<-c("(Top - Bottom) >")
          
          Diff_Top_bottom<-ggplot()+
            geom_line(data = Diff_Top_bottom_pred, aes(x=x,y=value, colour=factor(Diff_Top_bottom_pred$GeneID)), size = 1) +
            geom_point(data = Diff_Top_bottom_perc, aes(x=x,y=value, colour=Diff_Top_bottom_perc$GeneID)) +
            labs(title=paste(what,topperc,sep=""))
          
          
        }else{
          Diff_Top_bottom<-ggplot()+
            labs(title=paste("No significant Top-Bottom >" ,topperc,"%","\n","has been found", sep=""))
        }
        
        print(Diff_Top_bottom)
      }
      
      else{
        conc<- sigConc()
        
        pred.names <- sigPredNames()
        final.Names <- finalNames()
        
        topperc<-30 #difference in % between top and bottom
        data_merged_2 <- dataMerge2()
        diffinter<- data_merged_2[(data_merged_2$predX1 -data_merged_2[,paste("predX",input$chans,sep = "")]) > topperc & data_merged_2$predX1 <= 100, ]
        
        
        
        
        if(nrow(diffinter)>0){
          Diff_Top_bottom_pred<-shape_for_ggplot_pred(diffinter,log2(conc),pred.names)
          Diff_Top_bottom_perc<-shape_for_ggplot_perc(diffinter,log2(conc),final.Names)
          what<-c("(Top - Bottom) >")
          
          Diff_Top_bottom<-ggplot()+
            geom_line(data = Diff_Top_bottom_pred, aes(x=x,y=value, colour=factor(Diff_Top_bottom_pred$GeneID)), size = 1) +
            geom_point(data = Diff_Top_bottom_perc, aes(x=x,y=value, colour=Diff_Top_bottom_perc$GeneID)) +
            labs(title=paste(what,topperc,sep=""))
          
          
        }else{
          Diff_Top_bottom<-ggplot()+
            labs(title=paste("No significant Top-Bottom >" ,topperc,"%","\n","has been found", sep=""))
        }
        
        print(Diff_Top_bottom)
      }
    }
    
    
  })
  
  output$Slope_pl <- renderPlot({
    
    if(input$modtyp != 'sigmoid'){
      plot.new()
      legend('topleft', c("Linear fit applied, no sigmoidal plots available"),bty = 'n')
    }else{
      
      if(input$datype == 'intensity'){
        top<-15 #max prot to plot
        
        conc<- sigConc()
        pred.names <- paste0('predX',1:(input$chans -1))
        final.Names <- paste0('rep1_C',0:(input$chans - 2))
        
        data_merged_2 <- dataMerge2()
        
        #Here make the subselections for using the ggplot functions SLOPE
        slope<-na.omit(data_merged_2[data_merged_2$SlopePval<0.05 ,])
        slope_ordered<-na.omit(slope[order(slope$SlopePval, decreasing = F),][1:top,])
        if(nrow(slope_ordered)>0){
          slope_pred<-shape_for_ggplot_pred(slope_ordered,log10(conc),pred.names)
          slope_perc<- shape_for_ggplot_perc(slope_ordered,log10(conc),final.Names)
          what<-c("Slope (p.val) ")
          Slope_pl<-ggplot()+
            geom_line(data = slope_pred, aes(x=x,y=value, colour=factor(slope_pred$GeneID)), size = 1) +
            geom_point(data = slope_perc, aes(x=x,y=value,colour=slope_perc$GeneID))+
            labs(title=paste(what,"Top",top,sep=""))
          
        }else{Slope_pl<-ggplot()+
          labs(title="No significant Sigmoidal Slope has been found")
        }
        print(Slope_pl)
      } else{ 
        top<-15 #max prot to plot
        
        conc<- sigConc()
        pred.names <- sigPredNames()
        final.Names <- finalNames()
        
        data_merged_2 <- dataMerge2()
        
        #Here make the subselections for using the ggplot functions SLOPE
        slope<-na.omit(data_merged_2[data_merged_2$SlopePval<0.05 ,])
        slope_ordered<-na.omit(slope[order(slope$SlopePval, decreasing = F),][1:top,])
        if(nrow(slope_ordered)>0){
          slope_pred<-shape_for_ggplot_pred(slope_ordered,log10(conc),pred.names)
          slope_perc<- shape_for_ggplot_perc(slope_ordered,log10(conc),final.Names)
          what<-c("Slope (p.val) ")
          Slope_pl<-ggplot()+
            geom_line(data = slope_pred, aes(x=x,y=value, colour=factor(slope_pred$GeneID)), size = 1) +
            geom_point(data = slope_perc, aes(x=x,y=value,colour=slope_perc$GeneID))+
            labs(title=paste(what,"Top",top,sep=""))
          
        }else{Slope_pl<-ggplot()+
          labs(title="No significant Sigmoidal Slope has been found")
        }
        print(Slope_pl)
      }
    }
    

  })
  
  output$RB50 <- renderPlot({
    
    if(input$modtyp != 'sigmoid'){
      plot.new()
      legend('topleft', c("Linear fit applied, no sigmoidal plots available"),bty = 'n')
    }else{
      if(input$datype =='intensity'){
        conc<- sigConc()
        
        pred.names <- paste0('predX',1:(input$chans -1))
        final.Names <- paste0('rep1_C',0:(input$chans - 2))
        
        top<-15 #max prot to plot
        
        data_merged_2 <- dataMerge2()
        
        # RB50<-na.omit(data_merged_2[data_merged_2$RB50Err < as.numeric(summary(data_merged_2$RB50Err)[5]) & data_merged_2$RB50Pval<0.05 & data_merged_2$predX1-data_merged_2$predX9 >0 & data_merged_2$predX1 <= 100,])
        RB50 <- data.frame(na.omit(data_merged_2[data_merged_2$RB50Err < as.numeric(summary(data_merged_2$RB50Err)[5]) & data_merged_2$RB50Pval < 0.05
                                                 & data_merged_2$predX1-data_merged_2[,paste0('predX',(input$chans - 1))] >0 & data_merged_2$predX1 <= 100,]))
        
        
        RB50_ordered<- na.omit(RB50[order(RB50$RB50Pval, decreasing = F),][1:top,])
        
        if(nrow(RB50_ordered)>0){
          RB50_pred<-shape_for_ggplot_pred(RB50_ordered,log10(conc),pred.names)
          RB50_perc<-shape_for_ggplot_perc(RB50_ordered,log10(conc),final.Names)
          what<-c("RB50 (p.val) ")
          RB50_pl<-ggplot()+
            geom_line(data = RB50_pred, aes(x=x,y=value, colour=factor(RB50_pred$GeneID)), size = 1) +
            geom_point(data = RB50_perc, aes(x=x,y=value,colour=RB50_perc$GeneID))+
            labs(title=paste(what,"Top",top,sep=""))
          print(RB50_pl)
        }else{
          RB50_pl<-ggplot()+
            labs(title="No significant RB50 has been found")
        print(RB50_pl)
        }
        
      }else {
        conc<- sigConc()
        pred.names <- sigPredNames()
        final.Names <- finalNames()
        top<-15 #max prot to plot
        
        data_merged_2 <- dataMerge2()
        
        # RB50<-na.omit(data_merged_2[data_merged_2$RB50Err < as.numeric(summary(data_merged_2$RB50Err)[5]) & data_merged_2$RB50Pval<0.05 & data_merged_2$predX1-data_merged_2$predX9 >0 & data_merged_2$predX1 <= 100,])
        RB50 <- data.frame(na.omit(data_merged_2[data_merged_2$RB50Err < as.numeric(summary(data_merged_2$RB50Err)[5]) & data_merged_2$RB50Pval < 0.05
                                                 & data_merged_2$predX1-data_merged_2$predX9 >0 & data_merged_2$predX1 <= 100,]))
        
        
        RB50_ordered<- na.omit(RB50[order(RB50$RB50Pval, decreasing = F),][1:top,])
        
        if(nrow(RB50_ordered)>0){
          RB50_pred<-shape_for_ggplot_pred(RB50_ordered,log10(conc),pred.names)
          RB50_perc<-shape_for_ggplot_perc(RB50_ordered,log10(conc),final.Names)
          what<-c("RB50 (p.val) ")
          RB50_pl<-ggplot()+
            geom_line(data = RB50_pred, aes(x=x,y=value, colour=factor(RB50_pred$GeneID)), size = 1) +
            geom_point(data = RB50_perc, aes(x=x,y=value,colour=RB50_perc$GeneID))+
            labs(title=paste(what,"Top",top,sep=""))
        }else{
          RB50_pl<-ggplot()+
            labs(title="No significant RB50 has been found")
        }
        print(RB50_pl)
      }
    }

  })
  

  
  # Linear elements of sigmoidal obselete 
  
  # output$Linear_pl1 <- renderPlot({
  #   
  #   if(input$modtyp != 'sigmoid'){
  #     plot.new()
  #     legend('topleft', c("Linear fit applied, no sigmoidal plots available"),bty = 'n')
  #   }else{
  #     conc<- sigConc()
  #     
  #     pred.names <- sigPredNames()
  #     final.Names <- finalNames()
  #     top<-15 #max prot to plot
  #     
  #     data_merged_2 <- dataMerge2()
  #     
  #     
  #     linere<-data_merged_2[data_merged_2$RB50Pval=="lm-fit:intercept.slope.quadratic" & data_merged_2$predX1 > data_merged_2[,paste("predX",input$chans,sep = "")] & data_merged_2$predX1 <= 100 & data_merged_2$Lower_LimitPval < 0.05,]
  #     
  #     if(nrow(linere)>0){
  #       lin_pred<-shape_for_ggplot_pred(linere,log10(conc),pred.names)
  #       lin_perc<-shape_for_ggplot_perc(linere,log10(conc),final.Names)
  #       what<-c("Lm model Slope Pval< 0.05 \n (failed dose-resp.) ")
  #       
  #       Linear_pl1<-ggplot()+
  #         geom_line(data = lin_pred, aes(x=x,y=value, colour=factor(lin_pred$GeneID)), size = 1) +
  #         geom_point(data = lin_perc, aes(x=x,y=value,colour=lin_perc$GeneID))+
  #         labs(title=paste(what,top,sep=""))
  #       
  #       
  #     }else{
  #       Linear_pl1<-ggplot()+
  #         labs(title="No significant Lm slope has been found")
  #     }
  #     
  #     print(Linear_pl1)
  #   }
  # })
  # 
  # output$Linear_pl2 <- renderPlot({
  #   
  #   if(input$modtyp != 'sigmoid'){
  #     plot.new()
  #     legend('topleft', c("Linear fit applied, no sigmoidal plots available"),bty = 'n')
  #   }else{
  #     conc<- sigConc()
  #     pred.names <- sigPredNames()
  #     final.Names <- finalNames()
  #     top<-15 #max prot to plot
  #     
  #     data_merged_2 <- dataMerge2()
  #     
  #     
  #     linere<-data_merged_2[data_merged_2$RB50Pval=="lm-fit:intercept.slope.quadratic" & data_merged_2$predX1 > data_merged_2[,paste("predX",input$chans,sep = "")] & data_merged_2$predX1 <= 100 & data_merged_2$Upper_LimitPval < 0.05,]
  #     
  #     if(nrow(linere)>0){
  #       lin_pred<-shape_for_ggplot_pred(linere,log10(conc),pred.names)
  #       lin_perc<-shape_for_ggplot_perc(linere,log10(conc),final.Names)
  #       what<-c("Lm model Quadratic Pval< 0.05 \n (failed dose-resp.) ")
  #       
  #       Linear_pl2<-ggplot()+
  #         geom_line(data = lin_pred, aes(x=x,y=value, colour=factor(lin_pred$GeneID)), size = 1) +
  #         geom_point(data = lin_perc, aes(x=x,y=value,colour=lin_perc$GeneID))+
  #         labs(title=paste(what,top,sep=""))
  #       
  #       
  #     }else{
  #       Linear_pl2<-ggplot()+
  #         labs(title="No significant Lm Quadratic has been found")
  #     }
  #     
  #     
  #     print(Linear_pl2)
  #   }
  #   
  # })
  # 
  output$testmerge <- DT::renderDataTable({
    
    # a<- intData()
    a<- dataMerge2()
    # a <- indexmatrix()
    
    # a <- rSu()
    
    
    
    if(input$modtyp == 'sigmoid'){
      if(input$incpd == TRUE){
        DT::datatable(data.frame(GeneID = a$GeneID, RB50 = a$RB50Coef, RB50pval = a$RB50Pval,Topminusbottom = a$Top_minus_min ,correctedRB50 = a$correctedRB50, depletionConst = a$depletionConstant, Kinase = a$Kinase),
                      options = list(scrollX = TRUE)  )
      }else{
        DT::datatable(data.frame(GeneID = a$GeneID, RB50 = a$RB50Coef, RB50pval = a$RB50Pval,Topminusbottom = a$Top_minus_min , Kinase = a$Kinase),
                      options = list(scrollX = TRUE)  )
      }
    
    } else{
      DT::datatable(data.frame(GeneID = a$GeneID.x, Intercept = signif(a$P.Value), Slope = signif(a$P.Value.x), Quadtratic = signif(a$P.Value.y), Kinase = a$Kinase),
                    options = list(scrollX = TRUE)  )
    }
    
    
    
  })
  
  
  output$kintab <- DT::renderDataTable({
    data.merged <- dataMerge()
    
    proteome<-as.vector(toupper(data.merged$GeneID))
    DT::datatable(data.frame(GeneID = intersect(proteome,kinome)))
  }
  
  )
  
  #################
  ## OUTPUT download
  #################
  
  output$peprmv<- downloadHandler(
    filename = function() {
      paste("removedPep", '.csv', sep='')
    },
    content = function(file) {
      write.csv(pepdwn(), file)
    }
  )
  
  output$report <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "report.html",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "report.Rmd")
      file.copy("report.Rmd", tempReport, overwrite = TRUE)
      
      # Set up parameters to pass to Rmd document
      params <- list(chans = input$chans, reps = input$reps, data = dataMerge(), data2 = dataMerge2(), channel = channels(), avthrsli = input$avthrssli,
                     pvalsli = input$pvalsli, indexmat = indexmatrix(), RSu = rSu(), finNam = finalNames(), datype = input$datype, sigmodin = input$modtyp, concen = sigConc(), sigPred = sigPredNames()
                      ,
                      vennip = input$venninput, kin = kinome1
                      # , upVenn <- uploadVenn()
                     )
      
      
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
  
  output$downloadData <- downloadHandler(
    
    filename = function() {
      paste(input$dataset, '.csv', sep='')
    },
    content = function(file) {
      write.csv(dataMerge2(), file)
    }
  )
  
  ## rest elements -- doesn't get rendered 
  
  output$testkd<- DT::renderDataTable({
    
    # kd<- rKd()
    # 
    # DT::datatable(data.frame(kd))
      # a<- dataMerge2()
       # a<- rSu()
      a<- intData()
    DT::datatable(a,
                  options = list(scrollX = TRUE) )
  })
  
  
  
})

# Run the application


shinyApp(ui = ui, server = server)