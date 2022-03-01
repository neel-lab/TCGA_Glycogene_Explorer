#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyTree)
library(DT)
#File paths:
webApp_data<-'./Data/webApp_data/'
webApp_functions<-'./Data/webApp_functions/'
for (f in list.files(webApp_functions,pattern='*.R',full.names=T)){
	source(f)
}
#Source Database:
source(file.path(webApp_functions,'TCGA_DB_functions/TCGA_DB_connect.R'))
# Define UI for application that draws a histogram

shinyUI(bootstrapPage(navbarPage("TCGA Glycogene Explorer",
    
  # Application title
   #titlePanel("TCGA Glycogene Explorer"),

  ######################################
  # t-SNE UI:
  ######################################
  
  tabPanel('Glycogene tSNE Plot',
           sidebarPanel(
             helpText('Explore TCGA data sets with t-SNE plots'), width = 3,
             
             #Select the Cancer Type:
	     ctype_name_select('tSNE_ctype',TCGA_DB,multiSelect=FALSE),
             #Select level for plotting
	     uiOutput('metaoption_ui'),
             #radioButtons('meta_level',label='Select Data Labeling',
             #             choices = list('Tumor vs Normal'=1,'Clinical Metadata'=2),
             #             selected = 1),
	     conditionalPanel(
			   'input.meta_level==2',
			    uiOutput('tSNE_metaSelect_UI')
	     ),
             #Select the Gene List:
	     #GlycoOnto_UI, just groups:
	     #Render GlycoOnto tree UI:
	     #shinyTree("tSNE_glycoOnto", stripes = TRUE, multiple = TRUE, animation = FALSE,checkbox=TRUE,themeIcons=FALSE,themeDots=FALSE),
	     strong('Select Glycogene set:'),
	     glycoOnto_selectUI("tSNE_glycoOnto",singleSelect=TRUE,checkbox=TRUE),
	     #glycoOnto_selectUI("tSNE_glycoOnto",singleSelect=TRUE,checkbox=TRUE),
	     #t-SNE parameters:
             strong('t-SNE Parameters'),
             br(),br(),
             radioButtons('perp_select',label = 'Perplexity Method',
                          choices = list('Calculate'=1,'Specify'=2),
                          selected = 1),
             conditionalPanel(condition = "input.perp_select==2",
                              numericInput('perp_num',label='Perplexity',value=0)),
             numericInput('max_iter',label = 'Iterations',value = 1000),
             radioButtons('PCA_dim',label = 'Initial PCA Dimensions',
                          choices = list('Calculate'=1,'Specify'=2)),
             conditionalPanel(condition = 'input.PCA_dim==2',
                              numericInput('PCA_initdim',label='PCA Initial Dimension',
                                           value=20,min=2)),
             actionButton('calc_button',label = 'Run tSNE'),br(),
             strong('About Glycogene tSNE Viewer:'),br(),
             p('For any user-defined glycogene class, users can create tSNE visualizations of cancer samples for
               analyzed cancer types.  Users can modify tSNE computation parameters, such as iteration number,
               perplexity (cluster density parameter), as well as the intermediate PCA dimensions used before 
               embedding in 2D space.  Users have the choice of either labeling the data at the \"Tumor vs Normal\"
               level, or can select clinical metadata from the TCGA.')
             
             
           ),mainPanel(plotOutput('cancer_lab_tsne',height = 850))
           ),
  ######################################
  # Differential Expression UI:
  ######################################
  tabPanel("DE Viewer",
    sidebarPanel(width=2,
        #Select the Cancer Type:
	ctype_name_select('DE_ctypes',TCGA_DB,multiSelect=TRUE),
	#uiOutput('DE_metaoption_ui'),
        glycoOnto_selectUI("DE_glycoOnto",singleSelect=TRUE,checkbox=TRUE),
	uiOutput('DE_geneSelectUI'),
	uiOutput('DE_cleargene_button'),
      strong('Analysis: '),
      p('The log2 fold change statistics plotted were computed using DESeq2, comparing the differences in 
        glycogene expression in tumor samples with respect to normal. The significance of the log2 fold 
        change values has undergone multiple testing correction using the Benjamini Hochberg method.')
    ),fluidRow(
      shiny::column(3,offset=2,
                    plotOutput('deBarPlot'))
    )
  ),
  ######################################
  # Glycogene Distribution UI:
  ######################################

  tabPanel('Glycogene Scatterplot Viewer',
           sidebarPanel(width=2,
		glycoOnto_selectUI("dist_glycoOnto",singleSelect=TRUE,checkbox=TRUE),
		uiOutput('dist_geneSelectUI'),
             shiny::radioButtons('plt_Type',label='Select Plot Type',
                                choices=c('Binned Bar','Scatter Plot'),selected = 'Scatter Plot'),
             p('Glycogene VST data from TCGA datasets for different cancers presented either as a scatter plot 
               with each dot representing one patient, or binned bars where every cancer type is discretized 
               into 300 bins.')
           ),
           mainPanel(
             
             plotOutput('glycoDist')
             
           )
        ),
  ######################################
  # GlycoEnzOnto Enrichment UI:
  ######################################
  tabPanel('Dysregulated Glycosylation Pathways',
           sidebarPanel(width=3,

             #Select the Cancer Type:
	     ctype_name_select('enrich_ctype',TCGA_DB,multiSelect=FALSE),
	     #Render the tabset of enrichment data for 
		# the given cancer type here:
		uiOutput("subtypeSelect"),
	        conditionalPanel(
	           	   "input.enrich_subtype_select!=''",
	           	    uiOutput('levelSelect')
	        ),
		conditionalPanel(
			   "input.enrich_subtype_level_select!=''",
			   uiOutput('universeSelect'),
		 ),
		conditionalPanel(
			   "input.enrich_subtype_level_select!=''",
			   uiOutput('directionSelect')
			   ),
             p('Significantly dysregulated glycogenes for each cancer type are mapped to specific pathways. 
               The radar plots show the -log10(pvalues) of pathways that are either up- or downregulated in each pathway, as determined by Fisher\'s Exact Test.
               Any radar points having a value extending beyond the grey circle indicates that pathway is enriched with a p-value less than or equal to 0.05.
               All pathay genes\' differential expression statistics are viewable in the tables below.  Significantly
               up and downregulated genes are highlighted in red and blue, respectively.')
           ),
           mainPanel(
		static_tree_ui("enrichTreeOut",singleSelect=F,checkbox=F)
		#static_tree_ui("enrichTreeOut",singleSelect=FALSE,checkbox=FALSE)
             )
           )
  
)))
