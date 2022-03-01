#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyTree)
library(RSQLite)
library(rjson)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(DT)
library(Rtsne)
library(gridExtra)
library(grid)
library(tidyverse)

#File paths:
webApp_data<-'./Data/webApp_data/'
webApp_functions<-'./Data/webApp_functions/'

#Initialize TCGA_DB:
source(file.path(webApp_functions,'TCGA_DB_functions/TCGA_DB_connect.R'))
message('TCGA_DB connection established')
#Loading "webApp_data":
#GlycoOnto structure:
glycoOnto_hierarchy<-fromJSON(file=file.path(webApp_data,'glycoOnto_hierarchy.json'))
#Pathway and Function gene lists:
#Loads "pathwayList", "functionList", and various coloring R objects:
for (f in list.files(file.path(webApp_data),pattern='*.rda',full.names=T)){
	load(f)
}
#Source webApp functions:
for (f in list.files(file.path(webApp_functions),pattern='*.R',recursive=T,full.names=T)){
	source(f)
}
SNFG_functionColors<-read.table(file='./Data/webApp_data/glycoOnto_functionColors.tsv',sep='\t',header=T) %>% 
	mutate(color=sapply(color,function(x) sub('^\\ ','',x)))
message('Finished sourcing dependencies')

print(glycoOnto_hierarchy)

#############################################
# TCGA Glycogene Explorer WebApp Server START
#############################################

shinyServer(function(input, output,session) {
  
  ######################################
  # t-SNE Logical Functions:
  ######################################
  #Metadata Select UI
  output$metaoption_ui<-renderUI({ctype_label_level_select('meta_level',input$tSNE_ctype,TCGA_DB)})
  message("Rendering metaUi")
  output$tSNE_metaSelect_UI<-renderUI({
	if (is.null(input$tSNE_ctype)){
		return(NULL)
	}
	else {
		ctype_subtype_select(input$tSNE_ctype,TCGA_DB)
	}
	  
  })
  
  #Glycogene set selection UI
  output$tSNE_glycoOnto<-shinyTree::renderTree({glycoOnto_hierarchy})

  #Making tSNE plot based on selections:
  output$cancer_lab_tsne<-renderPlot({
    #First gather data based on inputs:
    # "TCGA_get_tSNE_data gathers information
	  
    ### Wait for "calc" button to be pressed:
    input$calc_button

    ### Gene set selection:
    # Snippet below looks at which ontology groups 
    #  are selected and parses the glycogenes involved
    tree<-input$tSNE_glycoOnto
    req(tree)
    selected_paths<-get_selected(tree)
    geneSet<-get_geneList(selected_paths,pathwayList,functionList)

    ### Data parsing:
    # If there is no selection for a metadata level or there are 
    # not glycogenes selected, return NULL for the plot.
    # Otherwise, run the database "get" method for read counts
    #Convert ctype to pid:
    tSNE_ctype_pid<-get_ctype_project(input$tSNE_ctype,TCGA_DB)
    if (input$meta_level==2){
      #Ensure an actual selection was made:
      if (is.null(input$metaChoice) | length(geneSet)==0){
        return(NULL)
      } else {
        #Parse data 
        dta<-TCGA_get_tSNE_data(tSNE_ctype_pid,input$metaChoice,geneSet,TCGA_DB)
        ### Remove NA columns:
        dta<-dta %>% dplyr::filter(!is.na(subtype_level))
        labs<-dta %>% pull(subtype_level)
      }
    } else if (input$meta_level==1){
      if (length(geneSet)==0){
        return(NULL)
      } else {
        dta<-TCGA_get_tSNE_data(tSNE_ctype_pid,'shortLetterCode',geneSet,TCGA_DB)
        labs<-dta %>% pull(shortLetterCode)
      }
    }

    ### tSNE parameters:
    #Perplexity:
    if (as.numeric(input$perp_select)==1){
      perp<-dim(dta)[1]/5
    }
    else if (as.numeric(input$perp_select)==2){
      perp<-as.numeric(input$perp_num)
    }
    #PCA initial dims:
    if (is.null(input$PCA_initdim)==FALSE){
      initial_dims=as.numeric(input$PCA_initdim)
    }
    else {
      if (is.null(dim(dta))){
        #Means dta is 1d
        initial_dims=1
      } else {
      initial_dims=dim(dta)[2]/0.25
      }
    }
    ### 
    
    ### Running tSNE algorithm on gene count data:
    tsne_data<-Rtsne(X=as.data.frame(dta) %>% select(one_of(geneSet)) %>% as.matrix(),
		     dims=2,
		     initial_dims = initial_dims,
		     perplexity = perp,
                     max_iter = as.numeric(input$max_iter),
		     verbose = TRUE)

    ### Put tSNE results in data frame:
    tsne_df<-data.frame(x=tsne_data$Y[,1],
			y=tsne_data$Y[,2],
			label=labs)
    ### Create Plot:
    g<-ggplot(tsne_df) + geom_point(aes(x=x,y=y,colour=label),size=3) + labs(x='tSNE 1',y='tSNE 2',title=paste(as.character(input$cancer_type))) +
      theme_light() + 
      scale_color_discrete(name='Tissue Type') +
      theme(axis.text.x = element_text(size = 18),axis.text.y = element_text(size=20),
            axis.title.y = element_text(size=25),axis.title.x=element_text(size=25),legend.title = element_text(size=22.5),
            legend.text=element_text(size=18.5), legend.key.size = unit(0.65,units = 'cm'),title=element_text(size=25))
    return(g)
  })

  ######################################
  # DE BoxPlot Logical Functions:
  ######################################

  #Metadata Select UI
  #output$DE_metaoption_ui<-renderUI({ctype_label_level_select('DE_metaSelect_UI',input$DE_ctype,TCGA_DB)})
  
  #Glycogene set selection UI
  output$DE_glycoOnto<-glycoOnto_tree(glycoOnto_hierarchy)

  geneSelections<-reactiveVal(c())
  #Accumulating gene list to build custom gene lists
  observeEvent(input$DE_geneSet,{
  	newList<-c(geneSelections(),input$geneSelection)
  	geneSelections(unique(newList))
   })
  
  observeEvent(input$DE_cleargenes_button,{
  	geneSelections(c())
   })

  # Custom geneList UI:
  output$DE_geneSelectUI<-renderUI({
  	tree<-input$DE_glycoOnto
  	req(tree)
  	selected<-get_selected(tree)
  	#Get genelist:
  	geneList<-get_geneList(selected,pathwayList,functionList)
  	#Get already-selected genes:
  	preSelect_geneList<-geneSelections()[geneSelections() %in% geneList]
  	if (length(preSelect_geneList)>0){
  		selectInput("DE_geneSet","Select Genes:",choices=geneList,multiple=T,selected=preSelect_geneList)
  	} else {
  		selectInput("DE_geneSet","Select Genes:",choices=geneList,multiple=T)
  	}
  })
  output$DE_cleargene_button<-renderUI({
	  actionButton("de_cleargenes","Clear Genes")
  })

  ### Differential Expression barplot:
  output$deBarPlot<-renderPlot({
    # load('./Data/Glycogene_types_list.rda')
    #Process selections
    # print('Beginning')
    geneSet<-input$DE_geneSet
    ctypes<-input$DE_ctypes
    if (any(is.null(c(geneSet,ctypes)))){
	return(NULL)
    } else {
        DE_ctype_pids<-sapply(input$DE_ctypes,function(x) get_ctype_project(x,TCGA_DB))
	dta<-get_diffExp_data(DE_ctype_pids,geneSet,TCGA_DB)
    }
    #Map genes and colors into dataframe:
    dta<-SNFG_map_main(dta,functionList,SNFG_functionColors)
    #dta$gene_type<-sapply(dta$geneName,function(g) get_ggene_functions(g,Glycogene_types_list))
    #dta$gene_color<-sapply(dta$geneName, function(ggene) SNFG_mapping_function(ggene,functionList,SNFG_functionColors))
    #Order the data by gene name, then by function:
    #dta<-dta[order(dta$geneName,dta$Function),]
    # Create plot labels where the genes have functions in them:
    dta$labels<-dta$geneName
    # print('About To Plot...')
    p<-plot_scatter_DE_plot(dta,SNFG_functionColors)
    return(p)
  },height = 850,width=1150)


  ############################################
  # Glycogene Distribution Logical Functions:
  ############################################

  # Glycogene set selection UI
  output$dist_glycoOnto<-glycoOnto_tree(glycoOnto_hierarchy)
  # Custom geneList UI:
  output$dist_geneSelectUI<-renderUI({
  	tree<-input$dist_glycoOnto
  	req(tree)
  	selected<-get_selected(tree)
  	#Get genelist:
  	geneList<-get_geneList(selected,pathwayList,functionList)
	#Only one option allowed, return select:
	selectInput("dist_gene","Select Glycogene:",choices=geneList,multiple=F)
  })

  output$glycoDist<-renderPlot({
    geneSelection<-input$dist_gene
    plotSelection<-input$plt_Type
    if (is.null(geneSelection) | is.null(plotSelection)){
      return(NULL)
    } else {
	    #Gather VST count data for selected gene:
	    dta<-TCGA_get_geneCount_data(geneSelection,TCGA_DB)
    	    plot_gene_dists(dta,plotSelection)
    }
  },height = 850)

  ###############################
  # GlycoPathway Enrichment Logic
  ###############################
  
  #"enrichTabset" renders glycogene pathway and molecular function
  # enrichments in a tabset, where each tab is broken into different 
  # differential expression comparisions
  #  "build_enrich_tabset" is a function which constructs hierarhcial tabs
  #  which render different hierarchical trees of enri
  #  Arguments:
  #  enrichment_ctype: the cancer type selected to view enrichment results:
  #output$enrichmentTabset<-renderUI({enrichment_vars_tabset(input$enrich_ctype,TCGA_DB,build_enrich_tabset,glycoOnto_hierarchy)})


  #enrichment_ctype<-input$enrich_ctype
  output$subtypeSelect<-renderUI({
	  ctype_pid<-get_ctype_project(input$enrich_ctype,TCGA_DB)
	  selectInput("enrich_subtype_select","Select Subtype",choices=c('shortLetterCode',get_subtypes(ctype_pid,TCGA_DB)),multiple=F)
  })
  output$levelSelect<-renderUI({
	ctype_pid<-get_ctype_project(input$enrich_ctype,TCGA_DB)
	subtype_levels<-get_subtype_levels(ctype_pid,input$enrich_subtype_select,TCGA_DB)
	selectInput("enrich_subtype_level_select","Select Subtype Level",choices=subtype_levels,multiple=F)
  })

  output$universeSelect<-renderUI({
	radioButtons('universe_choice',label='Pick Universe',
		     choices=list('Glycogenes'='glycogene','Genomic'='genome')
		     )
  })
  output$directionSelect<-renderUI({
	radioButtons('direction_choice',label='Pick Direction',
		     choices=list('Upregulated'='up','Downregulated'='down')
		     )
  })
  output$enrichTreeOut<-renderTree({
	  #Search for input variables rendered by "enrichment_vars_tabset"
	  # to render enrichment statistics in a tree:
	  ctype<-input$enrich_ctype
          ctype_pid<-get_ctype_project(input$enrich_ctype,TCGA_DB)
	  if (any(sapply(c(input$enrich_subtype_level_select,input$universe_choice,input$direction_choice),is.null))){
		  return(NULL)
	  } else {
		  #Debugging:
		  #Gather Expression Data:
		  message(input$enrich_subtype_select)
		  if (input$enrich_subtype_select=='shortLetterCode'){
			  modSubtype<-sub(input$enrich_subtype_select,'',input$enrich_subtype_level_select)
		  } else {
			  modSubtype<-input$enrich_subtype_level_select
		  }
		  message(modSubtype)
		  enrich_dta<-get_enrichment_dta(ctype_pid,modSubtype,input$universe_choice,input$direction_choice,TCGA_DB)
		  message(enrich_dta)
		  #Build the Relabeled Enrichment Tree:
		  glycoOnto_enrich_hierarchy<-glycoOnto_hierarchy[names(glycoOnto_hierarchy)!='Compartment']
		  enrich_tree_stats<-build_enrich_tree(glycoOnto_enrich_hierarchy,enrich_dta)
		  return(enrich_tree_stats)
	  }
  })
  

  ###################################
  # END GlycoPathway Enrichment Logic
  ###################################
})
