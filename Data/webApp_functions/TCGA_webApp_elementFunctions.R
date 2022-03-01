library(shiny)
library(rjson)
library(dplyr)
source('./Data/webApp_functions/TCGA_DB_functions/count_data_getters.R')
source('./Data/webApp_functions/TCGA_DB_functions/ctype_name_getters.R')
glycoOntoStruct<-fromJSON(file='./Data/webApp_data/glycoOnto_hierarchy.json')
#Pathway and molecular function data:
load('./Data/webApp_data/pathwayList.rda')
load('./Data/webApp_data/functionList.rda')

#### Methods to generate shiny objects within WebApp:

#Full Cancer Name select UI:
ctype_name_select<-function(nme,TCGA_DB,multiSelect){
	#Run database query to get ctype names:
	ctypes<-dbGetQuery(TCGA_DB,
			   "SELECT DISTINCT TCGA_Diseases.disease_type FROM TCGA_Diseases WHERE TCGA_Diseases.project_id IN (SELECT DISTINCT ctype FROM TCGA_count_data)") %>% pull(disease_type)
	#Create selectInput element, and return:
	si<-selectInput(nme,
		 label = 'Select Cancer Type',
		 choices = ctypes,
		 multiple=multiSelect #True or False
	)
	return(si)
}

#Conditional radio button select UI for "tumor vs normal" and 
# cancer metadata labeling of tSNEs:
ctype_label_level_select<-function(nme,ctype_select,TCGA_DB){
	#This function returns radio button selection 
	# based on whether a ctype has any valid subtypes
	ctype_pid<-get_ctype_project(ctype_select,TCGA_DB)
	subtype_list<-get_subtypes(ctype_pid,TCGA_DB)
	if (is.na(subtype_list)){
		rb<-NULL
	} else {
		rb<-radioButtons(nme,label='Select Data Labeling',
			  choices = list('Tumor vs Normal'=1,'Clinical Metadata'=2),
			  selected = 1)
	}
	return(rb)
}

#Render a select UI when cancer metadata is available 
ctype_subtype_select<-function(ctype,TCGA_DB){
	#Renders UI for selecting cancer subtypes
	#Gather subtypes:
	ctype_pid<-get_ctype_project(ctype,TCGA_DB)
	ctype_subtype<-get_subtypes(ctype_pid,TCGA_DB)
	#Render selection ui
	selectInput(inputId = 'metaChoice',label='Select Metadata Variable',
		choices = ctype_subtype)
}


### GlycoOntology UI render:
# This function must run server-side
glycoOnto_selectUI<-function(tree_outputName,singleSelect=TRUE,checkbox=TRUE){
	treeUI<-shinyTree(tree_outputName, stripes = TRUE, multiple = singleSelect, animation = FALSE,checkbox=checkbox,themeIcons=FALSE,themeDots=FALSE)
	return(treeUI)
}

glycoOnto_tree<-function(glycoOntoStruct){
	tree<-shinyTree::renderTree(glycoOntoStruct)
}


get_geneList<-function(selections,pathwayList,functionList){
	genes<-c(
		 as.character(do.call(c,lapply(selections,function(x) pathwayList[[x]]))),
		 as.character(do.call(c,lapply(selections,function(x) functionList[[x]])))
		 )
	return(sort(genes))
}

#Returns a custom gene selection ui to create a custom panel of glycogenes
tree_geneSelect_ui<-function(inputTree,pathwayList,functionList,geneSelections){
	tree<-inputTree
	req(tree)
	selected<-get_selected(tree)
	#Get genelist:
	geneList<-get_geneList(selected,pathwayList,functionList)
	#Get already-selected genes:
	preSelect_geneList<-geneSelections()[geneSelections() %in% geneList]
	if (length(preSelect_geneList)>0){
		return(selectInput("geneSelection","Select Genes:",choices=geneList,multiple=T,selected=preSelect_geneList))
	} else {
		return(selectInput("geneSelection","Select Genes:",choices=geneList,multiple=T))
	}
}

#Make a static tree which renders the enrichment results:
static_tree_ui<-function(treeName,singleSelect=TRUE,checkbox=TRUE){
	treeUI<-shinyTree(treeName, stripes = TRUE, multiple = singleSelect, animation = FALSE,checkbox=checkbox,themeIcons=FALSE,themeDots=FALSE)
	return(treeUI)
}

#Enrichment Options UI:


build_enrich_tabset<-function(tabName,ctype_code,db){
	#Creates nested tab structure to 
	# parse glycogene enrichment results:
	#Gather subtype levels of tab name:
	subtype_levels<-get_subtype_levels(ctype_code,tabName,db)
	subtype_levels_pretty<-sapply(subtype_levels,function(x) gsub(tabName,'',x))

	tsp<-tabPanel(tabName,
			#Radio select for subtypes:
			radioButtons(paste(tabName,'subtype_choice',sep='_'),label="Select Level",
				     #choices=mapply(function(n,v) n=v,subtype_levels_pretty,subtype_levels)
				     choices=subtype_levels
				     ),
			radioButtons(paste(tabName,'universe_choice',sep='_'),label='Pick Universe',
				     choices=list('Glycogenes'='glycogene','Genomic'='genome')
				     ),
			radioButtons(paste(tabName,'direction_choice',sep='_'),label='Pick Direction',
				     choices=list('Upregulated'='up','Downregulated'='down')
				     )
			)
	return(tsp)
}

render_enrich_tree<-function(subtype,appInput,treeStruct){
	subtype_choices<-names(appInput)[grepl(paste('^',subtype,"\\_",sep=''),names(appInput))]
	if (any(sapply(subtype_choices,function(x) is.null(appInput[[x]])))){
		return(NULL)
	} else {
		#Build Enrichment tree:
		tree_data_with_stats<-function(treeStruct,ctype,db,direction,universe,variable,level){
		rendTree<-renderTree({tree_data_with_stats})
		return(rendTree)
	}
	}
}

enrichment_vars_tabset<-function(ctype,db,build_enrich_tabset,glycoOntoStruct){
	ctype_proj<-get_ctype_project(ctype,TCGA_DB)
	#Gather Unique subtypes:
	subtypes<-get_subtypes(ctype_proj,db)
	subtypes_pretty<-sapply(subtypes,function(x) gsub('\\_',' ',x))
	#Total tab sets:
	tabNames<-c('shortLetterCode',subtypes)
	prettyTabNames<-c('Tumor Tissue',subtypes_pretty)
	#Apply nested tabset structure:
	tbs<-unname(lapply(tabNames,function(x) build_enrich_tabset(x,ctype_proj,db)))
	tsp<-do.call(tabsetPanel,tbs)
	return(tsp)
}

