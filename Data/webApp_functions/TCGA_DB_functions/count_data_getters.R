# Getter functions for the "TCGA_count_data" 
# - Counts required for TSNE and Distribution plots
# Methods required
#
#
#
library(RSQLite)

get_geneSet_counts<-function(ctype,geneSet,db){
	#Construct query:
	# query:
	#  Count table view:
	viewName<-sub('\\-','\\_',ctype)
	viewName<-paste(viewName,'count',sep='_')
	#Fix ctype name
	ctype<-paste("'",ctype,"'",sep='')
	#  Gene List:
	gene_queryList<-paste('(',paste(sapply(geneSet,function(x) paste("'",x,"'",sep='')),collapse=','),')',sep='')
	ctype_part<-paste("ctype=",ctype,sep='')
	#Total:
	query<-paste("SELECT * FROM",viewName,"WHERE",ctype_part,"AND geneName IN",gene_queryList)
	
	#Execute Query
	res<-dbGetQuery(db,query)
	return(res)
}

get_geneSet_countMatrix<-function(ctype,geneSet,db){
	#This function creates expression matrices for use in the TCGA algorithm
	count_data<-get_geneSet_counts(ctype,geneSet,db) %>% 
		pivot_wider(names_from=geneName,values_from=expression,id_cols=c('patient_id'))
	return(count_data)

}

get_subtypeData<-function(ctype_code,subtype,db){
	ctype_code<-paste("'",ctype_code,"'",sep='')
	subtype<-paste("'",subtype,"'",sep='')
	query<-paste("SELECT * FROM TCGA_patient_metadata WHERE",paste("ctype",ctype_code,sep='='),"AND",paste("subtype",subtype,sep='='))
	message(query)
	res<-dbGetQuery(db,query) %>% select(patient_id,subtype_level)
	return(res)
}

get_tumor_normal_data<-function(ctype_code,db){
	ctype_code<-paste("'",ctype_code,"'",sep='')
	query<-paste("SELECT DISTINCT patient_id,shortLetterCode FROM TCGA_patient_metadata WHERE",paste("ctype",ctype_code,sep='='))
	res<-dbGetQuery(db,query)
	#Change labels
	res<-res %>%
		mutate(shortLetterCode=case_when(shortLetterCode=='T'~'Tumor',shortLetterCode=='N'~'Normal'))
	return(res)
}


countMatrix_metadata_merge<-function(countMatrix,ctype,subtype,db){
	#Merges "subtype" label information alongside gene expression data
	#Getting subtype info:
	if (subtype=='shortLetterCode'){
		subtype_dta<-get_tumor_normal_data(ctype,db)
	} else {
		#May want to add error control here:
		subtype_dta<-get_subtypeData(ctype,subtype,db)
	}
	#Merging:
	mergeData<-countMatrix %>% as.data.frame() %>%
		left_join(subtype_dta %>% as.data.frame(),by='patient_id')
	return(mergeData)
}

TCGA_get_tSNE_data<-function(ctype,subtype,geneSet,db){
	#Main method for gathering data in matrix format for use
	# in the t-SNE panel of the webapp.

	#Get count data for cancer type:
	countMatrix<-get_geneSet_countMatrix(ctype,geneSet,db)
	#Merge with metadata:
	countMatrix<-countMatrix_metadata_merge(countMatrix,ctype,subtype,db)
	return(countMatrix)
}


get_gene_counts<-function(geneName,db){
	#Construct query:
	# query:
	#  Gene List:
	geneName=paste("'",geneName,"'",sep='')
	gene_part<-paste("geneName=",geneName,sep="")
	#Total:
	query<-paste("SELECT * FROM TCGA_count_data WHERE",gene_part)
	
	#Execute Query
	res<-dbGetQuery(db,query)
	return(res)
}

TCGA_get_geneCount_data<-function(geneName,db){
	countData<-get_gene_counts(geneName,db)
	#Merge with Tumor-Normal Data:
	tumorNormal_dta<-do.call(rbind,lapply(unique(countData$ctype),function(ct){
				get_tumor_normal_data(ct,db)
			}))
	countData<-countData %>% left_join(tumorNormal_dta,by='patient_id')
	return(countData)
}

#Differential Expression Getter:
get_diffExp_data<-function(ctypes,geneSet,db){
	#Construct query:
	#Fix ctype name
	ctype_queryList<-paste('(',paste(sapply(ctypes,function(x) paste("'",x,"'",sep='')),collapse=','),')',sep='')
	#  Gene List:
	gene_queryList<-paste('(',paste(sapply(geneSet,function(x) paste("'",x,"'",sep='')),collapse=','),')',sep='')
	# Tumor v Normal data only:
	query<-paste("SELECT geneName,ctype,logFC,`adj.P.Val` FROM TCGA_diffExp_data WHERE variable='shortLetterCode' AND level='T' AND ctype IN",ctype_queryList,"AND geneName IN",gene_queryList)
	#Execute Query
	res<-dbGetQuery(db,query)
	return(res)
}

#Enrichment Data getter:
get_enrich_data<-function(ctype,db,direction,universe,variable,level){
	#Construct query:
	# query:
	#  Count table view:
	viewName<-sub('\\-','\\_',ctype)
	viewName<-paste(viewName,'enrich',sep='_')
	#Construct query:
	#Fix ctype name
	ctype<-paste("'",ctype,"'",sep='')
	ctype_part<-paste("ctype=",ctype,sep='')
	#Direction,universe,variable,level component:
	direction<-paste("'",direction,"'",sep='')
	direction_part<-paste("direction=",direction,sep='')
	universe<-paste("'",universe,"'",sep='')
	universe_part<-paste("universe=",universe,sep='')
	variable<-paste("'",level,"'",sep='')
	variable_part<-paste("universe=",level,sep='')
	level<-paste("'",level,"'",sep='')
	level_part<-paste("level=",level,sep='')
	#Total:
	# select based on direction and universe, only extracting 
	middlePart<-paste(c(direction_part,universe_part,level_part),collapse=" AND ")
	query<-paste("SELECT * FROM",viewName,"WHERE",middlePart,"AND",ctype_part)
	#Execute Query
	res<-dbGetQuery(db,query)
	return(res)
}

make_enrichment_tree<-function(tree,dta){
	#Get names at current tree level:
	curNames<-names(tree)
	#Update the names based on the enrichment in dta:
	newNames<-sapply(curNames,function(nm){
		pVal<-dta %>% dplyr::filter(pathways==nm) %>% 
			pull(adjP)
		pVal<-signif(pVal,3)
		newName<-paste(nm,paste('(adj-p-val=',pVal,')',sep=''))
		return(newName)
	})
	names(tree)<-newNames
	for (nm in names(tree)){
		if (class(tree[[nm]])=='list'){
			tree[[nm]]<-make_enrichment_tree(tree[[nm]],dta)
		}
	}
	return(tree)
}

tree_data_with_stats<-function(treeStruct,ctype,db,direction,universe,variable,level){
	#Gather Data:
	enrichDta<-get_enrich_data(ctype,db,direction,universe,variable,level)
	#Update tree data:
	enrichTree<-make_enrichment_tree(treeStruct,enrichDta)
	#This returns a nested list of pathway and function enrichments 
	# with p-values in them, rendered just like PANTHER or DAVID.
	return(enrichTree)
}
