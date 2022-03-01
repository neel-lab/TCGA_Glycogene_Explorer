library(RSQLite)

get_enrichment_dta<-function(ctype,subtype,universe,direction,db){
	#Construct Query
	ctype<-paste("'",ctype,"'",sep='')
	subtype<-paste("'",subtype,"'",sep='')
	universe<-paste("'",universe,"'",sep='')
	direction<-paste("'",direction,"'",sep='')
	query<-paste('SELECT * FROM TCGA_glycoPathway_enrichments WHERE',paste('ctype',ctype,sep='='),'AND',
		     paste('level',subtype,sep='='),'AND',
		     paste('universe',universe,sep='='),'AND',
		     paste('direction',direction,sep='='),
		     sep=' ')
	res<-dbGetQuery(db,query)
	return(res)
}

relabeler<-function(pthName,stat){
	return(paste(pthName,paste('( adj-p-val = ',as.character(formatC(stat,format='e',digits=3)),' )',sep='')))
}

build_enrich_tree<-function(glycoOnto_hierarchy,enrichDta){
	#Recursively rename 
	names(glycoOnto_hierarchy)<-sapply(names(glycoOnto_hierarchy),function(n){
		if (n %in% c('Glycogene','Pathway','Compartment')){
			return(n)
		}
		stat<-enrichDta %>% dplyr::filter(pathways==n) %>% pull(adjP)
		newName<-relabeler(n,stat)
		return(newName)
	})
	for (p in names(glycoOnto_hierarchy)){
		if (is.list(glycoOnto_hierarchy[[p]])){
			glycoOnto_hierarchy[[p]]<-build_enrich_tree(glycoOnto_hierarchy[[p]],enrichDta)

		} else {
			next
		}
	}
	return(glycoOnto_hierarchy)
}
