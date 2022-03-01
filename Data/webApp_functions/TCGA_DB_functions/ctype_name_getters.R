library(RSQLite)

get_ctype_names<-function(ctype_code,db){
	#Construct query:
	ctype_code<-paste("'",ctype_code,"'",sep='')
	query<-paste("SELECT DISTINCT disease_type FROM TCGA_Diseases WHERE project_id LIKE 'TCGA-%'")
	#Execute Query
	res<-dbGetQuery(db,query) %>% pull(disease_type)
	return(res)
}


get_ctype_name<-function(ctype_code,db){
	#Construct query:
	ctype_code<-paste("'",ctype_code,"'",sep='')
	query<-paste("SELECT disease_type FROM TCGA_Diseases WHERE",paste("project_id",ctype_code,sep="="))
	#Execute Query
	res<-dbGetQuery(db,query) %>% pull(disease_type)
	return(res)
}

get_ctype_project<-function(ctype_name,db){
	#Construct query:
	ctype_name<-paste("'",ctype_name,"'",sep='')
	query<-paste("SELECT project_id FROM TCGA_Diseases WHERE",paste("disease_type",ctype_name,sep="="))
	#Execute Query
	res<-dbGetQuery(db,query) %>% pull(project_id)
	return(res)
}

get_subtypes<-function(ctype_code,db){
	ctype_code<-paste("'",ctype_code,"'",sep='')
	query<-paste("SELECT DISTINCT subtype FROM TCGA_patient_metadata WHERE",paste("ctype",ctype_code,sep='='))
	res<-dbGetQuery(db,query) %>% pull(subtype)
	return(res)
}


get_subtype_levels<-function(ctype_code,subtype,db){
	ctype_code<-paste("'",ctype_code,"'",sep='')
	subtype<-paste("'",subtype,"'",sep='')
	query<-paste("SELECT DISTINCT level FROM TCGA_glycoPathway_enrichments WHERE",paste("ctype",ctype_code,sep='='),"AND",paste("variable",subtype,sep="="))
	res<-dbGetQuery(db,query) %>% pull(level)
	return(res)
}

get_subtype_data_name<-function(ctype_code,subtype,db){
	ctype_code<-paste("'",ctype_code,"'",sep='')
	subtype<-paste("'",subtype,"'",sep='')
	query<-paste("SELECT DISTINCT subtype FROM TCGA_patient_metadata WHERE",paste("ctype",ctype_code,sep='='),"AND",paste("subtype",subtype,sep='='))
	res<-dbGetQuery(db,query) %>% pull(subtype)
}
