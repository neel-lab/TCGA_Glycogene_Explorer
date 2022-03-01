library(RSQLite)
library(dplyr)


get_ctype_subTypes<-function(ctype_code,db){
	#Construct query:
	ctype_code<-paste("'",ctype_code,"'",sep='')
	query<-paste("SELECT disease_type FROM TCGA_Diseases WHERE",paste("project_id",ctype_code,sep="="))
	#Execute Query
	res<-dbGetQuery(db,query) %>% pull(disease_type)
	return(res)
}
