library(RSQLite)




get_diff_exp_dta<-function(ctype_list,geneList,db){
	#For a list of cancer types and glycogenes,
	# get differential expression data.
	viewName<-paste(ctype,'DEA',sep='_')
	viewName<-sub('\\-','\\_',viewName)
	#  Ctype List
	ctype_queryList<-paste('(',paste(sapply(ctype_list,function(x) paste("'",x,"'",sep='')),collapse=' '),')',sep='')
	#  Gene List:
	gene_queryList<-paste('(',paste(sapply(geneList,function(x) paste("'",x,"'",sep='')),collapse=' '),')',sep='')
	#Total:
	query<-paste("SELECT * FROM",viewName,"WHERE ctype IN",ctype_queryList,"AND geneName IN",gene_queryList)
	#Execute Query
	res<-dbGetQuery(db,query)
	return(res)
}
