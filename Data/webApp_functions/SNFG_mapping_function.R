get_ggene_functions<-function(ggene,functionList){
	return(names(functionList)[sapply(names(functionList),function(x) ggene %in% functionList[[x]])])
}

top_geneFunction<-function(ggene,functionList){
	return(get_ggene_functions(ggene,functionList)[1])
}

SNFG_mapping_function<-function(ggene,functionList,SNFG_functionColors){
  ggene_functions<-get_ggene_functions(ggene,functionList) 
  colorVec<-SNFG_functionColors %>% 
	 dplyr::filter(ggene_function %in% ggene_functions) %>% 
	 pull(color) 
  #If has multiple, just return the first one:
  return(unique(colorVec)[1])
}

SNFG_map_main<-function(df,functionList,SNFG_functionColors){
	df<- df %>% 
		mutate(gene_color=sapply(geneName,function(x) SNFG_mapping_function(x,functionList,SNFG_functionColors)),
		       gene_function=sapply(geneName,function(x) top_geneFunction(x,functionList))
		       )
	return(df)
}
