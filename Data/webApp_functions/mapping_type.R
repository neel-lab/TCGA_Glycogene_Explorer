mapping_type<-function(gene_list,Glycogene_types_list){
  result<-apply(as.matrix(sapply(gene_list,function(y) sapply(names(Glycogene_types_list),function(x) y %in% Glycogene_types_list[[x]]))),2,function(x) names(x[x==TRUE])[1] )
  result[is.na(result)==TRUE]='Unknown'
  return(result)
}