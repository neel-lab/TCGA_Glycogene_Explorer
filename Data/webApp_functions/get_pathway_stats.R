get_pathway_stats<-function(PathList,all_DE_data_pst_list_sub,cancer_selection,pathway){
  library(DT)
  #1. Parse relevant DE data:
  dta=all_DE_data_pst_list_sub[[cancer_selection]]
  dta$log2FC<-as.numeric(as.character(dta$log2FC))
  dta$padj<-as.numeric(as.character(dta$padj))
  
  #2. Find all significantly up and downregulated genes, store in a list:
  dta_parse<-dta[row.names(dta) %in% PathList[[pathway]],]
  # dta_parse<-dta_parse[dta_parse$padj<=0.05,]
  dta_parse<-dta_parse[order(-dta_parse$log2FC),]
  dta_parse<-cbind(row.names(dta_parse),dta_parse)
  colnames(dta_parse)[1]<-'Glycogenes'

  return(dta_parse)
}