parse_dysreg_path_fisher<-function(PathList,all_DE_data_pst_list_sub,cancer_selection){
  #1. Parse relevant DE data:
  dta=all_DE_data_pst_list_sub[[cancer_selection]]
  dta$log2FC<-as.numeric(as.character(dta$log2FC))
  dta$padj<-as.numeric(as.character(dta$padj))
  dta$padj[is.na(dta$padj)==TRUE]=1
  #2. Find all significantly up and downregulated genes, store in a list:
  upreg_geneList<-c()
  downreg_geneList<-c()
  for (g in row.names(dta)){
    if (dta[g,'log2FC']<0 && dta[g,'padj']<=0.05){
      downreg_geneList<-c(downreg_geneList,g)
    } else if (dta[g,'log2FC']>0 && dta[g,'padj']<=0.05){
      upreg_geneList<-c(upreg_geneList,g)
    }
  }
  
  
  #3. Do the enrichment test:
  
  #Store total gene space in vector:
  geneSpace<-unique(unlist(PathList))
  
  dysreg_mat<-cbind( c(rep('up',length(PathList)),rep('down',length(PathList))),
                    c(TCGA_fisher_pathways(all_DE_data_pst_list_sub,cancer_selection,PathList,geneSpace,'up'),
                      TCGA_fisher_pathways(all_DE_data_pst_list_sub,cancer_selection,PathList,geneSpace,'down')) )
  dysreg_mat<-cbind(row.names(dysreg_mat),dysreg_mat)
  row.names(dysreg_mat)<-seq(1,dim(dysreg_mat)[1])
  dysreg_mat<-as.data.frame(dysreg_mat)
  colnames(dysreg_mat)<-c('Pathway','Dysregulation','log10_pValue')
  dysreg_mat$log10_pValue<-(-log10(as.numeric(as.character(dysreg_mat$log10_pValue))))
  return(dysreg_mat)
}