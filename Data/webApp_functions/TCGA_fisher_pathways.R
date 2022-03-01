TCGA_fisher_pathways<-function(all_DE_data_pst_list_sub,ctype,pathGenes,geneSpace,direction='up'){
  #First, parse DE data
  dta<-all_DE_data_pst_list_sub[[ctype]]
  dta<-dta[row.names(dta) %in% geneSpace,]
  
  fisher_stats<-sapply(names(pathGenes),function(x){
    if (direction=='up'){
      pathDysreg<-dim(dta[(row.names(dta) %in% pathGenes[[x]]) & (as.numeric(as.character(dta$log2FC))>0) & (as.numeric(as.character(dta$padj))<=0.05),])[1]
      totalDysreg<-dim(dta[(as.numeric(as.character(dta$log2FC))>0) & (as.numeric(as.character(dta$padj))<=0.05),])[1]
    } else if (direction=='down'){
      pathDysreg<-dim(dta[(row.names(dta) %in% pathGenes[[x]]) & (as.numeric(as.character(dta$log2FC))<0) & (as.numeric(as.character(dta$padj))<=0.05),])[1]
      totalDysreg<-dim(dta[(as.numeric(as.character(dta$log2FC))<0) & (as.numeric(as.character(dta$padj))<=0.05),])[1]
    }
    
    cTable<-t(matrix(c(pathDysreg,(length(pathGenes[[x]])-pathDysreg),(totalDysreg-pathDysreg),length(geneSpace)),ncol=2))
    return(fisher.test(cTable,alternative = 'greater')$p.value)
  })
  
  return(fisher_stats)
}