parse_dysreg_path<-function(PathList,all_DE_data_pst_list_sub,cancer_selection){
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
  #3. Compute the percentage of each glycogene pathway represented:
  cancer_upreg_glycan_tally<-matrix(nrow=1,ncol=length(names(PathList)),
                                    dimnames = list(x='',y=names(PathList)))
  cancer_downreg_glycan_tally<-matrix(nrow=1,ncol=length(names(PathList)),
                                      dimnames = list(x='',y=names(PathList)))
  for (p in names(PathList)){
    up_pct<-100*(sum(upreg_geneList %in% PathList[[p]])/length(PathList[[p]]))
    down_pct<-100*(sum(downreg_geneList %in% PathList[[p]])/length(PathList[[p]]))
    # message(paste(c,p,'up percent: ',up_pct))
    # message(paste(c,p,'down percent: ',down_pct))
    cancer_upreg_glycan_tally[1,p]<-up_pct
    cancer_downreg_glycan_tally[1,p]<-down_pct
  }
  
  #4. Manipulate Data
  
  cancer_upreg_glycan_tally<-as.data.frame(cancer_upreg_glycan_tally)
  cancer_downreg_glycan_tally<-as.data.frame(cancer_downreg_glycan_tally)
  
  cancer_upreg_glycan_tally$dysreg<-'up'
  cancer_downreg_glycan_tally$dysreg<-'down'
  
  cancer_dysreg_df<-rbind(cancer_upreg_glycan_tally,cancer_downreg_glycan_tally)
  cancer_dysreg_df_melt<-melt(cancer_dysreg_df)
  cancer_dysreg_df_melt$dysreg<-factor(cancer_dysreg_df_melt$dysreg,levels = c('up','down'))
  # cancer_dysreg_df_melt$variable[cancer_dysreg_df_melt$variable=='BloodGroup_pathway']='Blood_Group_pathway'
  cancer_dysreg_df_melt$variable<-sapply(cancer_dysreg_df_melt$variable,
                                         function(x) gsub(pattern = '\\_',replacement = '\n',x))
  return(cancer_dysreg_df_melt)
}