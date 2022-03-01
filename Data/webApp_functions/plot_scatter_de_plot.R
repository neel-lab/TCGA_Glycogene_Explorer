plot_scatter_DE_plot<-function(dta,SNFG_functionColors){
  #source('./Data/SNFG_mapping_function.R')
  #Genes that have only one disease associated with it:
  dta$sig<-factor(sapply(dta$`adj.P.Val`,function(x) ifelse(x<=0.05,TRUE,FALSE)),levels = c('FALSE','TRUE'))

  ggene_types<-SNFG_functionColors %>% dplyr::filter(ggene_function %in% unique(dta$gene_function))
  ggene_types_vec<-setNames(ggene_types$ggene_function,ggene_types$color)
  message(ggene_types_vec)

  textColor<-unique(dta$gene_color)
  shapevec<-c('TRUE'=16,'FALSE'=17)
  p<-ggplot() +  
    geom_boxplot(data=dta,aes(x=geneName,y=logFC,fill=gene_function),alpha=0.55) +
    #Significant Data:
    geom_point(data=dta,aes(x=geneName,y=logFC,
                            shape=sig),color='black',alpha=0.9,size=6.75) +
    geom_point(data=dta,aes(x=geneName,y=logFC,colour=ctype,
                                            shape=sig),alpha=0.9,size=5)+
    # geom_point(data=dta[dta$`adj.P.Val`>0.05,],aes(x=gene,y=log2FC,colour=disease,size=-log10(padj)),shape=17,alpha=0.9) +
    # geom_point(data=dta,aes(x=gene,y=log2FC,colour=disease),size=4.5,alpha=0.9)+
    # geom_jitter(data=dta[dta$geneName %in% multiples,],aes(x=gene,y=log2FC,colour=disease),width = 0.2,size=4.5,alpha=0.9) + 
    geom_hline(yintercept = 0,color='black',size=1.25) +
    scale_fill_manual(name='Glycogene Type',values=ggene_types_vec) +
    scale_color_discrete(name='Cancer Types') + 
    labs(title=paste('Tumor vs Normal Glycogene Differential Expression')) +
    scale_shape_discrete(guide=FALSE) +
    scale_shape_manual(values = shapevec,name='Significant?') +
    guides(colour=guide_legend(override.aes = list(size=5,shape=16)),shape=guide_legend(override.aes = list(size=5))) +
    theme_light() + 
    theme(axis.text.x = element_text(angle = 45,size = 21.5,hjust=1),axis.text.y = element_text(size=20),
          axis.title.y = element_text(size=25),axis.title.x=element_text(size=25),legend.title = element_text(size=22.5),
          legend.text=element_text(size=18.5), legend.key.size = unit(0.65,units = 'cm'),title=element_text(size=25))

  return(p)
}
