plot_gene_dists<-function(dta,opts){
  geneName<-unique(dta$geneName)
  dta$shortLetterCode<-factor(dta$shortLetterCode,levels=c('Normal','Tumor'),ordered=T)
  if (opts=='Binned Bar'){
  ggplot(dta,aes(x=ctype,y=expression,fill=shortLetterCode)) + stat_bin_2d(aes(fill=shortLetterCode),bins=300) + 
    theme_bw() + 
    scale_fill_manual(values=c('Tumor'='red','Normal'='blue'),name='Tissue Type') +
    #scale_fill_manual(name='Tissue Type') +
    labs(x='Cancer Types',y='VST',title=paste(geneName,'Binned Expression Across Cancers')) + 
    coord_flip() +
    theme(axis.text.x = element_text(size = 18),axis.text.y = element_text(size=20),
          axis.title.y = element_text(size=25),axis.title.x=element_text(size=25),legend.title = element_text(size=22.5),
          legend.text=element_text(size=18.5), legend.key.size = unit(0.65,units = 'cm'),title=element_text(size=25))
  } else if (opts=='Scatter Plot'){
    ggplot(dta,aes(x=ctype,y=expression,colour=shortLetterCode)) + geom_jitter(alpha=0.25)+ 
      theme_bw() + 
      scale_colour_manual(values=c('Tumor'='red','Normal'='blue'),name='Tissue Type') +
      labs(x='Cancer Types',y='VST',title=paste(geneName,'Expression Across Cancers')) + 
      guides(colour = guide_legend(override.aes = list(size=10,alpha=1))) +
      coord_flip() +
      theme(axis.text.x = element_text(size = 18),axis.text.y = element_text(size=20),
            axis.title.y = element_text(size=25),axis.title.x=element_text(size=25),legend.title = element_text(size=22.5),
            legend.text=element_text(size=18.5), legend.key.size = unit(0.65,units = 'cm'),title=element_text(size=25))
  }
}
