plot_scatter_DE_plot<-function(dta,SNFG_color_mapping_vector){
  source('./Data/SNFG_mapping_function.R')
  #Genes that have only one disease associated with it:
  singles_inds<-which(sapply(unique(dta$gene), function(x) length(unique(dta[dta$gene==x,'disease'])))==1)
  multiples_inds<-seq(1,length(unique(dta$gene)))[-singles_inds]
  singles<-unique(dta$gene)[singles_inds]
  multiples<-unique(dta$gene)[multiples_inds]
  dta$sig<-factor(sapply(dta$padj,function(x) ifelse(x<=0.05,'padj <= 0.05','padj > 0.05')),
                  levels = c('padj <= 0.05','padj > 0.05'))
  shapevec<-c('padj <= 0.05'=16,'padj > 0.05'=17)
  agl<-ifelse(dim(dta)[1]>10,90,45)
  p<-ggplot() +  
    geom_boxplot(data=dta,aes(x=labels,y=log2FC,fill=gene_type),alpha=0.55) +
    #Significant Data:
    geom_point(data=dta,aes(x=labels,y=log2FC,
                            shape=sig),color='black',alpha=0.9,size=6.75) +
    geom_point(data=dta,aes(x=labels,y=log2FC,colour=disease,
                                            shape=sig),alpha=0.9,size=5)+
    # geom_point(data=dta[dta$padj>0.05,],aes(x=gene,y=log2FC,colour=disease,size=-log10(padj)),shape=17,alpha=0.9) +
    # geom_point(data=dta,aes(x=gene,y=log2FC,colour=disease),size=4.5,alpha=0.9)+
    # geom_jitter(data=dta[dta$gene %in% multiples,],aes(x=gene,y=log2FC,colour=disease),width = 0.2,size=4.5,alpha=0.9) + 
    geom_hline(yintercept = 0,color='black',size=1.25) +
    scale_fill_manual(name='Glycogene Type',values=SNFG_color_mapping_vector) +
    scale_color_discrete(name='Cancer Types') + 
    labs(title=paste('Tumor vs Normal Glycogene Differential Expression'),x='Gene [Function]') +
    # scale_size_continuous(name='Significance (-log10(padj))',range = c(3,10)) + 
    scale_shape_discrete(guide=FALSE) +
    scale_shape_manual(values = shapevec,name='Statistical Significance') +
    guides(colour=guide_legend(override.aes = list(size=5,shape=16)),shape=guide_legend(override.aes = list(size=5))) +
    labs(y='log2(Tumor/Normal)') +
    theme_light() + 
    theme(axis.text.x = element_text(angle = agl,size = 15,hjust=1,color = unique(dta[,c('gene','gene_color')])$gene_color),axis.text.y = element_text(size=20),
          axis.title.y = element_text(size=25),axis.title.x=element_text(size=25),legend.title = element_text(size=22.5),
          legend.text=element_text(size=18.5), legend.key.size = unit(0.65,units = 'cm'),
          title=element_text(size=25))
  
  # labelPos<-data.frame(xstart=c(),xend=c(),ypos=c(),lincol=c())
  # labelSize<-ifelse(length(unique(dta$gene))>5,5,8)
  # stagger=TRUE
  # for (f in unique(dta$Function)){
  #   #Find indecies of the function:
  #   xinds<-which(dta$Function==f)
  #   #Use gene Indeces:
  #   xstart<-dta[min(xinds),'gene']
  #   xend<-dta[max(xinds),'gene']
  #   # xstart<-min(xinds)
  #   # xend<-max(xinds)
  #   #Provide the gene names for the given function
  #   x=dta[which(dta$Function==f),'gene']
  #   
  #   ypos<-max(dta$log2FC)*1.45
  #   lincol<-unique(dta[which(dta$Function==f),'gene_color'])[1]
  #   if (xstart!=xend){
  #     textloc<-mean(c(xstart,xend))
  #   } else {
  #     textloc<-xstart
  #   }
  #   labelPos<-rbind(labelPos,cbind(xstart,xend,ypos,lincol))
  #   # labelPos$xstart<-c(labelPos$xstart,xstart)
  #   # labelPos$xend<-c(labelPos$xend,xend)
  #   # labelPos$ypos<-c(labelPos$ypos,ypos)
  #   # labelPos$lincol<-c(labelPos$lincol,lincol)
  #   if (stagger==TRUE){
  #     staggerMult=1.075
  #     stagger=FALSE
  #   } else {
  #     staggerMult=0.95
  #     stagger=TRUE
  #   }
  #   p<-p + 
  #     # annotate("segment",x=x,y=rep(ypos,length(x)),color=lincol,size=5,family="")
  #     #annotate("segment",x=xstart,xend=xend,y=ypos,yend=ypos,color=lincol,size=5,family="") +
  #     # geom_segment(aes(x=xstart,xend=xend,y=(ypos*staggerMult)-0.5,yend=(ypos*staggerMult)-0.5),color=lincol,size=5) +
  #         #annotate("text",x=textloc,y=ypos*1.25,label=f,color=lincol,size=7,family="")
  #         annotate("text",x=textloc,y=ypos*staggerMult,label=f,color=lincol,size=labelSize,family="")
  # }
  # # p<-p + geom_segment(aes(x=labelPos$xstart,xend=labelPos$xend,y=labelPos$ypos,
  # #                 yend=labelPos$ypos),color=labelPos$lincol,size=2.5)
  
  return(p)
  
}