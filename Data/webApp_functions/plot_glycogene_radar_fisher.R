plot_glycogene_radar<-function(dysreg_df,cancer_selection){
  
  number_of_bar<-length(unique(dysreg_df$Pathway))
  angle<-90 - (360 * (seq(1,number_of_bar)-0.5)) /(number_of_bar)
  label_dta<-data.frame(hjust=ifelse(angle< (-90),0,1),angle=ifelse(angle < -90, angle+180, angle))
  # label_dta<-data.frame(hjust=ifelse(angle< (-90),1,0),angle=angle)
  label_dta$label<-unique(dysreg_df[order(dysreg_df$Pathway),]$Pathway)
  label_dta<-cbind(dysreg_df,label_dta)
  
  
  #5. Radar Plotting:
  p<-ggplot(dysreg_df,aes(x=Pathway,y=log10_pValue,colour=Dysregulation,fill=Dysregulation,group=Dysregulation)) + 
    geom_point(size=3.75) + 
    geom_hline(yintercept = (-log10(0.05)),alpha=0.35,size=3.5) +
    geom_polygon(size=1,alpha=0.25) + 
    scale_x_discrete() + 
    theme_light() + 
    scale_color_manual(values=c('red','blue'),name='Dysregulation') + 
    guides(fill=FALSE) +
    # scale_fill_manual(values=c('red','blue')) + 
    labs(title=paste(cancer_selection,'Dysregulated\nGlycosylation Pathways'),x='',y='-log10(pvalue)') + 
    # theme(axis.text.x = element_text(size = 12,hjust = 1),axis.text.y = element_text(size=18),
    #       axis.title.y = element_text(size=25),axis.title.x=element_text(size=25),legend.title = element_text(size=22.5),
    #       legend.text=element_text(size=18.5), legend.key.size = unit(0.65,units = 'cm'),title=element_text(size=25)) +
    theme(axis.text.x = element_blank(),axis.text.y = element_text(size=18),
          axis.title.y = element_text(size=25),axis.title.x=element_text(size=25),legend.title = element_text(size=22.5),
          legend.text=element_text(size=18.5), legend.key.size = unit(0.65,units = 'cm'),title=element_text(size=25)) +
    geom_text(data=label_dta,aes(x=label,y=rep(1.25*max(dysreg_df$log10_pValue),dim(label_dta)[1]),label=label),
              angle=label_dta$angle,color='black',size=4.5) +
    coord_polar() 
  return(p)
}