plot_glycogene_radar<-function(dysreg_df,cancer_selection){
  #5. Radar Plotting:
  p<-ggplot(dysreg_df,aes(x=variable,y=value,colour=dysreg,fill=dysreg,group=dysreg)) + 
    geom_point(size=3.75) + 
    geom_polygon(size=1,alpha=0.25) + 
    scale_x_discrete() + 
    theme_light() + 
    ylim(0,100) + 
    scale_color_manual(values=c('red','blue'),name='Dysregulation') + 
    guides(fill=FALSE) +
    # scale_fill_manual(values=c('red','blue')) + 
    labs(title=paste(cancer_selection,'Dysregulated\nGlycosylation Pathways'),x='',y='') + 
    theme(axis.text.x = element_text(size = 17,hjust = 1),axis.text.y = element_text(size=18),
          axis.title.y = element_text(size=25),axis.title.x=element_text(size=25),legend.title = element_text(size=22.5),
          legend.text=element_text(size=18.5), legend.key.size = unit(0.65,units = 'cm'),title=element_text(size=25)) +
    coord_polar() 
  return(p)
}