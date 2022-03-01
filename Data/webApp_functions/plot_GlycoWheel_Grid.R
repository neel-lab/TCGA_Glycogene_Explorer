plot_glycoWheel<-function(all_DE_data_pst_list_sub,geneClass,Glycogene_types_list){
  #Parses the gene classes and creates glycoWheels
  library(stringr)
  mapVec<-sapply(names(Glycogene_types_list),function(x){
    if (is.na(str_extract(string=x,pattern = 'transfer'))){
      return(paste(x,'Wheel',sep = ''))
    } else{
      return(gsub(pattern = '(transfer)?ase',replacement = 'Wheel',x=x))
    }
  })
  
  wheelData<-lapply(names(all_DE_data_pst_list_sub),
                        function(x) all_DE_data_pst_list_sub[[x]][row.names(all_DE_data_pst_list_sub[[x]]) %in% Glycogene_types_list[[geneClass]],])
  names(wheelData)<-names(all_DE_data_pst_list_sub)
  #Create Wheels:
  wheelList<-list()
  for (n in names(wheelData)){
    dta<-wheelData[[n]]
    dta$log2FC<-as.numeric(as.character(dta$log2FC))
    dta$padj<-as.numeric(as.character(dta$padj))
    dta<-dta[order(row.names(dta)),]
    cls<-ifelse(dta$log2FC<0,'Down','Up')
    dta$dysreg<-factor(cls,levels = c('Up','Down','Not Significant'))
    mapVec<-list('Up'='red','Down'='blue','Not Significant'='grey')
    cols<-sapply(dta$dysreg,function(x) mapVec[[x]])
    
    number_of_bar<-nrow(dta)
    angle<-90 - (360 * (seq(1,dim(dta)[1])-0.5)) /(number_of_bar)
    label_dta<-data.frame(hjust=ifelse(angle< (-90),0,1),angle=ifelse(angle < -90, angle+180, angle))
    # label_dta<-data.frame(hjust=ifelse(angle< (-90),1,0),angle=angle)
    row.names(label_dta)<-row.names(dta)
    label_dta<-cbind(dta,label_dta)
    
    p<-ggplot(dta) + geom_bar(aes(x=row.names(dta),y=log2FC,fill=dysreg),stat='identity') +
      labs(title=n,y='log2FC',x='') +  
      scale_y_continuous(limits = c(-7.5,7.5),position = 'left') + 
      theme_light() +
      guides(color=FALSE) + 
      theme(axis.text.x = element_blank(),legend.position = 'none',axis.text.y=element_text(size=15),
            title=element_text(size=20))+
      # theme(axis.title=element_text(size=10),axis.text=element_text(size = 7,hjust = 1),
      #       title=element_text(size=12),legend.position = 'none') +
      coord_polar(start = 0,direction = 1) +
      geom_text(data=label_dta,aes(x=row.names(dta),y=rep(7.5,dim(dta)[1]),label=row.names(label_dta)),angle=label_dta$angle)
    
    wheelList[[n]]<-p
  }
  #Order higher expressing cancers first in grid
  DE_lfcsum<-sapply(names(wheelData),function(x) sum(as.numeric(as.character(wheelData[[x]][,'log2FC']))))
  order_to_plot<-names(DE_lfcsum[order(-DE_lfcsum)])
  wheelList_ordered<-lapply(order_to_plot,function(x) wheelList[[x]])
  # names(wheelList_ordered)<-order_to_plot
  
  #Save plots in pairs :
  outputPlotList<-list()
  layout_matrix<-matrix(nrow=1,ncol = 2,seq(1,2))
  ind=1
  for (i in seq(1,length(wheelList_ordered),by = 2)){
    wheelPartList<-list(wheelList_ordered[[i]],wheelList_ordered[[i+1]])
    grobs<-sapply(wheelPartList,function(x) ggplotGrob(x))
    outputPlotList[[ind]]<-arrangeGrob(wheelPartList,grobs=grobs,layout_matrix=layout_matrix)
    ind=ind+1
  }
  return(outputPlotList)
  
  # #Break Down the plots into multiple grobs
  # #TitleGrob
  # grobs1<-sapply(wheelList_ordered[1:10],function(x) ggplotGrob(x))
  # layout_matrix<-matrix(nrow=2,ncol=10,seq(1,10))
  # wheelGrid1<-arrangeGrob(wheelList_ordered[1:10],grobs=grobs1,layout_matrix=layout_matrix)
  # #SecondGrob
  # grobs2<-sapply(wheelList_ordered[11:20],function(x) ggplotGrob(x))
  # wheelGrid2<-arrangeGrob(wheelList_ordered[11:20],grobs=grobs2,layout_matrix = layout_matrix)
  # #FinalGrob
  # grobs3<-sapply(wheelList_ordered[21:22],function(x) ggplotGrob(x))
  # wheelGrid3<-arrangeGrob(wheelList_ordered[21:22],grobs=grobs3,layout_matrix=matrix(nrow=1,ncol = 2,seq(1,2)))
  # 
  # return(list('TitleBlock'=wheelGrid1,'MiddleBlock'=wheelGrid2,'LastBlock'=wheelGrid3,'PlotHeader'=mapVec[geneClass]))
  
}