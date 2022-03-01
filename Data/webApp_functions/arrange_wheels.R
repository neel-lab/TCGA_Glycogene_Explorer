arrange_wheels<-function(plotList,geneClass,Glycogene_types_list){
  mapVec<-sapply(names(Glycogene_types_list),function(x){
    if (is.na(str_extract(string=x,pattern = 'transfer'))){
      return(paste(x,'Wheel',sep = ''))
    } else{
      return(gsub(pattern = '(transfer)?ase',replacement = 'Wheel',x=x))
    }
  })
  output_grids<-lapply(1:length(plotList),function(x){
    if (x==1){
      tlt<-mapVec[geneClass]
      renderPlot({ grid.arrange(plotList[[x]],
                   top=textGrob(tlt,gp=gpar(fontsize=30)))
      })
    } else {
      renderPlot({
      grid.arrange(plotList[[x]])
      })
    }
    
    
  })
  do.call(tagList,output_grids)
  return(output_grids)
}