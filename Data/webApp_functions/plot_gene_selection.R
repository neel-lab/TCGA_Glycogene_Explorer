plot_gene_selection<-function(data_parse){
  p<-ggplot(data_parse) + geom_jitter(aes(x=variable,y=value,colour=definition),width=0.25,alpha=0.55,size=3) + 
    theme_light()
  return(p)
  
}