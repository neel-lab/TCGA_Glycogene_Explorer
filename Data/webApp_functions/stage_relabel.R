stage_relabel <-function(data){
# Simplify the stage labeling
# Takes TCGA dataframes with tumor stage info and relabels with simpler labeling
	 stage_subset_labels<-c('a','b','c')
	 stagelabels<-as.character(data$tumor_stage)
	 stagelabels_test<-stagelabels[stagelabels!='not reported']
	 relabel = FALSE
	 for (stage in as.character(stagelabels_test)){
	 	label_elts<-unlist(strsplit(as.character(stage),' '))
	 	#Check to see if stage label has "stage" then "ii*"
	 	#If subset labels are detected, mark to relabel
	 	if (length(label_elts) > 1){
	 		stagenum<-label_elts[2]
	 		logical<-sapply(stage_subset_labels,function(x) grepl(x,stagenum))
	 		if (sum(logical) >=1){
	 			relabel = TRUE
	 			#message(paste(disease,'marked for relabeling'))
	 			break }
	 		
	 	}
	 }
	 if (relabel==TRUE){
	 	#Create new stage names
	 	splitting_stages<-sapply(stagelabels,function(x) unlist(strsplit(x,' '))[2])
	 	condensing_stages<-c()
	 	for (stage in splitting_stages){
	 		subset_val<-sapply(stage_subset_labels,function(x) grepl(x,stage))
	 		subset_val<-sum(subset_val)
	 		if ((nchar(stage)>=2) & (subset_val>=1)){
	 			new_stage<-substr(stage,start=1,stop=(nchar(stage)-1))
	 			condensing_stages<-c(condensing_stages,new_stage)
	 		}
	 		else if ((nchar(stage)<2) & (subset_val>=1)){
	 			new_stage<-substr(stage,start=1,stop=1)
	 			condensing_stages<-c(condensing_stages,new_stage)
	 		}
	 		else {
	 			condensing_stages<-c(condensing_stages,stage)
	 		}
	 	}
	 	condensing_stages[condensing_stages=='reported']='not reported'
	 	#For the data with stages corresponding to normal tissue, relabel the stage as "Normal"
	 	#Second row of tissue_and_stage will be modified with appropriate labels
	 	tissue_and_stage<-cbind(as.character(data$shortLetterCode),as.character(condensing_stages))
	 	for (row in 1:dim(tissue_and_stage)[1]){
	 	  if (tissue_and_stage[row,1]=='NT'){
	 	    tissue_and_stage[row,2]='Normal'
	 	  }
	 	}
	 	data$tumor_stage<-as.factor(tissue_and_stage[,2])
	 }
	#DONE RELABELING
	return(data)
}
