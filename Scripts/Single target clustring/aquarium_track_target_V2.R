###Analysing single target from aquarium trials
source("Scripts/Single target clustring/clustring_functions.R")


###Apply EMD functions
##Read data frames
#1. FR data
survey_data <- read.csv("Data/Single target analysis/Aquarium/FR data to EMD.csv", stringsAsFactors=FALSE)
#2. target clustring metadata
target_metadata <- read.csv("Data/Single target analysis/Aquarium/FR_cluster_index.csv", stringsAsFactors=FALSE)
#target_metadata$Target_index=as.character(target_metadata$Target_index)

##Subset 'survey_data' to include only tarets with mutiple records 
survey_data_sub=survey_data[survey_data$Target_index %in% unique(target_metadata$Target_index),]

##Take first 100 rows
survey_data_sub=survey_data_sub[1:100,]
##Apply functions
#With absolut values
Aq_test2_100rows=emd_distance_matrix(survey_data_sub,percent_in=0.9)
colnames(Aq_test2_100rows)=c("target1","target2","emd")
Aq_test2_100rows$target1=as.numeric(Aq_test2_100rows$target1)
Aq_test2_100rows$target2=as.numeric(Aq_test2_100rows$target2)
#With relative freq trend
Aq_test2_100rows_rel=emd_distance_matrix(survey_data_sub,percent_in=0.9,relative=T)
colnames(Aq_test2_100rows_rel)=c("target1","target2","emd")
Aq_test2_100rows_rel$target1=as.numeric(Aq_test2_100rows_rel$target1)
Aq_test2_100rows_rel$target2=as.numeric(Aq_test2_100rows_rel$target2)


##Merge with "Region_ID"
plot_pairs_from_output=function(emd_data,n_plot=1,relative_in=F){
  target_metadata=target_metadata[target_metadata$Target_index %in% emd_data$target1,]
  region_unique=unique(target_metadata$Region_ID)
  for (r in n_plot){
    print(r)
    region_c=region_unique[r]
    #get the target_id for this region_c
    targets_ids=target_metadata[target_metadata$Region_ID==region_c,"Target_index"]
    if(length(targets_ids)>1){
      #subset emd data
      region_emd=emd_data[emd_data$target1 %in% targets_ids & emd_data$target2 %in% targets_ids,]
      
      region_emd_targets=unique(c(region_emd$target1,region_emd$target2))
      ##Plot pairs
      print(plot_pair(survey_data_sub,region_emd_targets,emd_data,relative_in))}else{}
  }}

n_plot_in=c(18:20)
plot_pairs_from_output(Aq_test2_100rows,n_plot=n_plot_in)
plot_pairs_from_output(Aq_test2_100rows_rel,n_plot=n_plot_in,relative_in=T)

##matrix distance
dist_matrix=function(distance_df){
  #create new matrix
  vals<-sort(unique(c(as.character(distance_df$target1), as.character(distance_df$target2))))
  nm <- matrix(NA, nrow=length(vals), ncol=length(vals), dimnames=list(vals, vals))
  diag(nm) <- 0   
  
  #fill
  distance_df$target1=as.character(distance_df$target1)
  distance_df$target2=as.character(distance_df$target2)
  nm[as.matrix(distance_df[, 1:2])] <- distance_df[,3]
  #fill reversed
  nm[as.matrix(distance_df[, 2:1])] <- distance_df[,3]
  nm
  return(nm)
}
##Apply function
emd_matrix=dist_matrix(Aq_test2_100rows)


##cluster with NDS
# Classical MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name

d <- emd_matrix # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
fit # view results
##ggplot
fit_df=data.frame(fit$points)
fit_df$Target_index=row.names(fit[["points"]])
fit_df=merge(fit_df,survey_data_sub[,c("Target_index","Depth")],by="Target_index")
target_metadata1=target_metadata
target_metadata1$Target_index=as.character(target_metadata1$Target_index)
fit_df=merge(fit_df,target_metadata[,c("Target_index","Region_ID")],by="Target_index")
p <- ggplot(fit_df, aes(X1,X2))
p + geom_point()+ geom_text(aes(X1, X2, label=row.names(fit[["points"]])))
p+ geom_point(aes(colour = Region_ID),cex=5) +
  scale_colour_gradientn(colours = terrain.colors(10))
