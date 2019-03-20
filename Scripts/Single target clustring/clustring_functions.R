###Load libraries
library(readr)
library(data.table)
library(tidyr)
library(dplyr)
library(emdist)
library(statip)
library(utils)
library(ggfortify)
library(ecodist)
library(plotly)
library(utils)
library(memisc)

##Find targets with high SD in response
sd_response=function(data_in,percent=0.9){
  data_in$SD=NA
  for (i in 1:nrow(data_in)){
    data_in[i,"SD"]=sd(data_in[i,8:(ncol(data_in)-1)])
  }
  data_sd=data_in%>%
    dplyr::select(Target_index,SD)%>%
    dplyr::filter(SD < quantile(SD, percent))
  
  
  return(data_sd)
}

##Function to prepare single target to earth mover distance format
sample_hist=function(data_in,sample_id,relative=F){
  sample_data=data_in %>%
    dplyr::filter(Target_index ==  sample_id) %>%
    dplyr::select(colnames(data_in)[8:ncol(data_in)])
  ##Transform data- first column is freq and second is response
  freq_numeric=as.numeric(unlist(lapply(strsplit(colnames(sample_data), "\\X"), "[", 2)))
  sample_out=data.frame(freq=freq_numeric,res=t(sample_data[1,]))
  ##Round to freq bins of 1 HRTZ
  sample_out$freq_bin=round(sample_out$freq,digits=0)
  sample_out=sample_out%>%
    dplyr::filter(freq_bin>=50 & freq_bin<=85)
  ##Get mean for each freq bin
  sample_bin=sample_out %>%
    dplyr::group_by(freq_bin) %>%
    dplyr::summarise(mean = mean(X1))
  if (isTRUE(relative)){
    sample_bin$mean=sample_bin$mean/min(sample_bin$mean)
  }else{}
  #Add weight for emd
  sample_bin$weight=1
  matrix_out=data.matrix(sample_bin[,c("weight","mean")])
  return(matrix_out)
  
}


##Apply earth mover distance
emd_freq_dist=function(data_in=survey_data_sub,sampleA,sampleB,relative_in=F){
  ##prepare the samples
  sample_1=sample_hist(data_in,sampleA,relative_in)
  sample_2=sample_hist(data_in,sampleB,relative_in)
  ##Run emd
  emd_dist=emd(sample_1,sample_2)
  # h1=unname(sample_1[,2])
  # h2=unname(sample_2[,2])
  # hellinger_dist=hellinger(h1,h2)
  return(emd_dist)
}


##Calculate distance between all sample pairs
emd_distance_matrix=function(data_in,percent_in=0.9,relative_in=F){
  ##Filter by standard deviation
  sd_df=sd_response(data_in,percent_in)
  ##Unique samples
  samples_u=unique( sd_df$Target_index)
  grid_samples=data.frame(t(combn(samples_u,2)))
  print(head(grid_samples))
  grid_samples$emd=NA
  ##Apply emd distance to each sample pair
  for (i in 1:nrow(grid_samples)){
    print(nrow(grid_samples)-i)
    grid_samples[i,"emd"]=emd_freq_dist(data_in,grid_samples[i,"X1"],grid_samples[i,"X2"],relative_in)
  }
  grid_samples$X1=as.character(grid_samples$X1)
  grid_samples$X2=as.character(grid_samples$X2)
  return(grid_samples)
}

###plot pair of histograms
hist_data=function(data_in,sample_id,relative){
  sample_data=data_in %>%
    dplyr::filter(Target_index ==  sample_id) %>%
    dplyr::select(colnames(data_in)[8:ncol(data_in)])
  ##Transform data- first column is freq and second is response
  freq_numeric=as.numeric(unlist(lapply(strsplit(colnames(sample_data), "\\X"), "[", 2)))
  sample_out=data.frame(freq=freq_numeric,res=t(sample_data[1,]))
  ##Round to freq bins of 1 HRTZ
  sample_out$freq_bin=round(sample_out$freq,digits=0)
  sample_out=sample_out%>%
    dplyr::filter(freq_bin>=50 & freq_bin<=85)
  ##Get mean for each freq bin
  sample_bin=sample_out %>%
    dplyr::group_by(freq_bin) %>%
    dplyr::summarise(mean = mean(X1))
  sample_bin$target=sample_id
  if (isTRUE(relative)){
    sample_bin$mean=sample_bin$mean/min(sample_bin$mean)
  }else{}
  return(sample_bin)
}


##Plot multiple histograms by target index
plot_pair=function(data_in,targets,emd_data,relative_in=F){
  hist_pair=data.frame()
  hist_sd=data.frame()
  for(t in targets){
    hist1=hist_data(data_in,t,relative_in)
    hist_sd_c=data.frame(sd=sd(hist1$mean))
    hist_sd_c$target=t
    hist_pair=rbind(hist_pair,hist1)
    hist_sd=rbind(hist_sd,hist_sd_c)}
  p <- ggplot(hist_pair, aes(freq_bin, mean))
  p= p + geom_point(aes(colour = factor(target)))
  #emd data
  emd_sub=emd_data[emd_data$X1 %in% targets & emd_data$X2 %in% targets,]
  return(list(p,hist_sd,emd_sub))
}

##Transform to distance matrix format
dist_matrix=function(distance_df){
  #create new matrix
  vals<-sort(unique(c(as.character(distance_df$X1), as.character(distance_df$X2))))
  nm <- matrix(NA, nrow=length(vals), ncol=length(vals), dimnames=list(vals, vals))
  diag(nm) <- 0   
  
  #fill
  nm[as.matrix(distance_df[, 1:2])] <- distance_df[,3]
  #fill reversed
  nm[as.matrix(distance_df[, 2:1])] <- distance_df[,3]
  nm
  return(nm)
}

