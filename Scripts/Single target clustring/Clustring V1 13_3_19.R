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
###Load acostic survey single target table
survey_data <- read.csv(unz("Data/Single target analysis/Hikaru samples.zip", "Hikaru samples/three group data 13_3_19_csv.csv"), header = TRUE,
                 sep = ",") 

set.seed(201)
survey_data_sub=survey_data %>%
  filter(Target_index %in% sample(x=survey_data$Target_index, size=100, replace = FALSE))

##Find targets with high SD in response
sd_response=function(data_in,percent=0.9){
  data_in$SD=NA
  for (i in 1:nrow(data_in)){
    data_in[i,"SD"]=sd(data_in[i,8:(ncol(data_in)-1)])
  }
  data_sd=data_in%>%
    select(Target_index,SD)%>%
    filter(SD < quantile(SD, percent))
  
  
  return(data_sd)
}

sd_df=sd_response(survey_data_sub)
summary(sd_df$SD)


##Function to prepare single target to earth mover distance format
sample_hist=function(data_in,sample_id,relative=F){
  sample_data=data_in %>%
    filter(Target_index ==  sample_id) %>%
    select(colnames(data_in)[8:ncol(data_in)])
  ##Transform data- first column is freq and second is response
  freq_numeric=as.numeric(unlist(lapply(strsplit(colnames(sample_data), "\\X"), "[", 2)))
  sample_out=data.frame(freq=freq_numeric,res=t(sample_data[1,]))
  ##Round to freq bins of 1 HRTZ
  sample_out$freq_bin=round(sample_out$freq,digits=0)
  sample_out=sample_out%>%
    filter(freq_bin>=50 & freq_bin<=85)
  ##Get mean for each freq bin
  sample_bin=sample_out %>%
    group_by(freq_bin) %>%
    summarise(mean = mean(X1))
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


# emd_distance_matrix_dplyr=function(data_in){
#   ##Unique samples
#   samples_u=unique(data_in$Target_index)
#   grid_samples=data.frame(t(combn(samples_u,2)))
#   #print(head(grid_samples))
#   ##Apply emd distance to each sample pair
#   data_out=grid_samples %>% 
#     rowwise() %>% 
#     mutate(emd=emd_freq_dist(data_in,X1,X2)) %>% 
#     select(X1,X2,emd)
#   grid_samples$X1=as.character(grid_samples$X1)
#   grid_samples$X2=as.character(grid_samples$X2)
#   return(grid_samples)
# }


# emd_distance_matrix_docall=function(data_in){
#   ##Unique samples
#   samples_u=unique(data_in$Target_index)
#   grid_samples=data.frame(t(combn(samples_u,2)))
#   #print(head(grid_samples))
#   ##Apply emd distance to each sample pair
#   grid_samples_out=mapply(emd_freq_dist, grid_samples$X1,grid_samples$X2)
#   grid_samples$X1=as.character(grid_samples$X1)
#   grid_samples$X2=as.character(grid_samples$X2)
#   return(grid_samples)
# }

##TEMP
data_in=survey_data_sub
sampleA=grid_samples[i,"X1"]
sampleB=grid_samples[i,"X2"]
##TEMP  
  
hikaru_100_SD0.7=emd_distance_matrix(survey_data_sub,percent_in=0.7)

###plot pair of histograms
hist_data=function(data_in,sample_id,relative){
  sample_data=data_in %>%
    filter(Target_index ==  sample_id) %>%
    select(colnames(data_in)[8:ncol(data_in)])
  ##Transform data- first column is freq and second is response
  freq_numeric=as.numeric(unlist(lapply(strsplit(colnames(sample_data), "\\X"), "[", 2)))
  sample_out=data.frame(freq=freq_numeric,res=t(sample_data[1,]))
  ##Round to freq bins of 1 HRTZ
  sample_out$freq_bin=round(sample_out$freq,digits=0)
  sample_out=sample_out%>%
    filter(freq_bin>=50 & freq_bin<=85)
  ##Get mean for each freq bin
  sample_bin=sample_out %>%
    group_by(freq_bin) %>%
    summarise(mean = mean(X1))
  sample_bin$target=sample_id
  if (isTRUE(relative)){
    sample_bin$mean=sample_bin$mean/min(sample_bin$mean)
  }else{}
  return(sample_bin)
}
##
hist_data(data_in,5,relative=T)

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

#similar pair
plot_pair(survey_data_sub,c(1183,2811),hikaru_100)
#dissimilar
plot_pair(survey_data_sub,c(745,3141,3415),hikaru_100)

plot_pair(survey_data_sub,c(478,518,625,502,1949,1947,3141),hikaru_100)



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

##Apply function
emd_matrix=dist_matrix(hikaru_100_SD0.9)


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
p <- ggplot(fit_df, aes(X1,X2))
p + geom_point()+ geom_text(aes(X1, X2, label=row.names(fit[["points"]])))
p+ geom_point(aes(colour = Depth),cex=5) +
  scale_colour_gradientn(colours = terrain.colors(10))

###Add mean intensity
data_int_f=function(data_in){
  data_freq=data_in[,8:ncol(data_in)]
  data_mean=rowMeans(data_freq)
  data_out=cbind(data_in$Target_index,data_mean)
  colnames(data_out)=c("Target_index","mean_int")
  data_out=data.frame(data_out)
  data_out$rel_int=data_out$mean_int/min(data_out$mean_int)
  return(data_out)
}

data_mean_int=data_int_f(survey_data_sub)

##add to distance matrix
fit_df_3d=merge(fit_df,data_mean_int,by="Target_index")

###3D ploting
p_3d <- plot_ly(fit_df_3d, x = ~X1, y = ~X2, z = ~rel_int,
             marker = list(color = ~Depth, colorscale = c('#FFE1A1', '#683531'), showscale = TRUE)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'X1'),
                      yaxis = list(title = 'X2'),
                      zaxis = list(title = 'Size')))


p_3d
#####Bray curtis distance

##Prepare data
##Function to prepare single target to earth mover distance format
sample_hist_BC=function(data_in){
  data_out=data.frame()
  for (i in unique(data_in$Target_index)){
    print(i)
  sample_data=data_in %>%
    filter(Target_index ==  i) %>%
    select(colnames(data_in)[8:ncol(data_in)])
  ##Transform data- first column is freq and second is response
  sample_out=data.frame(freq=as.numeric(colnames(sample_data)),res=t(sample_data[1,]))
  ##Round to freq bins of 1 HRTZ
  sample_out$freq_bin=round(sample_out$freq,digits=0)
  ##Get mean for each freq bin
  sample_bin=sample_out %>%
    group_by(freq_bin) %>%
    summarise(mean = mean(X1))
  data_row=data.frame(t(sample_bin$mean))
  colnames(data_row)=sample_bin$freq_bin
  data_row$Target_index=i
  data_out=rbind(data_out,data_row)
}
  return(data_out)
  
}

BC_data=sample_hist_BC(survey_data_sub)


##Run bray curtis
BC_out=bcdist(BC_data)
fit <- cmdscale(BC_out,eig=TRUE, k=2) # k is the number of dim
fit # view results
##ggplot
fit_df=data.frame(fit$points)
#fit_df$Target_index=row.names(fit[["points"]])
#fit_df=merge(fit_df,survey_data_sub[,c("Target_index","Depth")],by="Target_index")
p <- ggplot(fit_df, aes(X1,X2))
p + geom_point()
p+ geom_point(aes(colour = Depth),cex=5) +
  scale_colour_gradientn(colours = terrain.colors(10))
