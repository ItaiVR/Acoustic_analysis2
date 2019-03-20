##Import FR metadata
FR.metadata <- read.csv("~/kinneret modeling/Acoustic_analysis2/Data/Single target analysis/Aquarium/FR metadata.csv", stringsAsFactors=FALSE)
##Import data on track to target index
targets.metadata <- read.csv("~/kinneret modeling/Acoustic_analysis2/Data/Single target analysis/Aquarium/targets metadata.csv", stringsAsFactors=FALSE)

###Merge data to identify target index to FR metadata

# ##Split time column in 'targets.metadata'
# targets.metadata$Ping_time_n=unlist(lapply(strsplit(as.character(targets.metadata$Ping_time), "\\."), "[", 1))
# targets.metadata$Ping_mili_n=as.numeric(unlist(lapply(strsplit(as.character(targets.metadata$Ping_time), "\\."), "[", 2)))/10
# 
# ##Create merge column
# targets.metadata$merge_time=paste(targets.metadata$Ping_time_n,targets.metadata$Ping_mili_n,sep=".")
# FR.metadata$merge_time=paste(paste(" ",FR.metadata$Ping_time,sep=""),FR.metadata$Ping_milliseconds,sep=".")

##Merge
#round range
FR.metadata$Range=round(FR.metadata$Range,digits=4)
targets.metadata$Target_range=round(targets.metadata$Target_range,digit=4)
FR.metadata_n=merge(FR.metadata,targets.metadata[,c("Ping_number","Target_range","Region_ID")],by.x=c("Ping_index","Range"),by.y=c("Ping_number","Target_range"),all.y=F)

##Write to file
write.csv(FR.metadata_n,"~/kinneret modeling/Acoustic_analysis2/Data/Single target analysis/Aquarium/FR_cluster_index.csv",row.names = F)

