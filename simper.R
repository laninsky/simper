#simper_file <- (data_file,expl_file,group_name) {
#}
data_file <- "fishies"
expl_file <- "fishtype"
group_name <- "Type"

# Simper on data with two groups
sp_data <- read.table(data_file,stringsAsFactors=FALSE)
gp_data <- read.table(expl_file,stringsAsFactors=FALSE)

zerorows <- NULL
for (i in 1:(dim(sp_data)[1])) {
  if (sum(sp_data[i,1:(dim(sp_data)[2])])==0) {
    zerorows <- c(zerorows,i)
  }
}

if (length(zerorows)>0) {
  sp_data <- sp_data[-zerorows,]
  gp_data <- sp_data[-zerorows,]
}

sp_output <- matrix(NA,ncol=7,nrow=(dim(sp_data)[2]))
sp_output[,1] <- names(sp_data)
gp_col <- which(names(gp_data)==group_name)
gp_names <- unique(gp_data[,gp_col])

sp_data_1 <- sp_data[(which(gp_data[,gp_col]==gp_names[1])),]
sp_data_2 <- sp_data[(which(gp_data[,gp_col]==gp_names[2])),]
sp_data_1 <- data.matrix(sp_data_1)
sp_data_2 <- data.matrix(sp_data_2)

#Rewrite this to take into account the number of processors, print out as an Rscript to run in parallel.
for (i in 1:(dim(sp_data)[2])) {
  temp_sp <- NULL
  for (j in 1:(dim(sp_data_1)[1])) {
  bray_curt <- abs(sp_data_2[,i]-sp_data_1[j,i])/(rowSums(sp_data_2)+sum(sp_data_1[j,]))
  bray_curt <- bray_curt[!is.na(bray_curt)]
  temp_sp <- c(temp_sp,bray_curt)
  }
  sp_output[i,4] <- mean(temp_sp,na.rm=TRUE)
  sp_output[i,5] <- mean(temp_sp,na.rm=TRUE)/sd(temp_sp,na.rm=TRUE)
  sp_output[i,2] <- mean(sp_data_1[,i],na.rm=TRUE)
  sp_output[i,3] <- mean(sp_data_2[,i],na.rm=TRUE)   
  print(paste("done with species",i))
  flush.console()
}

sp_output <- sp_output[order(sp_output[,5], decreasing=TRUE),]
sp_output <- sp_output[order(sp_output[,4], decreasing=TRUE),]

for (i in 1:(dim(sp_data)[2])) {
  sp_output[i,6] <- as.numeric(sp_output[i,4])/sum(as.numeric(sp_output[,4]))*100
  if(i!=1) {
    sp_output[i,7] <- as.numeric(sp_output[(i-1),7])+as.numeric(sp_output[i,6])
    } else {
    sp_output[i,7] <- sp_output[i,6]
  }
}  

sp_output_name <- c("Name",paste(gp_names[1],".Av.Abund",sep=""),paste(gp_names[2],"Av.Abund",sep=""),"Av.Diss","Diss/SD", "Contrib%", "Cum.%")
sp_output <- rbind(sp_output_name,sp_output)
write.table(sp_output,"species_simper_contributions.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)

