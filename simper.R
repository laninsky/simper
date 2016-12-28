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
gp_row <- which(names(gp_data)==group_name)
gp_names <- unique(gp_data[,gp_row])

for (i in 1:(dim(sp_data)[2])) {
  temp_sp <- NULL
  for (j in 1:((dim(sp_data)[1])-1)) {
  print(paste("up to row",j,"of",dim(sp_data)[1],"in total"))
  flush.console()
    for (k in (j+1):(dim(sp_data)[1])) {
      if(gp_data[j,gp_row]==gp_data[k,gp_row]) {
        next
      } else {
        bray_curt <- abs(sp_data[j,i]-sp_data[k,i])/(sp_data[j,i]+sp_data[k,i])
        temp_sp <- c(temp_sp,bray_curt)
      }
    }
  }
  sp_output[i,4] <- mean(temp_sp,na.rm=TRUE)
  sp_output[i,5] <- mean(temp_sp,na.rm=TRUE)/sd(temp_sp,na.rm=TRUE)
  sp_output[i,2] <- mean(sp_data[(which(gp_data[,1]==gp_names[1])),i],na.rm=TRUE)
  sp_output[i,3] <- mean(sp_data[(which(gp_data[,1]==gp_names[2])),i],na.rm=TRUE)   
  print(paste("done with species",i))
  flush.console()
}




