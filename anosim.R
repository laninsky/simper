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

sp_output <- NULL
gp_col <- which(names(gp_data)==group_name)
gp_names <- unique(gp_data[,gp_col])

sp_data_1 <- sp_data[(which(gp_data[,gp_col]==gp_names[1])),]
sp_data_2 <- sp_data[(which(gp_data[,gp_col]==gp_names[2])),]
sp_data_1 <- data.matrix(sp_data_1)
sp_data_2 <- data.matrix(sp_data_2)

for (j in 1:(dim(sp_data_1)[1])) {
  temp_sp <- rowSums(abs(sweep(sp_data_2, 2, sp_data_1[j,], "-")))/(rowSums(sp_data_2)+sum(sp_data_1[j,]))
  temp_sp <- matrix(temp_sp,ncol=1)
  sp_output <- rbind(sp_output,temp_sp)
 }
  
  
