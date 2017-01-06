#Change these to the names of your various files. Data_file and expl_file should be a matrix of species abundances and a dataframe of explanatory variables, respectively. Group name is the column in the expl_file that has the group you want to perform the anosim over. 
data_file <- "fishies"
expl_file <- "fishtype"
group_name <- "Type"
permutations <- 9999

# Anosim on data with two groups
sp_data <- read.table(data_file,stringsAsFactors=FALSE,header=TRUE)
gp_data <- read.table(expl_file,stringsAsFactors=FALSE,header=TRUE)

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

sp_output_between <- NULL
sp_output_within1 <- NULL
sp_output_within2 <- NULL
gp_col <- which(names(gp_data)==group_name)
gp_names <- unique(gp_data[,gp_col])

sp_data_1 <- sp_data[(which(gp_data[,gp_col]==gp_names[1])),]
sp_data_2 <- sp_data[(which(gp_data[,gp_col]==gp_names[2])),]
sp_data_1 <- data.matrix(sp_data_1)
sp_data_2 <- data.matrix(sp_data_2)

for (j in 1:(dim(sp_data_1)[1])) {
  temp_sp <- rowSums(abs(sweep(sp_data_2, 2, sp_data_1[j,], "-")))/(rowSums(sp_data_2)+sum(sp_data_1[j,]))
  temp_sp <- matrix(temp_sp,ncol=1)
  sp_output_between <- rbind(sp_output_between,temp_sp)
 }
write.table(sp_output_between,"sp_output_between",quote=FALSE,col.names=FALSE,row.names=FALSE)
rm(sp_output_between)

for (j in 1:((dim(sp_data_1)[1])-1)) {
  if(!(is.null(dim(sp_data_1[-1:-j,])[1]))) {
    temp_sp <- rowSums(abs(sweep(sp_data_1[-1:-j,], 2, sp_data_1[j,], "-")))/(rowSums(sp_data_1[-1:-j,])+sum(sp_data_1[j,]))
    temp_sp <- matrix(temp_sp,ncol=1)
    sp_output_within1 <- rbind(sp_output_within1,temp_sp)
    } else {
    temp_sp <- sum(abs(sp_data_1[-1:-j,]-sp_data_1[j,]))/(sum(sp_data_1[-1:-j,])+sum(sp_data_1[j,]))
    temp_sp <- matrix(temp_sp,ncol=1)
    sp_output_within1 <- rbind(sp_output_within1,temp_sp)
    }
 }
write.table(sp_output_within1,"sp_output_within1",quote=FALSE,col.names=FALSE,row.names=FALSE)
rm(sp_output_within1)

for (j in 1:(dim(sp_data_2)[1])) {
  if(!(is.null(dim(sp_data_2[-1:-j,])[1]))) {
    temp_sp <- rowSums(abs(sweep(sp_data_2[-1:-j,], 2, sp_data_2[j,], "-")))/(rowSums(sp_data_2[-1:-j,])+sum(sp_data_2[j,]))
    temp_sp <- matrix(temp_sp,ncol=1)
    sp_output_within2 <- rbind(sp_output_within2,temp_sp)
    } else {
    temp_sp <- sum(abs(sp_data_2[-1:-j,]-sp_data_2[j,]))/(sum(sp_data_2[-1:-j,])+sum(sp_data_2[j,]))
    temp_sp <- matrix(temp_sp,ncol=1)
    sp_output_within2 <- rbind(sp_output_within2,temp_sp)
    }
 }

temp1 <- rep("within",(dim(sp_output_within2)[1]))
sp_output_within2 <- cbind(temp1,sp_output_within2)  

rm(sp_data)
rm(sp_data_1)
rm(sp_data_2)
rm(zerorows)
rm(temp_sp)
rm(data_file)
rm(expl_file)

sp_output_within1 <- as.matrix(read.table("sp_output_within1",stringsAsFactors=FALSE))
temp1 <- rep("within",(dim(sp_output_within1)[1]))
sp_output_within1 <- cbind(temp1,sp_output_within1)

sp_output_within1 <- rbind(sp_output_within1,sp_output_within2)
rm(sp_output_within2)

sp_output_between <- as.matrix(read.table("sp_output_between",stringsAsFactors=FALSE))
temp1 <- rep("between",(dim(sp_output_between)[1]))
sp_output_between <- cbind(temp1,sp_output_between)

sp_output_within1 <- rbind(sp_output_within1,sp_output_between)
rm(sp_output_between)

dist_ranks <- rank(sp_output_within1[,2],ties.method="average")
sp_output_within1 <- cbind(sp_output_within1,dist_ranks)
write.table(sp_output_within1,"ranked_distances",quote=FALSE,col.names=FALSE,row.names=FALSE)

no_of_samples <- dim(gp_data)[1]
M <- ((no_of_samples)*(no_of_samples-1))/2

obs_R <- ((mean(as.numeric(sp_output_within1[(which(sp_output_within1[,1]=="between")),3])))-(mean(as.numeric(sp_output_within1[(which(sp_output_within1[,1]=="within")),3]))))/(M/2)

sp_1_group <- sum(gp_data[,gp_col]==gp_names[1])
sp_2_group <- sum(gp_data[,gp_col]==gp_names[2])

within1_count <- choose(sp_1_group,2)
within2_count <- choose(sp_2_group,2)
between_count <- sp_1_group*sp_2_group

within1_ref <- unlist(combn(sp_1_group,2))
between_ref <- matrix(ncol=2,nrow=between_count)
between_ref[,1] <- rep(1:sp_1_group, each = (between_count/sp_1_group)) 
between_ref[,2] <- rep(1:sp_2_group, (between_count/sp_2_group))
within2_ref <- unlist(combn(sp_2_group,2))

sim_R_matrix <- matrix(NA,ncol=1,nrow=permutations)

for (i in 1:permutations) {
  coords <- NULL
  new_sp_1 <- sample(no_of_samples, sp_1_group, replace=FALSE)
  sp_1_cat <- new_sp_1[new_sp_1<=sp_1_group]
  sp_2_cat <- new_sp_1[new_sp_1>sp_1_group]
  sp_2_cat <- sp_2_cat-sp_1_group
  if(length(sp_1_cat)>=1) {
  coords <- which((within1_ref[1,] %in% sp_1_cat)!=(within1_ref[2,] %in% sp_1_cat))#
    if(length(sp_2_cat)>=1) {
      coords <- c(coords,(which((between_ref[,1] %in% sp_1_cat)!=(between_ref[,2] %in% sp_2_cat))+within1_count+within2_count))
    }
  }
  if(length(sp_2_cat)>=1) {
    coords <- c(coords,(which((within2_ref[1,] %in% sp_2_cat)!=(within2_ref[2,] %in% sp_2_cat))+within1_count))
  }
  sim_R_matrix[i,1] <- ((mean(as.numeric(sp_output_within1[coords,3])))-(mean(as.numeric(sp_output_within1[-coords,3]))))/(M/2)
  print(paste("done with permutation",i))
  flush.console()
}

pvalue <- sum(sim_R_matrix>=obs_R)/permutations

print(paste("observed R value of",obs_R,"between",paste(gp_names,collapse=" & "),"with a p-value of",pvalue,"based on",permutations,"permutations"))
