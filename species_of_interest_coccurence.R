#Reading in the fish (or whatever!) presence/absence file - your species of interest should be in column one, with the other species in the remaining columns - there should be labels for columns (not rows)
# Change the file name to whatever you are reading in (make sure it is csv)
rocky <- read.csv(file="Rockstar_winter.csv",sep=",", stringsAsFactors=FALSE,header=FALSE)

#Defining some stuff that we'll refer to in the control loop below
maxcolumn <- dim(rocky)[2]
maxrows <- (dim(rocky)[1])-1
rock_results_final <- NULL

#Reworking the data into a matrix so it can be manipulated
rockmatrix <- as.numeric(as.matrix(rocky[2:(dim(rocky)[1]),1:maxcolumn]))
dim(rockmatrix) <- c((dim(rocky)[1])-1,maxcolumn)

# So this looks evil, but basically it is saying, compare the 1st column, to every column 2 through k
# For the comparison to each column, look down the rows. At each jth row, see whether our fishies of
# interest were present or not ##if(rockmatrix[j,1]==1)##. If they are present, then we see whether the
# fishies in column k were also present ##if(rockmatrix[j,1]==rockmatrix[j,k]##
# if they are we count them using that count function
k <- 1
for (k in 2:maxcolumn) {
  count <- 0
  for (j in 1:maxrows) {
    if(is.na(rockmatrix[j,1])) {
      break
    }
    if(rockmatrix[j,1]==1) {
      if(rockmatrix[j,1]==rockmatrix[j,k]) {
      count <- count + 1 }
      }
    }
# For our results, we spit out our fish of interest, the fishie we are comparing to (k), and how many times they co-occured in total (j)
  rock_results_temp <- cbind(rocky[1,1],rocky[1,k],count)
  rock_results_final <- rbind(rock_results_final,rock_results_temp)
}

write.csv(rock_results_final,file="rock_results.final.csv")
