#Reading in the rockfish co-occurance file - your first species of interest should be in column one - there should be labels for columns (not rows)
# Change the file number to whatever you are reading in (make sure it is csv)
rocky <- read.csv(file="Rockstar_winter.csv",sep=",", stringsAsFactors=FALSE,header=FALSE)

#Defining some stuff that we'll refer to in the control loop below
maxcolumn <- dim(rocky)[2]
maxrows <- (dim(rocky)[1])-1
results_rocky <- NULL
total_results_rocky <- NULL
rock_results_final <- NULL

#Reworking the data into a matrix so it can be manipulated
rockmatrix <- as.numeric(as.matrix(rocky[2:(dim(rocky)[1]),1:maxcolumn]))
dim(rockmatrix) <- c((dim(rocky)[1])-1,maxcolumn)

# So this looks evil, but basically it is saying, compare the ith column, to every column 1 through k
# For the comparison to each kth column, look down the rows. At each jth row, see whether our fishies
# were present or not ##if(rockmatrix[j,i]==1)##. If they are present, then we see whether the
# fishies in column k were also present ##if(rockmatrix[j,i]==rockmatrix[j,k]##
# if they are we count them using that count function
### I've just convereted this whole thing to for rather than while functions so double check results ###
for (i in 1:maxcolumn) {
for (k in (i+1):maxcolumn) {
count <- 0
for (j in 1:maxrows) {
if(rockmatrix[j,i]==1) {
if(rockmatrix[j,i]==rockmatrix[j,k]) {
count <- count + 1 }
}}
# For our results, we spit out the fishie we are interested in (i), the fishie we are comparing to (k), and how many times they co-occured in total (j)
# Then we up our counter for k, and go on to the next fish we want to compare i to.
rock_results_temp <- cbind(rocky[1,i],rocky[1,k],count)
rock_results_final <- rbind(rock_results_final,rock_results_temp)
}
# Then after cycling through all the fishie comparisons for fishie i, we up the i counter and go on to the next fish to do detailed comparisons for
}

write.csv(rock_results_final,file="rock_results.final.csv")
