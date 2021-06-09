library('Matrix')

# Get args
args = commandArgs(trailingOnly=TRUE)

# Transform .RDS in .MTX
temp <- readRDS(args[1])
writeMM(temp, args[2])
write.csv(rownames(temp), args[3])
write.csv(colnames(temp), args[4])
