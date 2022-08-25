library(policytree)

data  <- read.csv("IFLS_orig.txt",header=TRUE)
x <- subset(data,select=-c(A,Y,scores.DML,scores.DR,scores.CF))
nrows  <- dim(data)[1]
gammas <- cbind(data[,"scores.DML"],numeric(nrows))
tree  <- policy_tree(x,gammas,2)
tree
