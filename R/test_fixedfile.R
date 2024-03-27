library(fastpolicytree)
library(policytree)

compare <- function(filename, actions, depth)
{
    
    data  <- as.matrix(read.table(filename, header=TRUE))
    nrows  <- dim(data)[1]
    ncols  <- dim(data)[2]
    p <- ncols-actions
    X <- data[,1:p]
    gamma  <- data[,(p+1):ncols]
    #print(gamma)
    times <- list()
    policies <- vector("list",2)
    trees <- vector("list",2)
    rewards  <- vector("list",2)
    for (method in 1:2)
    {
        if (method == 1 )
            time <- system.time(tree <- policy_tree(X, gamma, depth))
        else
            time <- system.time(tree <- fastpolicytree(X, gamma, depth))
        time <- as.vector(time)[3]
        policy <- predict(tree, X)
        #print(policy)
        reward  <- sum(gamma[cbind(1:nrows,policy)])
        times  <- c(times,time)
        trees[[method]]  <- tree
        policies[[method]]  <- policy
        rewards[[method]]  <- reward
        print(tree)
        print(reward)
    }
        
    if (!identical(policies[[1]],policies[[2]]))
    {
        print("Not identical")
        #print(trees[[1]])
        #print(trees[[2]])
        #print(policies[[1]])
        #print(policies[[2]])
        colnames(X) <- as.character(1:p)
        dat <-cbind(X,gamma)
        write.table(dat,file="tmpdat.txt",row.names=FALSE)
        quit("no")
    }
    
    cat(filename, actions, depth, times[[1]], times[[2]], "\n")
}

nactions  <- 2
depth  <- 3

args = commandArgs(trailingOnly = TRUE)
filename  <- args[1]
if(length(args) >= 2)
    nactions  <- as.integer(args[2])
if(length(args) >= 3)
    depth  <- as.integer(args[3])
    
compare(filename, nactions, depth)





