library(fastpolicytree)
library(policytree)

compare <- function(s, n, p, actions, depth)
{
        set.seed(s)
        X <- matrix(sample(0:1, n*p, replace = T, prob = c(.5,.5)), nrow=n, ncol=p)
        W <- sample(seq(1:actions), n, replace = TRUE)-1
        Y <- X[,1] + X[,2] * (W == 1) + X[,3] * (W == max(actions)) + runif(n, min=0)

        cf <- grf::causal_forest(X, Y, W)
        gamma <- double_robust_scores(cf)

        times <- list()
        policies <- vector("list",3)
        trees <- vector("list",3)
        for (method in 1:2)
        {
            if (method == 1 )
                time <- system.time(tree <- policy_tree(X, gamma, depth))
            else
                time <- system.time(tree <- fastpolicytree(X, gamma, depth))
            time <- as.vector(time)[3]
            policy <- predict(tree, X)
            times  <- c(times,time)
            trees[[method]]  <- tree
            policies[[method]]  <- policy
            ##print(policy)
            ##print(tree)
        }
        
        if (!identical(policies[[1]],policies[[2]]))
        {
                print(trees[[1]])
                cat("method: ", method-1, "\n")
                print(trees[[method+1]])
                colnames(X) <- as.character(1:p)
                dat <-cbind(X,gamma)
                write.table(dat,file="tmpdat.txt",row.names=FALSE)
                quit("no")
        }

        cat(s, n, p, actions, depth, times[[1]], times[[2]], "\n")
    }

nseeds  <- 20
n <- 500
p  <- 30
nactions  <- 2
depth  <- 3

args = commandArgs(trailingOnly = TRUE)
if(length(args) >= 1)
    nseeds  <- as.integer(args[1])
if(length(args) >= 2)
    n  <- as.integer(args[2])
if(length(args) >= 3)
    p  <- as.integer(args[3])
if(length(args) >= 4)
    nactions  <- as.integer(args[4])
if(length(args) >= 5)
    depth  <- as.integer(args[5])
    
for (s in 1:nseeds)
    compare(s, n, p, nactions, depth)





