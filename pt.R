#!/usr/bin/env Rscript

library(policytree)

args = commandArgs(trailingOnly=TRUE)
datafile  <- args[1]
nactions  <- as.integer(args[2])
if (length(args) > 2) {
    depth  <- as.integer(args[3])
} else depth <- 2
data  <- read.table(datafile,header=TRUE)
nrows  <- dim(data)[1]
ncols  <- dim(data)[2]
x <- data[,(1:(ncols-nactions))]
gammas <- data[,((ncols-nactions+1):ncols)]
tree  <- policy_tree(x,gammas,depth,TRUE)
tree
