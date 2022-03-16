#!/usr/bin/env Rscript

# script name: ultram-dic_tree_loop.R
# run: Rscript ultram-dic_tree_loop.R

#Installing packages
install.packages("phytools")
library(phytools)

install.packages("phangorn")
library(phangorn)

install.packages("ape")
library(ape)

install.packages("picante")
library(picante)

# load all the trees
myTrees <- list.files(pattern="*order")
# loop in mytrees
for (t in myTrees) {
  #read tree
  tree <- read.tree(t)
  #do ultrametric tree
  ult.tree<-chronos(tree, lambda=0)
  # do dicotomous
  tree.dic <- multi2di(ult.tree, random=TRUE) 
  # save ultrametric tree to a pdf
  pdf(paste0(t,"_ultra_dic.pdf"), width = 8, height = 11)
  plot <- plotTree(tree.dic)
  dev.off()
  # save new ultrametric tree to a file
  write.tree(tree.dic, file = paste0(t,".ultra.dic.tree"), append = FALSE, digits = 10, tree.names = FALSE)
}

