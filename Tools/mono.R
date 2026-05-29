# !/usr/bin/Rscript
# Author: Tao Xiong
# Date: 2025.04.16
# This script is used to check the monophyly of each clade


# 修改zoom函数
# 修改 ape::zoom 使其在左侧大图也显示 tip labels
# 修改 ape::zoom，左侧大图固定 cex = 0.5，右侧子图使用用户传入的 cex
assignInNamespace("zoom",
  function(phy, focus, subtree = FALSE, col = rainbow, ...) {
    if (!is.list(focus)) focus <- list(focus)
    n <- length(focus)
    for (i in 1:n) if (is.character(focus[[i]])) 
        focus[[i]] <- which(phy$tip.label %in% focus[[i]])
    if (is.function(col)) {
        col <- if (deparse(substitute(col)) == "grey") grey(1:n/n) else col(n)
    }
    ext <- vector("list", n)
    for (i in 1:n) 
        ext[[i]] <- drop.tip(phy, phy$tip.label[-focus[[i]]], 
                             subtree = subtree, rooted = TRUE)
    nc <- round(sqrt(n)) + 1
    nr <- ceiling(sqrt(n))
    M <- matrix(0, nr, nc)
    x <- c(rep(1, nr), 2:(n + 1))
    M[1:length(x)] <- x
    layout(M, c(1, rep(3/(nc - 1), nc - 1)))
    
    # 边颜色：高亮树枝
    colo <- rep("black", dim(phy$edge)[1])
    for (i in 1:n) colo[which.edge(phy, focus[[i]])] <- col[i]
    
    # tip 标签颜色：属于焦点组的标签使用与边相同的颜色
    tip.col <- rep("black", length(phy$tip.label))
    for (i in 1:n) {
        tip.col[focus[[i]]] <- col[i]
    }
    
    # 左侧大图：固定 cex=0.5，保留用户的其他参数，并加入 tip.color
    dot_args <- list(...)
    dot_args$cex <- NULL  # 移除用户可能的 cex 避免冲突
    left_args <- c(list(phy, edge.color = colo, show.tip.label = TRUE, 
                        cex = 0.5, tip.color = tip.col), dot_args)
    do.call(plot.phylo, left_args)
    
    # 右侧子图：完全使用用户参数
    for (i in 1:n) {
        right_args <- c(list(ext[[i]], edge.color = col[i]), list(...))
        do.call(plot.phylo, right_args)
    }
  }, "ape")


file<-commandArgs(TRUE)

library(ape)
library(dplyr)
library(phytools)

if(!dir.exists("mono")){
  dir.create("mono")
}

#read tree
tree <- read.tree(file)
tip <- as.data.frame(tree$tip.label)
colnames(tip) <- c("tiplables")

#get the subfamily
tip$subf <- sapply(strsplit(as.character(tip$tiplables), "_"),`[`,1)

#get the subfamily name list
subfamily <- sapply(strsplit(as.character(tip$tiplables), "_"),`[`,1) %>% unique()


df <- data.frame()


for(name in subfamily){
  tips <- tip$tiplables[tip$subf == name]
  mono <- is.monophyletic(tree,tips,reroot = !is.rooted(tree))
  if (mono == FALSE && length(tips)>2 ){
    print(paste0(name," is a not a monophyletic clade"))
    pdf(paste("./mono/",name,".pdf",sep = ""),height = 100,width = 100)
    is.monophyletic(tree, tips,reroot = !is.rooted(tree),plot=T,cex=10,lwd=50)
    dev.off()
    write.table(tips,paste("./mono/",name,".txt",sep = ""),row.names = F,sep = "\t",quote =F)
  }
}
