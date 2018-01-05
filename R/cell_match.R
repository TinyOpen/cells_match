########### import data ##########
pairs <- read.table('data/A-B pair.txt')
a.genes <- read.table('data/A-gene-counts.txt', sep = "\t", header = T, row.names = 1)
b.genes <- read.table('data/B-gene-counts.txt', sep = "\t", header = T, row.names = 1)


##########  count  #################
a.genes.count <- a.genes[-c(1,2),]
a.genes.count <- apply(a.genes.count, 2, as.numeric)

b.genes.count <- b.genes[-c(1,2),]
b.genes.count <- apply(b.genes.count, 2, as.numeric)

row.names(a.genes.count) <- row.names(a.genes)[-c(1:2)]
row.names(b.genes.count) <- row.names(b.genes)[-c(1:2)]

cell.name <- colnames(a.genes.count)
a.gene.name <- row.names(a.genes.count)
b.gene.name <- row.names(b.genes.count)
n.cell <- dim(a.genes.count)[2]

colnames(a.genes.count) <- 1:n.cell
colnames(b.genes.count) <- 1:n.cell


#dim(a.genes.count)
#[1]  682 2734

########### group #################
group.ind <- b.genes[c(1,2),]
big.group.ind <- as.factor(as.vector(unlist(b.genes[1,])))
small.group.ind <- as.factor(as.vector(unlist(b.genes[2,])))

group.ind <- data.frame(big.group.ind = big.group.ind, small.group.ind = small.group.ind)
big.group.levels <- levels(big.group.ind)
small.group.levels <- levels(small.group.ind)


########### pairs ###############
class(pairs)
pairs <- as.character(unlist(pairs))
pairs.list <- strsplit(pairs, split="-")

# 2557 pairs

##########  Cell Match #############

a <- lapply(pairs.list, FUN = function(x){
  gene.name <- unlist(x)
  gene.name.a <- gene.name[1]
  gene.name.b <- gene.name[2]
  
  if(gene.name.b %in% b.gene.name & gene.name.a %in% a.gene.name){
    a.count <- a.genes.count[gene.name.a,]
    b.count <- b.genes.count[gene.name.b,]
    match <- as.numeric(a.count!=0 & b.count!=0)
    return(match)
  } else {
    match = "NULL"
    return(match)
  }
})

pair.nonexit <- which(a == "NULL")
pairs.list[60] # "ALOX5AP" "ALOX5"  
pairs.list[116] # "AREGB" "EGFR" 

which(gene.name.b == "EGFR") # integer(0)
which(gene.name.b == "ALOX5") # integer(0)
which(gene.name.a == "AREGB") # integer(0)
which(gene.name.a == "ALOX5AP") # integer(0)


## test pairs 1: "A2M"  "LRP1"
a.genes.count["A2M", 6] # 58
b.genes.count["LRP1", 6] # 1

gene.name <- unlist(pairs.list[2])

length(a)
length(pairs.list) # 2557
length(pair.nonexit) # 235
match.all <- a[-pair.nonexit]
length(match.all) # 2557 - 235 = 2322


cell.bet.mat <- sapply(1:n.cell, FUN = function(a.ind){
  cat("+")
  ### calculate cell a.ind's relationship with other cells 
  #a.2 <- sapply(pairs.list[-pair.nonexit], FUN = function(x){
  a.2 <- sapply(pairs.list[1:50], FUN = function(x){
    gene.name <- unlist(x)
    gene.name.a <- gene.name[1]
    gene.name.b <- gene.name[2]
    
    a.count <- a.genes.count[gene.name.a,]
    b.count <- b.genes.count[gene.name.b, a.ind]
    match <- as.numeric(b.count!=0 & a.count!=0)
    return(match)
  }) 
  ### dim(a.2) # 2734 cells * 2322 gene pairs
  
  a.ind.sum <- apply(a.2, 1, sum) # length(a.ind.mean [1] 2734
  return(a.ind.sum)
})


################## Big group ###################
### n * m ####
n.big <- length(big.group.levels)
big.bet.mat <- matrix(0, n.big, n.big)

for(i in 1:n.big){
  
  cell.i.ind <- which(group.ind[,1] == big.group.levels[i])
  
  for(j in 1:n.big){
    cell.j.ind <- which(group.ind[,1] == big.group.levels[j])
    sub.cell.bet.mat <- cell.bet.mat[cell.i.ind, cell.j.ind]
    big.bet.mat[i,j] <- sum(apply(sub.cell.bet.mat, 2, sum))
  }
}

library(d3heatmap)
d3heatmap(big.bet.mat, dendrogram = "none")

################# Small group ####################


n.small <- length(small.group.levels)
small.bet.mat <- matrix(0, n.small, n.small)

for(i in 1:n.small){
  cell.i.ind <- which(group.ind[,2] == small.group.levels[i])
  
  for(j in 1:n.small){
    cell.j.ind <- which(group.ind[,2] == small.group.levels[j])
    sub.cell.bet.mat <- cell.bet.mat[cell.i.ind, cell.j.ind]
    small.bet.mat[i,j] <- sum(apply(sub.cell.bet.mat, 2, sum))
  }
}

d3heatmap(small.bet.mat, dendrogram = "none")

################ pairs x group ################
big.bet.mat.pair <- matrix(0, 2322, n.big * n.big)

m = 1
for(i in 1:n.big){
  cell.i.ind <- which(group.ind[,1] == big.group.levels[i])
  
  for(j in 1:n.big){
    cat("+", m)
    cell.j.ind <- which(group.ind[,1] == big.group.levels[j])

    a <- sapply(pairs.list[-pair.nonexit], FUN = function(x){
      gene.name <- unlist(x)
      gene.name.a <- gene.name[1]
      gene.name.b <- gene.name[2]
      
      a.count <- a.genes.count[gene.name.a, cell.i.ind] # a 的基因数据
      b.count <- b.genes.count[gene.name.b, cell.j.ind] # b 的基因数据
      
      a.count <- as.numeric(a.count > 0)
      b.count <- as.numeric(b.count > 0)
      
      match <- sum(matrix(a.count, ncol = 1) %*% matrix(b.count, nrow = 1)) 
      return(match)
    })
    big.bet.mat.pair[,m] <- a
    m = m+1
  }
}


small.bet.mat.pair <- matrix(0, 2322, n.small * n.small)

m = 1
for(i in 1:n.small){
  cell.i.ind <- which(group.ind[,2] == small.group.levels[i])
  
  for(j in 1:n.small){
    cat("+", m)
    cell.j.ind <- which(group.ind[,2] == small.group.levels[j])

    a <- sapply(pairs.list[-pair.nonexit], FUN = function(x){
      gene.name <- unlist(x)
      gene.name.a <- gene.name[1]
      gene.name.b <- gene.name[2]
      
      a.count <- a.genes.count[gene.name.a, cell.i.ind] # a 的基因数据
      b.count <- b.genes.count[gene.name.b, cell.j.ind] # b 的基因数据
      
      a.count <- as.numeric(a.count > 0)
      b.count <- as.numeric(b.count > 0)
      
      match <- sum(matrix(a.count, ncol = 1) %*% matrix(b.count, nrow = 1)) 
      return(match)
    })
    small.bet.mat.pair[,m] <- a
    m = m+1
  }
}


colnames(big.bet.mat.pair) <- paste( 
                                    rep(big.group.levels, each = 15), 
                                    "-", 
                                    rep(big.group.levels, 15), 
                                    sep = "")

colnames(small.bet.mat.pair) <- paste(rep(small.group.levels, each = 33), 
                                      "-", rep(small.group.levels, 33), sep = "")

rownames(small.bet.mat.pair) <- pairs[-pair.nonexit]
rownames(big.bet.mat.pair) <- pairs[-pair.nonexit]

write.csv(small.bet.mat.pair, "small_bet_mat_pair.csv")
write.csv(big.bet.mat.pair, "big_bet_mat_pair.csv")


big.factor <- rep(0, 15 * 15)
m = 1
for(i in 1:15){
  cell.i.ind <- which(group.ind[,1] == big.group.levels[i])
  
  for(j in 1:15){
    cell.j.ind <- which(group.ind[,1] == big.group.levels[j])
    big.factor[m] <- length(cell.i.ind) * length(cell.j.ind)
    m = m + 1
  }
}
big.factor <- as.matrix(big.factor, ncol = 1)
big.bet.mat.pair.2 <- apply(big.bet.mat.pair, 1,
                              FUN = function(x) {
                                xx <- as.matrix(x, ncol = 1)
                                round(xx/big.factor, 4)} 
)

small.factor <- rep(0, n.small * n.small)
m = 1
for(i in 1:n.small){
  cell.i.ind <- which(group.ind[,2] == small.group.levels[i])
  
  for(j in 1:n.small){
    cell.j.ind <- which(group.ind[,2] == small.group.levels[j])
    small.factor[m] <- length(cell.i.ind) * length(cell.j.ind)
    m <- m + 1
  }
}
small.factor <- as.matrix(small.factor, ncol = 1)
small.bet.mat.pair.2 <- apply(small.bet.mat.pair, 1,
                              FUN = function(x) {
                                xx <- as.matrix(x, ncol = 1)
                                round(xx/small.factor, 4)} 
                              )
small.bet.mat.pair.2 <- t(small.bet.mat.pair.2)
big.bet.mat.pair.2 <- t(big.bet.mat.pair.2)



colnames(big.bet.mat.pair.2) <- paste( 
  rep(big.group.levels, each = 15), 
  "-", 
  rep(big.group.levels, 15), 
  sep = "")

colnames(small.bet.mat.pair.2) <- paste(rep(small.group.levels, each = 33), 
                                      "-", rep(small.group.levels, 33), sep = "")

write.csv(small.bet.mat.pair.2, "small_bet_mat_pair_2.csv")
write.csv(big.bet.mat.pair.2, "big_bet_mat_pair_2.csv")

matrix(1:4, 2) / matrix(3:6,2)
