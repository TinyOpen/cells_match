## calculate the CPM
cpm.ori <- read.table('all genes for CPM and cutoff.txt', sep = "\t", header = T, row.names = 1)
dim(cpm.ori)

counts <- cpm.ori[-c(1,2),]
counts <- apply(counts, 2, as.numeric)

count.cpm <- apply(counts, 2, FUN = function(x){
  x/sum(x)* 1000000
})
dim(count.cpm)

count.cpm.log <- apply(count.cpm, 2, FUN = function(x){
  ifelse(log2(x) > 5, 1,0) 
})

dim(count.cpm.log)
write.csv(count.cpm.log, "count_cpm_log_compare_5.csv")

remove(counts)
remove(count.cpm)
remove(cpm.ori)
remove(small.bet.mat.pair)
remove(small.bet.mat.pair.2)
remove(big.bet.mat.pair.2)
remove(big.bet.mat.pair)

##### match A and B

head(a.gene.name)
head(row.names(cpm.ori))

all.gene.name <- row.names(cpm.ori)[-c(1,2)]
length(all.gene.name)
length(which(all.gene.name %in% a.gene.name))

a.genes.cpm <- count.cpm.log[which(all.gene.name %in% a.gene.name),]
row.names(a.genes.cpm) <- all.gene.name[which(all.gene.name %in% a.gene.name)]
dim(a.genes.cpm)
dim(a.genes.count)

b.genes.cpm <- count.cpm.log[which(all.gene.name %in% b.gene.name),]
row.names(b.genes.cpm) <- all.gene.name[which(all.gene.name %in% b.gene.name)]

write.csv(b.genes.cpm, "b_genes_cpm.csv")
write.csv(a.genes.cpm, "a_genes_cpm.csv")

save(a.genes.count, file = "a_genes_count.RData")
save(b.genes.count, file = "b_genes_count.RData")
a.genes.count <- a.genes.cpm
b.genes.count <- b.genes.cpm
################ 
