## What this script does:
## Checks to see if pairs of genes can linearly separate BRACA1 and BRACA2 samples using LP solvers. 
###-----------------------## LOAD THE LIBRARIES ##-------------------------------------------------------------------------------------------------------------------------
library(RCurl) 
library(MASS)
library(lpSolve)
###-----------------------## GET AND PREPROCESS THE DATA ##-------------------------------------------------------------------------------------------------------------------------
## Get the microarray data
url <- getURL("http://yeast.ime.usp.br/~ronaldo/ibi5031-2013/data_breast_cancer.csv")
data <- read.csv(text = url, header = T, row.names=1)

## Filter to include only the BRACA1/BRACA2 patients
braca <- c("P01","P02","P03","P04","P05","P06","P07","P08","P09","P10","P18","P19","P20","P21","P22")
data <- data[braca]

## Get the classes
labels_plot <- data[3227,]
labels <- data[3227,]
for (i in 1:ncol(labels)) {
  if (labels[i] == "BRACA1") {
    labels[i] = +1
  }else {
    labels[i] = -1}}

## Transform the classes in a matrix of vectors
vector_y <- as.data.frame(labels)
vector_y <- t(vector_y)

# Pre-process the data
data <- data[1:3226,]
data_mod <- matrix(as.numeric(as.matrix(data)), nrow = nrow(data))
colnames(data_mod) <- colnames(data)
rownames(data_mod) <- rownames(data)
colunas <- colnames(data_mod)

## Get all possible 2 by 2 gene combinations
gene_list <- rownames(data_mod)
gene_combinations <- combn(gene_list, 2)

## Create a table with the results of linear programming solver 
solver_ls_genes.df <- matrix(0, nrow = 11, ncol = 6)
solver_ls_genes.df <- as.data.frame(solver_ls_genes.df)
colnames(solver_ls_genes.df) <- c("Gene1", "Gene2", "Weight1", "Weight2", "Bias", "Min_Tau")

###-----------------------## FUNCTIONS TO CHECK IF THE GENE PAIRS ARE LINEARLY SEPARABLE ##-------------------------------------------------------------------------------------------------------------------------
## Função para gerar a matriz G
matrix_g <- function(dados, gene1, gene2, Y) {
  sample_size <- ncol(dados)
  rows <- sample_size + 1
  gene1 <- gene_combinations[1,i]
  gene2 <- gene_combinations[2,i]
  G <- matrix(1, nrow = rows, ncol= 4)
  for (j in 1:length(dados[gene1,])) {
    G[j,2] <- (-1 * Y[j]) * dados[gene1,j]
  }
  for (k in 1:length(dados[gene2,])) {
    G[k,3] <- (-1 * Y[k]) * dados[gene2,k]
    G[k,1] <- (-1 * Y[k])
  }
  G[,4] <- -1
  G[rows,1:3] <- 0
  return(G)
}

## Função para gerar a matrix h
matrix_h <- function(dados) {
  sample_size <- ncol(dados)
  h <- matrix(-1, nrow = sample_size + 1, ncol= 1)
  h[sample_size + 1,] = 0
  return(h)
}

## Função para gerar a matrix c
matrix_c <- function() {
  c <- as.matrix(c(0, 0, 0, 1), ncol=1)
  return(c)
}

## Função do solver
lp_solver <- function(c, h, G) {
  solution <- solveLP(cvec = c, bvec = h, Amat = G, maximum = F, const.dir =  rep("<=", 16))$solution
  print(solution)
}

###-----------------------## RUN THE ANALYSIS ##-------------------------------------------------------------------------------------------------------------------------
## For each gene pair combination
paired_genes <- 0
genes_0_counts <- 0
ls_genes <- 0
c <- matrix_c()
for (i in 1:ncol(gene_combinations)) {
  paired_genes <- paired_genes + 1
  
  ## Criar a matriz G
  G <- matrix_g(data_mod, gene_combinations[1,i], gene_combinations[2,i], vector_y)
      
  ## Criar a matriz h
  h <- matrix_h(data_mod)
      
  ## Resolver o LP
  solution <- lp_solver(c, h, G)
      
  if ( solution[4] == 0) {
    ls_genes <- ls_genes + 1
    print(paste("##Pair is linearly separable:##", ls_genes))
    ## Inserir os dados na tabela
    solver_ls_genes.df[ls_genes,]$Gene1 <- gene_combinations[1,i]
    solver_ls_genes.df[ls_genes,]$Gene2 <- gene_combinations[2,i]
    solver_ls_genes.df[ls_genes,]$Weight1 <- solution[1]
    solver_ls_genes.df[ls_genes,]$Weight2 <- solution[2]
    solver_ls_genes.df[ls_genes,]$Bias <- solution[3]
    solver_ls_genes.df[ls_genes,]$Min_Tau <- solution[4]
  }
  print(paste("Number of pairs evaluated:", paired_genes))
}

## Criar a tabela em txt
write.table(solver_ls_genes.df, "Linear.Separable.Gene.Pairs.txt", quote = F, row.names = F, col.names = T, sep = "\t")
