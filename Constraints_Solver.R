## What this script does:
## Checks to see if pairs of genes can linearly separate BRACA1 and BRACA2 samples using LP solvers. 
-----------------------## CARREGAR BIBLIOTECAS ##-------------------------------------------------------------------------------------------------------------------------
library(RCurl)
library(MASS)
library(lpSolve)
-----------------------## PEGAR E PREPARAR OS DADOS ##-------------------------------------------------------------------------------------------------------------------------
## Pegar o arquivo
url <- getURL("http://yeast.ime.usp.br/~ronaldo/ibi5031-2013/data_breast_cancer.csv")
data <- read.csv(text = url, header = T, row.names=1)

## Pegar somente pacientes BRACA1/BRACA2
braca <- c("P01","P02","P03","P04","P05","P06","P07","P08","P09","P10","P18","P19","P20","P21","P22")
data <- data[braca]

## Pegar as classes
labels_plot <- data[3227,]
labels <- data[3227,]
for (i in 1:ncol(labels)) {
  if (labels[i] == "BRACA1") {
    labels[i] = +1
  }else {
    labels[i] = -1}}
h
## Transformar classes em matriz de vetores
vector_y <- as.data.frame(labels)
vector_y <- t(vector_y)

# Pre-processar os dados
data <- data[1:3226,]
data_mod <- matrix(as.numeric(as.matrix(data)), nrow = nrow(data))
colnames(data_mod) <- colnames(data)
rownames(data_mod) <- rownames(data)
colunas <- colnames(data_mod)

## Pegar todas combinações de genes 2 a 2 
gene_list <- rownames(data_mod)
gene_combinations <- combn(gene_list, 2)

## Montar a tabela de resultados do MSE
gene_combinations.df <- matrix(0, nrow = 11, ncol = 6)
gene_combinations.df <- as.data.frame(gene_combinations.df)
colnames(gene_combinations.df) <- c("Gene1", "Gene2", "bad_points", "bias", "weight1", "weight2")

## Montar a tabela de resultados do solver para LP
solver_ls_genes.df <- matrix(0, nrow = 11, ncol = 6)
solver_ls_genes.df <- as.data.frame(solver_ls_genes.df)
colnames(solver_ls_genes.df) <- c("Gene1", "Gene2", "Weight1", "Weight2", "Bias", "Min_Tau")
-----------------------## FUNÇÕES CLASSIFICADOR MSE ##-------------------------------------------------------------------------------------------------------------------------
## Função para gerar a matrix X
matrix_vector_genes <- function(dados, gene1, gene2) {
  sample_size <- ncol(dados)
  gene1 <- gene_combinations[1,i]
  gene2 <- gene_combinations[2,i]
  X <- matrix(1, nrow = sample_size, ncol= 3)
  X[,2] <- dados[gene1,]
  X[,3] <- dados[gene2,]
  return(X)
}

## Função para calcular o MSE
mse_calc <- function(X, Y) {
  a <- ginv(t(X) %*% X) %*% t(X) %*% Y ## a is a vector containing bias, weight1 and weight2
  return(a)
}

## Função para contar pontos mal classificados
count_missclassifications <- function(a, X, Y) {
  sample_size <- nrow(X)
  counts <- 0  
  for (i in 1:sample_size) {
    sin_sample <- (a[1] + a[2] * X[i,2] + a[3] * X[i,3])
    if ((sin_sample > 0 & Y[i] != 1) |  ## se é maior que 0, e o label for diferente de 1, então conte como malclassificado
      (sin_sample < 0 & Y[i] == 1)) {   ## se é menor que 0, e o label for igual a 1, então conte como malclassificado também
      counts <- counts + 1
    }
  }
  return(counts)
}
-----------------------## FUNÇÕES VERIFICAR SE GENES SAO LS ##-------------------------------------------------------------------------------------------------------------------------
## Função para gerar a matriz G
matrix_g <- function(dados, gene1, gene2, Y) {
  sample_size <- ncol(dados)
  rows <- sample_size + 1
  gene1 <- gene_combinations[1,i]
  gene2 <- gene_combinations[2,i]
  G <- matrix(-1, nrow = rows, ncol= 4)
  for (j in 1:length(dados[gene1,])) {
    G[j,1] <- (-1 * Y[j]) * dados[gene1,j]
  }
  for (k in 1:length(dados[gene2,])) {
    G[k,2] <- (-1 * Y[k]) * dados[gene2,k]
    G[k,3] <- (-1 * Y[k])
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
  h <- t(h)
  return(h)
}

## Função para gerar a matrix c
matrix_c <- function() {
  c <- c(0, 0, 0, 1)
  return(c)
}

## Função do solver
lp_solver <- function(c, G, h) {
  dir <- rep("<=", 16)
  solution <- lp("min", c, G, dir, h)$solution
  print(solution)
}

-----------------------## RODAR AS ANÁLISES ##-------------------------------------------------------------------------------------------------------------------------
## Para cada combinação de genes
paired_genes <- 0
genes_0_counts <- 0
ls_genes <- 0
c <- matrix_c()
for (i in 1:ncol(gene_combinations)) {
  paired_genes <- paired_genes + 1
  
  ## Gerar a matrix X
  #X <- matrix_vector_genes(data_mod, gene_combinations[1,i], gene_combinations[2,i])
  
  ## Calcular o MSE
  #a <- mse_calc(X, vector_y)
  
  ## Calcular o número de pontos malclassificados
  #bad_points <- count_missclassifications(a, X, vector_y)
    
  #if (bad_points == 0) {
   # genes_0_counts <- genes_0_counts + 1
    ## Inserir os resultados na tabela
    #gene_combinations.df[genes_0_counts,]$Gene1 <- gene_combinations[1,i]
    #gene_combinations.df[genes_0_counts,]$Gene2 <- gene_combinations[2,i]
    #gene_combinations.df[genes_0_counts,]$bad_points <- bad_points
    #gene_combinations.df[genes_0_counts,]$bias <- a[1]
    #gene_combinations.df[genes_0_counts,]$weight1 <- a[2]
    #gene_combinations.df[genes_0_counts,]$weight2 <- a[3]
    #}
  #else {
    ## Criar a matriz G
    G <- matrix_g(data_mod, gene_combinations[1,i], gene_combinations[2,i], vector_y)
      
    ## Criar a matriz h
    h <- matrix_h(data_mod)
      
    ## Resolver o LP
    solution <- lp_solver(c, G, h)
      
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
  print(paste("Genes with 0 points missclassified:",genes_0_counts))
}

rm(solver_ls_genes.df)
gene_combinations.df <- gene_combinations.df[1:10,]

## Criar a tabela em txt
write.table(gene_combinations.df, "Linear.Separable.Gene.Pairs.txt", quote = F, row.names = F, col.names = T, sep = "\t")
-----------------------## PLOTAR OS RESULTADOS ##-------------------------------------------------------------------------------------------------------------------------
rownames(labels_plot) <- "Group"
labels_plot <- t(labels_plot)
for (i in 1:nrow(gene_combinations.df)) {
  gene1 <- paste(gene_combinations.df[i,1], "$", sep = "")
  gene2 <- paste(gene_combinations.df[i,2], "$", sep = "")
  x <- grep(gene1, rownames(data_mod))
  y <- grep(gene2, rownames(data_mod))
  pair_plot <- rbind(data_mod[y,], data_mod[x,])
  pair_plot <- t(pair_plot)
  pair_plot <- cbind(pair_plot, labels_plot)
  pair_plot <- as.data.frame(pair_plot)
  b <- gene_combinations.df[i,4]
  w1 <- gene_combinations.df[i,5]
  w2 <- gene_combinations.df[i,6]
  pair_plot$V1 <- as.numeric(as.character(pair_plot$V1))
  pair_plot$V2 <- as.numeric(as.character(pair_plot$V2))
  minRange <- min(pair_plot$V1,pair_plot$V2)
  maxRange <- max(pair_plot$V1,pair_plot$V2)
  jpeg(paste(i, "_", gene_combinations.df[i,1], "_", gene_combinations.df[i,2], ".jpg", sep = ""))
  plot(pair_plot$V1, pair_plot$V2, col = pair_plot$Group, pch = 19, xlab= gene_combinations.df[i,1], ylab= gene_combinations.df[i,2])
  curve(expr = ((-b -w2 * x) / w1), from = minRange, to = maxRange, add = TRUE, col = "blue", lty = 2)
  legend("bottomright", pch=c(19,19), col=c("black", "red"), c("BRACA1", "BRACA2"), bty="n",  box.col="black", box.lwd = 0)
  dev.off()
}


