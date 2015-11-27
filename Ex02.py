import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
import cvxopt as opt
import itertools 

#Importar os dados
fileUrl = 'http://yeast.ime.usp.br/~ronaldo/ibi5031-2013/data_breast_cancer.csv'
cancerData = pd.read_csv(fileUrl)
X = cancerData.values[:,1:23]
X = X[0:len(X)-1].astype(np.float)
Y = cancerData.values[3226:3227,1:23]
BRACA2 = np.argwhere(Y[0,:] == 'BRACA2')[:,0]
BRACA1 = np.argwhere(Y[0,:] == 'BRACA1')[:,0]

#Criar as labels, 1 BRACA1, -1 BRACA2
d1 = len(X[0,BRACA1])
d2 = len(X[0,BRACA2])	
labels_braca1 = np.ones(d1) 
labels_braca2 = np.ones(d2)

##Pegar todas as combinacoes de genes par a par 
all_genes = cancerData.values[1:3226,0]
all_combinations = itertools.combinations(all_genes, 2)

#Montar tabela com os resultados
ls_gene_pairs = open("ls_gene_pairs.txt", "w")
ls_gene_pairs.write("Gene1")
ls_gene_pairs.write("\t")
ls_gene_pairs.write("Gene2")
ls_gene_pairs.write("\t")
ls_gene_pairs.write("w1")
ls_gene_pairs.write("\t")
ls_gene_pairs.write("w2")
ls_gene_pairs.write("\t")
ls_gene_pairs.write("b")
ls_gene_pairs.write("\t")
ls_gene_pairs.write("Min_Tau")
ls_gene_pairs.write("\n")

counter = 0
ls_genes_counts= 0

#Rodar o solver para todos os possiveis pares 
for combination in all_combinations:
	counter += 1
	print "Evaluating gene pair: " + str(counter)
	gene1_index = cancerData[cancerData["Genes"] == combination[0]].index.tolist()
	gene2_index = cancerData[cancerData["Genes"] == combination[1]].index.tolist()
	gene1_braca1 = X[gene1_index,BRACA1].astype(np.float)[:]
	gene2_braca1 = X[gene2_index,BRACA1].astype(np.float)[:]
	gene1_braca2 = X[gene1_index,BRACA2].astype(np.float)[:]
	gene2_braca2 = X[gene2_index,BRACA2].astype(np.float)[:]
	braca1 = np.asmatrix([-labels_braca1,-gene1_braca1,-gene2_braca1,-labels_braca1])
	braca2 = np.asmatrix([labels_braca2,gene1_braca2,gene2_braca2,-labels_braca2])
	braca1.transpose()
	braca2.transpose()
	dim = 2
	A3 = np.hstack([np.zeros(1), np.zeros(dim), -np.ones(1)])
	AA = np.vstack([braca1.transpose(),braca2.transpose(), A3])
	G = opt.matrix(AA)
	n = len(AA)
	np.ones(n).transpose()
	h = opt.matrix(-np.ones(n).transpose())
	h[n-1] = h[n-1]+1
	c = np.zeros(AA.shape[1]-1)
	c = np.asmatrix(np.hstack([c,np.ones(1)])).transpose()
	c = opt.matrix(c)
	sol = opt.solvers.lp(c,G,h)
	sol = (sol['x'])
	if sol[3] <= 0.01 and sol[3] >= 0:
		ls_genes_counts += 1
		print "Gene pair found, in total: " + str(ls_genes_counts) + " genes found"
		ls_gene_pairs.write(combination[0])
		ls_gene_pairs.write("\t")
		ls_gene_pairs.write(combination[1])
		ls_gene_pairs.write("\t")
		ls_gene_pairs.write(str(sol[1]))
		ls_gene_pairs.write("\t")
		ls_gene_pairs.write(str(sol[2]))
		ls_gene_pairs.write("\t")
		ls_gene_pairs.write(str(sol[0]))
		ls_gene_pairs.write("\t")
		ls_gene_pairs.write(str(sol[3]))
		ls_gene_pairs.write("\n")

ls_gene_pairs.close()
