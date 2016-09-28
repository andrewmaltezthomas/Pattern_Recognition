import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import math
import cvxopt as opt
import itertools
from pylab import *
from matplotlib.pyplot import savefig

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

#Importar os genes linearmente separaveis 
ls_genes = open("ls_gene_pairs.txt", "r")
lines_lp = ls_genes.readlines()[1:]

#Criar a tabela da maxima margem  
maximum_margin_genes = open("maximum_margin_gene_pairs.txt", "w")
maximum_margin_genes.write("Gene1")
maximum_margin_genes.write("\t")
maximum_margin_genes.write("Gene2")
maximum_margin_genes.write("\t")
maximum_margin_genes.write("w1")
maximum_margin_genes.write("\t")
maximum_margin_genes.write("w2")
maximum_margin_genes.write("\t")
maximum_margin_genes.write("b")
maximum_margin_genes.write("\t")
maximum_margin_genes.write("Maximum_Margin")
maximum_margin_genes.write("\n")

#Calcular a margem para cada par LS
counter = 0
for line in lines_lp:
	counter += 1
	line = line.rstrip()
        line = line.split("\t")
	print "Evaluating LS gene pair: " + str(counter)
        gene1_index = cancerData[cancerData["Genes"] == line[0]].index.tolist()
        gene2_index = cancerData[cancerData["Genes"] == line[1]].index.tolist() 
        gene1_braca1 = X[gene1_index,BRACA1].astype(np.float)[:]
        gene2_braca1 = X[gene2_index,BRACA1].astype(np.float)[:]
        gene1_braca2 = X[gene1_index,BRACA2].astype(np.float)[:]
        gene2_braca2 = X[gene2_index,BRACA2].astype(np.float)[:]
        braca1 = np.asmatrix([-labels_braca1,-gene1_braca1,-gene2_braca1])
        braca2 = np.asmatrix([labels_braca2,gene1_braca2,gene2_braca2])
        A1 = braca1.transpose()
        A2 = braca2.transpose()
        # Matriz G
	G = opt.matrix(np.vstack([A1,A2]))
        d = 3
        # Matriz P 
	P = opt.matrix(np.identity(d))
	P[0,0] = 0.0
	# Matriz q 
	q = opt.matrix(0.0,(d,1))
        h1 = np.asmatrix(np.ones(d1)).transpose()
        h2 = np.asmatrix(np.ones(d2)).transpose()
        # Matriz h 
	h = opt.matrix(np.vstack([-h1,-h2]))
        sol = opt.solvers.qp (P, q, G, h, None, None)
	sol = sol['x']
        b = sol[0]
        w1 = sol[1]
        w2 = sol[2]
       	Margin = 1.0 / math.sqrt(w1 ** 2 + w2 ** 2)
        maximum_margin_genes.write(line[0])
        maximum_margin_genes.write("\t")
       	maximum_margin_genes.write(line[1])
        maximum_margin_genes.write("\t")
        maximum_margin_genes.write(str(w1))
        maximum_margin_genes.write("\t")
        maximum_margin_genes.write(str(w2))
        maximum_margin_genes.write("\t")
        maximum_margin_genes.write(str(b))
        maximum_margin_genes.write("\t")
        maximum_margin_genes.write(str(Margin))
        maximum_margin_genes.write("\n")

maximum_margin_genes.close()
ls_genes.close()

#Ordenando os pares pelo tamanho da margem
df = pd.read_csv("maximum_margin_gene_pairs.txt", sep ="\t")
df = df.sort('Maximum_Margin', ascending = False)
df_100 = df[0:100]
df_10 = df[0:10]
df_100.to_csv('Top_100_Maximum_Margin.csv', index=False, sep= "\t")

#Plotar os resultados para os 10 pares com maior margem
k = 0
for index, row in df_10.iterrows():
        k += 1
        figure(k)
        i = cancerData[cancerData["Genes"] == row[0]].index.tolist()
        j = cancerData[cancerData["Genes"] == row[1]].index.tolist()
        i = int(i[0])
        j = int(j[0])
        xmin = math.floor(float(min([np.min(X[i,BRACA1]),np.min(X[i,BRACA2])])))
        xmax = math.ceil(float(max([np.max(X[j,BRACA1]),np.max(X[j,BRACA2])])))
        xspace = np.linspace(xmin,xmax)
        bias = float(row[4])
        w1 = float(row[2])
        w2 = float(row[3])
        a = - (w1 / w2)
        b = - (bias / w2)
        plt.scatter(list(X[i,BRACA1]),list(X[j,BRACA1]), color='r')
        plt.scatter(list(X[i,BRACA2]),list(X[j,BRACA2]), color='b')
        plt.legend(['BRACA1', 'BRACA2'])
        yspace = a * xspace + b
        yyspace = a * xspace + b - 1.0 / w2
        yyyspace = a * xspace + b + 1.0 / w2
        plt.plot(xspace, yyspace, 'k--')
        plt.plot(xspace, yspace, 'k-')import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import math
import cvxopt as opt
import itertools
from pylab import *
from matplotlib.pyplot import savefig

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

#Importar os genes linearmente separaveis 
ls_genes = open("ls_gene_pairs.txt", "r")
lines_lp = ls_genes.readlines()[1:]

#Criar a tabela da maxima margem  
maximum_margin_genes = open("maximum_margin_gene_pairs.txt", "w")
maximum_margin_genes.write("Gene1")
maximum_margin_genes.write("\t")
maximum_margin_genes.write("Gene2")
maximum_margin_genes.write("\t")
maximum_margin_genes.write("w1")
maximum_margin_genes.write("\t")
maximum_margin_genes.write("w2")
maximum_margin_genes.write("\t")
maximum_margin_genes.write("b")
maximum_margin_genes.write("\t")
maximum_margin_genes.write("Maximum_Margin")
maximum_margin_genes.write("\n")

#Calcular a margem para cada par LS
counter = 0
for line in lines_lp:
	counter += 1
	line = line.rstrip()
        line = line.split("\t")
	print "Evaluating LS gene pair: " + str(counter)
        gene1_index = cancerData[cancerData["Genes"] == line[0]].index.tolist()
        gene2_index = cancerData[cancerData["Genes"] == line[1]].index.tolist() 
        gene1_braca1 = X[gene1_index,BRACA1].astype(np.float)[:]
        gene2_braca1 = X[gene2_index,BRACA1].astype(np.float)[:]
        gene1_braca2 = X[gene1_index,BRACA2].astype(np.float)[:]
        gene2_braca2 = X[gene2_index,BRACA2].astype(np.float)[:]
        braca1 = np.asmatrix([-labels_braca1,-gene1_braca1,-gene2_braca1])
        braca2 = np.asmatrix([labels_braca2,gene1_braca2,gene2_braca2])
        A1 = braca1.transpose()
        A2 = braca2.transpose()
        # Matriz G
	G = opt.matrix(np.vstack([A1,A2]))
        d = 3
        # Matriz P 
	P = opt.matrix(np.identity(d))
	P[0,0] = 0.0
	# Matriz q 
	q = opt.matrix(0.0,(d,1))
        h1 = np.asmatrix(np.ones(d1)).transpose()
        h2 = np.asmatrix(np.ones(d2)).transpose()
        # Matriz h 
	h = opt.matrix(np.vstack([-h1,-h2]))
        sol = opt.solvers.qp (P, q, G, h, None, None)
	sol = sol['x']
        b = sol[0]
        w1 = sol[1]
        w2 = sol[2]
       	Margin = 1.0 / math.sqrt(w1 ** 2 + w2 ** 2)
        maximum_margin_genes.write(line[0])
        maximum_margin_genes.write("\t")
       	maximum_margin_genes.write(line[1])
        maximum_margin_genes.write("\t")
        maximum_margin_genes.write(str(w1))
        maximum_margin_genes.write("\t")
        maximum_margin_genes.write(str(w2))
        maximum_margin_genes.write("\t")
        maximum_margin_genes.write(str(b))
        maximum_margin_genes.write("\t")
        maximum_margin_genes.write(str(Margin))
        maximum_margin_genes.write("\n")

maximum_margin_genes.close()
ls_genes.close()

#Ordenando os pares pelo tamanho da margem
df = pd.read_csv("maximum_margin_gene_pairs.txt", sep ="\t")
df = df.sort('Maximum_Margin', ascending = False)
df_100 = df[0:100]
df_10 = df[0:10]
df_100.to_csv('Top_100_Maximum_Margin.csv', index=False, sep= "\t")

#Plotar os resultados para os 10 pares com maior margem
k = 0
for index, row in df_10.iterrows():
        k += 1
        figure(k)
        i = cancerData[cancerData["Genes"] == row[0]].index.tolist()
        j = cancerData[cancerData["Genes"] == row[1]].index.tolist()
        i = int(i[0])
        j = int(j[0])
        xmin = math.floor(float(min([np.min(X[i,BRACA1]),np.min(X[i,BRACA2])])))
        xmax = math.ceil(float(max([np.max(X[j,BRACA1]),np.max(X[j,BRACA2])])))
        xspace = np.linspace(xmin,xmax)
        bias = float(row[4])
        w1 = float(row[2])
        w2 = float(row[3])
        a = - (w1 / w2)
        b = - (bias / w2)
        plt.scatter(list(X[i,BRACA1]),list(X[j,BRACA1]), color='r')
        plt.scatter(list(X[i,BRACA2]),list(X[j,BRACA2]), color='b')
        plt.legend(['BRACA1', 'BRACA2'])
        yspace = a * xspace + b
        yyspace = a * xspace + b - 1.0 / w2
        yyyspace = a * xspace + b + 1.0 / w2
        plt.plot(xspace, yyspace, 'k--')
        plt.plot(xspace, yspace, 'k-')
        plt.plot(xspace, yyyspace, 'k--')
        plt.title("Maximum Linear Margin Scatter Plot")
        plt.xlabel(row[0])
        plt.ylabel(row[1])
        plt.savefig('Maximum_Margin_Pair_' + str(k) + '.png')
        plt.plot(xspace, yyyspace, 'k--')
        plt.title("Maximum Linear Margin Scatter Plot")
        plt.xlabel(row[0])
        plt.ylabel(row[1])
        plt.savefig('Maximum_Margin_Pair_' + str(k) + '.png')
