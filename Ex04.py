from sklearn import cross_validation
from sklearn import svm
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
import operator

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
labels_braca2 = -np.ones(d2)

#Importar os genes linearmente separaveis 
ls_genes = open("ls_gene_pairs.txt", "r")
lines_lp = ls_genes.readlines()[1:]

#Custos
costs = [1e-3,1e-2,1e-1,1,1e1,1e2,1e3]

#Criar a tabela do soft margem  
soft_margin_costs = open("soft_margin_costs.txt", "w")
soft_margin_costs.write("Gene1")
soft_margin_costs.write("\t")
soft_margin_costs.write("Gene2")
soft_margin_costs.write("\t")
soft_margin_costs.write("Best_Cost")
soft_margin_costs.write("\t")
soft_margin_costs.write("Accuracy_CV_LOO")
soft_margin_costs.write("\t")
soft_margin_costs.write("StDev_Accuracy")
soft_margin_costs.write("\t")
soft_margin_costs.write("Margin")
soft_margin_costs.write("\t")
soft_margin_costs.write("Number_of_Support_Vectors")
soft_margin_costs.write("\t")
soft_margin_costs.write("w1")
soft_margin_costs.write("\t")
soft_margin_costs.write("w2")
soft_margin_costs.write("\t")
soft_margin_costs.write("b")
soft_margin_costs.write("\n")

#Para cada par de genes LS, ver o custo com melhor score no CV, em caso de empate
# usar o custo com maior margem
counter = 0
for line in lines_lp:
	counter += 1
	line = line.rstrip()
	line = line.split("\t")
	#Dicionarios
	costs_results = {}
	margin_results = {}
	std_dev_support_vectors_weights = {}
	best_scores = list()
	best_cost = None
	gene1_index = cancerData[cancerData["Genes"] == line[0]].index.tolist()
	gene2_index = cancerData[cancerData["Genes"] == line[1]].index.tolist()
	braca1 = np.zeros(shape=(len(X[gene1_index,BRACA1]),2))
	braca2 = np.zeros(shape=(len(X[gene1_index,BRACA2]),2))
	l = 0
	q = 0
	#Montando os vetores gene1,gene2 para cada amostra
	for z in range(0,len(X[gene1_index,BRACA1])):
		features_braca1 = np.array([X[gene1_index,BRACA1].astype(np.float)[z],X[gene2_index,BRACA1].astype(np.float)[z]])
		braca1[l] = features_braca1
		l += 1 
	for k in range(0,len(X[gene1_index,BRACA2])):
		features_braca2 = np.array([X[gene1_index,BRACA2].astype(np.float)[k],X[gene2_index,BRACA2].astype(np.float)[k]])
		braca2[q] = features_braca2
		q += 1 
	#Montando os vetores de caracteristicas e classes
	features = np.vstack([braca1,braca2])
	labels = np.concatenate([labels_braca1,labels_braca2])
	for c in costs:
		print "Evaluating LS gene pair: " + str(counter)
		print "Using cost: " + str(c)
		#Rodando o SVM linear e o cross-validation leave-one-out
		kf = cross_validation.KFold(15, n_folds=15)
		classifier = svm.SVC(C=c, kernel = "linear")
		#Pegando os scores do cv 
		scores = cross_validation.cross_val_score(classifier, features, labels, cv=kf)
		print("Accuracy: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))
		#Calculando a margem para esse score
		classifier.fit(features, labels)
		vector_w = classifier.coef_.tolist()
		w1 = vector_w[0][0]
		w2 = vector_w[0][1]
		b = classifier.intercept_[0]
		Margin = 1.0 / math.sqrt(w1 ** 2 + w2 ** 2)
		#Pegando o numero de vetores de supporte
		support_vectors = len(classifier.support_vectors_)
		#Inserindo os valores de std dev do score, numero de vetores de suporte, w1, w2 e b no dicionario
		if std_dev_support_vectors_weights.has_key(str(c)):
			std_dev_support_vectors_weights[str(c)] = (float(scores.std()),support_vectors,w1,w2,b)
		else:
			std_dev_support_vectors_weights[str(c)] = (float(scores.std()),support_vectors,w1,w2,b)
		#Inserindo os scores no dicionario
		if costs_results.has_key(str(c)):
			costs_results[str(c)] = float(scores.mean())
		else:
			costs_results[str(c)] = float(scores.mean())
		#Inserindo as margens no dicionario
		if margin_results.has_key(str(c)):
			margin_results[str(c)] = Margin
		else:
			margin_results[str(c)] = Margin
	#Ordenando os custos pelo maior score
	sorted_costs_results = sorted(costs_results.items(), key=operator.itemgetter(1), reverse=True)
	#Vendo se ha empate dos scores entre os custos
	if sorted_costs_results[0][1] == sorted_costs_results[1][1]:
		best_scores.append(sorted_costs_results[0][0])
		best_scores.append(sorted_costs_results[1][0])
	else:
		best_cost = sorted_costs_results[0][0] # se nao ha empate	
		if sorted_costs_results[1][1] == sorted_costs_results[2][1]: # se ha varios empates dos scores
			best_scores.append(sorted_costs_results[2][0])
			if sorted_costs_results[2][1] == sorted_costs_results[3][1]:
				best_scores.append(sorted_costs_results[3][0])
				if sorted_costs_results[3][1] == sorted_costs_results[4][1]:
					best_scores.append(sorted_costs_results[4][0])
	print sorted_costs_results
	print margin_results
	#Se ha empate nos scores do custo, pegar o custo com maior margem
	if best_cost is None:
		margins_top_scores = {your_key:margin_results[your_key] for your_key in best_scores}
		best_cost = max(margins_top_scores.iteritems(), key=operator.itemgetter(1))[0]
		print best_scores
		print margins_top_scores
	print "BEST COST: " + best_cost
	#Adicionando os resultados ao arquivo
	soft_margin_costs.write(line[0])
	soft_margin_costs.write("\t")
	soft_margin_costs.write(line[1])
	soft_margin_costs.write("\t")
	soft_margin_costs.write(best_cost)
	soft_margin_costs.write("\t")
	soft_margin_costs.write(str(costs_results[best_cost]))
	soft_margin_costs.write("\t")
	soft_margin_costs.write(str(std_dev_support_vectors_weights[best_cost][0]))
	soft_margin_costs.write("\t")
	soft_margin_costs.write(str(margin_results[best_cost]))
	soft_margin_costs.write("\t")
	soft_margin_costs.write(str(std_dev_support_vectors_weights[best_cost][1]))
	soft_margin_costs.write("\t")
	soft_margin_costs.write(str(std_dev_support_vectors_weights[best_cost][2]))
	soft_margin_costs.write("\t")
	soft_margin_costs.write(str(std_dev_support_vectors_weights[best_cost][3]))
	soft_margin_costs.write("\t")
	soft_margin_costs.write(str(std_dev_support_vectors_weights[best_cost][4]))
	soft_margin_costs.write("\n")

ls_genes.close()
soft_margin_costs.close()

#Ordenando os pares pelo tamanho da margem
df = pd.read_csv("soft_margin_costs.txt", sep ="\t")
df = df.sort('Margin', ascending = False)
df_100 = df[0:100]
df_10 = df[0:10]
df_100.to_csv('Top_100_Soft_Margin_Costs.txt', index=False, sep= "\t")

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
	bias = float(row[9])
	w1 = float(row[7])
	w2 = float(row[8])
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
	plt.title("Soft Linear Margin Scatter Plot")
	plt.xlabel(row[0])
	plt.ylabel(row[1])
	plt.savefig('Soft_Margin_Pair_' + str(k) + '.png')

