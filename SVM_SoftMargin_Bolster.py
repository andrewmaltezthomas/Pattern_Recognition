from sklearn import cross_validation
from sklearn import svm
from sklearn import preprocessing
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
soft_margin_costs = open("bolstered_error_costs.txt", "w")
soft_margin_costs.write("Gene1")
soft_margin_costs.write("\t")
soft_margin_costs.write("Gene2")
soft_margin_costs.write("\t")
soft_margin_costs.write("Best_Cost")
soft_margin_costs.write("\t")
soft_margin_costs.write("Accuracy")
soft_margin_costs.write("\t")
soft_margin_costs.write("Bolstered_Error")
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

#Para cada par de genes LS, ver o custo com menor Bolstered Error, em caso de empate
# usar o custo com maior margem
counter = 0
for line in lines_lp:
	counter += 1
	line = line.rstrip()
	line = line.split("\t")
	#Dicionarios
	costs_results = {}
	margin_results = {}
	support_vectors_weights = {}
	mean_bolstered_errors = {}
	least_bolstered_errors = []
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
	features = np.concatenate([braca1,braca2])
	labels = np.concatenate([labels_braca1,labels_braca2])
	
	# Rodando as analises para cada custo
	for c in costs:
		print "Evaluating LS gene pair: " + str(counter)
		print "Using cost: " + str(c)
		bolstered_errors = []
		
		#Rodando o SVM linear e pegando a margem
		classifier = svm.SVC(C=c, kernel = "linear")
		classifier.fit(features, labels)
		vector_w = classifier.coef_.tolist()
		w1 = vector_w[0][0]
		w2 = vector_w[0][1]
		b = classifier.intercept_[0]
		Margin = 1.0 / np.sqrt(np.sum(classifier.coef_ ** 2))
		support_vectors = len(classifier.support_vectors_)
		
		#Calculando o erro de bolster
		predictions = classifier.predict(features)
		accuracy = 0.0 
		error = np.zeros(15)
		z = np.zeros(15)
		for i in range(0,15):
			z[i] = np.float(((np.matrix(vector_w[0]) * np.matrix(features[i]).transpose()) + b) / math.sqrt(w1 ** 2 + w2 ** 2))
			if predictions[i] == labels[i]:
				error[i]= 0.5 - math.sqrt(2) * math.erf(z[i]/4*0.6)
				accuracy += 1.0
			else:
				error[i]=0.5-2*math.erf(z[i]/4*0.6) 
		mean_bolster_error = sum(error)/len(error)
		accuracy = (accuracy / len(labels) * 100)
		
		#Inserindo os valores nos dicionarios
		if mean_bolstered_errors.has_key(str(c)):
			mean_bolstered_errors[str(c)] = (mean_bolster_error)
		else:
			mean_bolstered_errors[str(c)] = (mean_bolster_error)
		if support_vectors_weights.has_key(str(c)):
			support_vectors_weights[str(c)] = (support_vectors,w1,w2,b,accuracy)
		else:
			support_vectors_weights[str(c)] = (support_vectors,w1,w2,b,accuracy)
		
		#Inserindo as margens no dicionario
		if margin_results.has_key(str(c)):
			margin_results[str(c)] = Margin
		else:
			margin_results[str(c)] = Margin
	
	#Ordenando os custos pelo maior score
	sorted_mean_bolstered_errors = sorted(mean_bolstered_errors.items(), key=operator.itemgetter(1))
	
	#Vendo se ha empate dos bolstered errors entre os custos
	if sorted_mean_bolstered_errors[0][1] == sorted_mean_bolstered_errors[1][1]:
		least_bolstered_errors.append(sorted_mean_bolstered_errors[0][0])
		least_bolstered_errors.append(sorted_mean_bolstered_errors[1][0])
	else:
		best_cost = sorted_mean_bolstered_errors[0][0] # se nao ha empate	
		if sorted_mean_bolstered_errors[1][1] == sorted_mean_bolstered_errors[2][1]: # se ha varios empates dos scores
			least_bolstered_errors.append(sorted_mean_bolstered_errors[2][0])
			if sorted_mean_bolstered_errors[2][1] == sorted_mean_bolstered_errors[3][1]:
				least_bolstered_errors.append(sorted_mean_bolstered_errors[3][0])
				if sorted_mean_bolstered_errors[3][1] == sorted_mean_bolstered_errors[4][1]:
					least_bolstered_errors.append(sorted_mean_bolstered_errors[4][0])
	print sorted_mean_bolstered_errors
	
	#Se ha empate nos bolstered errors do custo, pegar o custo com maior margem
	if best_cost is None:
		margins_top_scores = {your_key:margin_results[your_key] for your_key in least_bolstered_errors}
		best_cost = max(margins_top_scores.iteritems(), key=operator.itemgetter(1))[0]
		print least_bolstered_errors
		print margins_top_scores
	print "BEST COST: " + best_cost
	
	#Adicionando os resultados ao arquivo
	soft_margin_costs.write(line[0])
	soft_margin_costs.write("\t")
	soft_margin_costs.write(line[1])
	soft_margin_costs.write("\t")
	soft_margin_costs.write(best_cost)
	soft_margin_costs.write("\t")
	soft_margin_costs.write(str(support_vectors_weights[best_cost][4]))
	soft_margin_costs.write("\t")
	soft_margin_costs.write(str(mean_bolstered_errors[best_cost]))
	soft_margin_costs.write("\t")
	soft_margin_costs.write(str(margin_results[best_cost]))
	soft_margin_costs.write("\t")
	soft_margin_costs.write(str(support_vectors_weights[best_cost][0]))
	soft_margin_costs.write("\t")
	soft_margin_costs.write(str(support_vectors_weights[best_cost][1]))
	soft_margin_costs.write("\t")
	soft_margin_costs.write(str(support_vectors_weights[best_cost][2]))
	soft_margin_costs.write("\t")
	soft_margin_costs.write(str(support_vectors_weights[best_cost][3]))
	soft_margin_costs.write("\n")
ls_genes.close()
soft_margin_costs.close()

#Ordenando os pares pelo tamanho da margem
df = pd.read_csv("bolstered_error_costs.txt", sep ="\t")
df = df.sort('Margin', ascending = False)
df_100 = df[0:100]
df_10 = df[0:10]
df_100.to_csv('Top_100_Maximum_Margin_Bolstered_Error_Costs.txt', index=False, sep= "\t")

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
	print X[i,BRACA1]
	print X[j,BRACA1]
	print X[i,BRACA2]
	print X[j,BRACA2]
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
	plt.savefig('Margin_Bolstered_Error_Pair_' + str(k) + '.png')

df = pd.read_csv("bolstered_error_costs.txt", sep ="\t")
df = df.sort('Bolstered_Error')
df_100 = df[0:100]
df_10 = df[0:10]
df_100.to_csv('Top_100_Minimum_Bolstered_Error_Costs.txt', index=False, sep= "\t")

#Plotar os resultados para os 10 pares com menor bolstered error 
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
	print X[i,BRACA1]
        print X[j,BRACA1]
        print X[i,BRACA2]
        print X[j,BRACA2]
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
	plt.savefig('Minimum_Bolstered_Error_Pair_' + str(k) + '.png')
