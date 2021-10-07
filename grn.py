import pandas as pd
import numpy as np
import csv
import os
import sys
import scipy as sp
import scipy.sparse as sprs
import scipy.spatial
import scipy.sparse.linalg
from scipy.sparse import csc_matrix
import igraph

path=os.getcwd() # Working directory
os.chdir(path) # Moving to working directory
file_name1 = sys.argv[-1]

# READ FILE TWO COLUMNS FUNCTION
def Read_Two_Column_File(file_name):
    with open(file_name, 'r') as f_input:
        csv_input = csv.reader(f_input, delimiter='\t', skipinitialspace=True)
        x = []
        y = []
        for cols in csv_input:
            x.append(cols[0])
            y.append(cols[1])

    return x, y

## PAGE RANK FUNCTION
# Reference: https://asajadi.github.io/fast-pagerank/
def pageRank(G, s = .85, maxerr = .001):
    '''
    Computes the pagerank for each of the n states.
    Used in webpage ranking and text summarization using unweighted
    or weighted transitions respectively.
    Args
    ----------
    G: matrix representing state transitions
       Gij can be a boolean or non negative real number representing the
       transition weight from state i to j.
    Kwargs
    ----------
    s: probability of following a transition. 1-s probability of teleporting
       to another state. Defaults to 0.85
    maxerr: if the sum of pageranks between iterations is bellow this we will
            have converged. Defaults to 0.001
    '''
    n = G.shape[0]

    # transform G into markov matrix M
    #M = csc_matrix(G,dtype=np.float)
    M = csc_matrix(G,dtype=float)
    rsums = np.array(M.sum(1))[:,0]
    ri, ci = M.nonzero()
    M.data /= rsums[ri]

    # bool array of sink states
    sink = rsums==0

    # Compute pagerank r until we converge
    ro, r = np.zeros(n), np.ones(n)
    while np.sum(np.abs(r-ro)) > maxerr:
        ro = r.copy()
        # calculate each pagerank at a time
        #for i in xrange(0,n):
        for i in range(0,n):
            # inlinks of state i
            Ii = np.array(M[:,i].todense())[:,0]
            # account for sink states
            Si = sink / float(n)
            # account for teleportation to state i
            Ti = np.ones(n) / float(n)

            r[i] = ro.dot( Ii*s + Si*s + Ti*(1-s) )

    # return normalized pagerank
    return r/sum(r)

def pagerank_power(A, p=0.85, max_iter=100,
                   tol=1e-06, personalize=None, reverse=False):
    ''' Calculates PageRank given a csr graph

    Inputs:
    -------
    A: a csr graph.
    p: damping factor
    max_iter: maximum number of iterations
    personlize: if not None, should be an array with the size of the nodes
                containing probability distributions.
                It will be normalized automatically.
    reverse: If true, returns the reversed-PageRank

    Returns:
    --------
    PageRank Scores for the nodes

    '''
    # In Moler's algorithm, $$G_{ij}$$ represents the existences of an edge
    # from node $$j$$ to $$i$$, while we have assumed the opposite!
    if reverse:
        A = A.T

    n, _ = A.shape
    r = np.asarray(A.sum(axis=1)).reshape(-1)

    k = r.nonzero()[0]

    D_1 = sprs.csr_matrix((1 / r[k], (k, k)), shape=(n, n))

    if personalize is None:
        personalize = np.ones(n)
    personalize = personalize.reshape(n, 1)
    s = (personalize / personalize.sum()) * n

    z_T = (((1 - p) * (r != 0) + (r == 0)) / n)[sp.newaxis, :]
    W = p * A.T @ D_1

    x = s
    oldx = np.zeros((n, 1))

    iteration = 0

    while sp.linalg.norm(x - oldx) > tol:
        oldx = x
        x = W @ x + s @ (z_T @ x)
        iteration += 1
        if iteration >= max_iter:
            break
    x = x / sum(x)

    return x.reshape(-1)

# Discretize columns using previus bin data:
def stDiscretizer(serie):
	mean = serie.mean()
	sd = serie.std()

	column = list(serie)

	#'F'
	e = mean + (2 * sd)
	#'E'
	d = mean + sd
	#'D'
	c = mean
	#'C'
	b = mean - sd
	#'B'
	a = mean - (2 * sd)
	#'A'
	
	result = []

	for x in column:
		if x < a:
			result.append('A')
		elif (x >= a) & (x < b):
			result.append('B')
		elif (x >= b) & (x < c):
			result.append('C')
		elif (x >= c) & (x < d):
			result.append('D')
		elif (x >= d) & (x < e):
			result.append('E')
		elif (x >= e):
			result.append('F')
		elif np. isnan(x):
			result.append('?')
		else:		
			print('Unespected error', str(mean), str(sd) , str(x))

	return result

def qDiscretizer(df):
	df_res = df.copy()
	for i in df_res.columns[0:len(df_res.columns)]:
		#print(i, 'Discretizing')
		#print('Normalizing ', i)
		df_res[i] = stDiscretizer(df_res[i])
	return df_res

# Reading input file and cocatenating the two columns
x, y = Read_Two_Column_File(file_name1)
c1 = x+y
c2 = y+x

# Generating Adjacency Matrix
df1 = pd.DataFrame({'Gene_A': c1,'Gene_B': c2,})
df1 = pd.crosstab(df1.Gene_A, df1.Gene_B)
idx = df1.columns.union(df1.index)
df1 = df1.reindex(index = idx, columns=idx, fill_value=0)
df1.to_excel('01_matrix_adjacences.xlsx')

# Computing Page_Rank
# Reference: https://asajadi.github.io/fast-pagerank/
G=df1.to_numpy()
PG = pagerank_power(G)
#df2 = pd.DataFrame(PG)
#df2.to_excel('page_rank.xlsx')

# Computing Degree
DG = G.sum(axis=1)
#df3 = pd.DataFrame(DG)
#df3.to_excel('degree.xlsx')

# Computing KNN
# Reference: http://olizardo.bol.ucla.edu/classes/soc-111/textbook/_book/5-7-sec-knn.html#sec:knn
n = G.shape[0]
g = igraph.Graph.Adjacency(G, directed=False) #Generating graph from adjacency matrix
g1=igraph.Graph.simplify(g)
# Method 1
#for i in range(0,n):
#    soma = 0
#    for j in range(0,n):
#        if G[i][j] > 0:
#            soma=soma+DG[j]
#    KNN[i]=soma/DG[i]
# Method 2
aux1=igraph.Graph.knn(g1)
aux2=np.asarray(aux1, dtype=object)
aux3=pd.DataFrame(aux2[0])
KNN=aux3.div(2)
#df4 = pd.DataFrame(KNN)
#df4.to_excel('knn.xlsx')

# Creating a file with topological parameters
a1=pd.DataFrame(KNN)
a2=pd.DataFrame(PG)
a3=pd.DataFrame(DG)
df1['GENE'] = df1.index #Create a column with df index
index = df1.iloc[:,-1:] #Getting the last column from df as dataframe
index = index.reset_index(drop=True)
df5 = pd.concat([index,a1,a2,a3], axis=1)
df5 = df5.set_index('GENE') #Tranform GENE column in to index
df5.columns.values[0] = 'KNN' #Renaming colmuns names
df5.columns.values[1] = 'PAGE_RANK'
df5.columns.values[2] = 'DEGREE'
df5.to_excel('02_network_parameters.xlsx')

#Discretizing data
df6 = qDiscretizer(df5)
#df6.to_excel('03_network_parameters_discretized.xlsx')
df6.to_csv('03_network_parameters_discretized.csv')

# Generating ARFF
CLASS = ["?" for x in range(n)]
df7 = pd.concat([index,pd.DataFrame(CLASS)], axis=1)
df7 = df7.set_index('GENE') #Tranform GENE column in to index
df7 = pd.concat([df6,df7], axis=1)
df7.columns.values[3] = 'CLASS'
dfaux = df7.iloc[:, [2, 1, 0, 3]]
df8=dfaux.values.tolist()

file_name='04_file_for_classification.arff'
file = open(file_name, 'w+', newline ='') 
with file:
    write = csv.writer(file) 
    write.writerows(df8)

def remove_empty_lines(filename):
    """Overwrite the file, removing empty lines and lines that contain only whitespace."""
    with open(filename, 'r+') as f:
        lines = f.readlines()
        f.seek(0)
        f.writelines(line for line in lines if line.strip())
        f.truncate()

remove_empty_lines(file_name)

def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)

txt1='@relation 03_Atributos_discretizados_Scerevisiae_RegulatoryNet-weka.filters.unsupervised.attribute.Remove-R1'
txt2='\n'
txt3='@attribute degree {D,C,F,E}'
txt4='@attribute page_rank {D,C,F,E}'
txt5='@attribute Average_nearest_neighbor_degree {B,D,C,E,F,A}'
txt6='@attribute class {regulators,targets}'
txt7='@data'

line_prepender(file_name,txt7)
line_prepender(file_name,txt2)
line_prepender(file_name,txt6)
line_prepender(file_name,txt5)
line_prepender(file_name,txt4)
line_prepender(file_name,txt3)
line_prepender(file_name,txt2)
line_prepender(file_name,txt1)

