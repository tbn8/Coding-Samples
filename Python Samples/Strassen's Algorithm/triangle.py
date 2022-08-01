import sys
import numpy as np
import time
import random
from scipy.special import comb
import matplotlib.pyplot as plt
from strassen import *

def makeGraph(n,p):
    #create random adjacency matrix
    graph=np.random.choice([0,1], size=(n,n), p=[1-p,p])
    #all diagonal entries should be 0 b/c there's no edge from a vertex to itself
    for i in range(n):
        graph[i,i]=0
    #make lower triangular matrix equal upper triangular
    i_lower = np.tril_indices(n, -1)
    graph[i_lower]=graph.T[i_lower]
    return graph

def countTriangles(A):
    n=A.shape[0]
    A_square=strassen(A,A,n,4)
    A_cube=strassen(A,A_square,n,4)
    return np.sum(np.diag(A_cube))/6

if __name__ == '__main__':
    # start1=time.time()
    n=1024
    pvalues=[0.01,0.02,0.03,0.04,0.05]
    num=[]
    expected=[]

    for i in range(len(pvalues)):
        p=pvalues[i]
        graph=makeGraph(n,p)
        num.append(countTriangles(graph))
        expected.append(comb(n,3)*(p**3))

    num=np.array(num)
    expected=np.array(expected)
    print(num,"\n",expected)
    # colors=np.array([])
    plt.scatter(num,expected,color='darkturquoise')

    for i,p in enumerate(pvalues):
        plt.annotate(p,(num[i],expected[i]),xytext=(num[i],expected[i]+200), ha='center',va='bottom',size=7)
    plt.title("Number of triangles in random graphs")
    plt.xlabel("Generated")
    plt.ylabel("Expected")
    # plt.show()
    plt.savefig('triangle.pdf')
