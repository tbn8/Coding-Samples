import sys
import datetime
import random
import math
import pandas as pd
from heap import *
from heapq import heappop, heappush, heapify

def KarmarkarKarp(list):
    #self-implementation of heap (very slow)
    if sys.argv[1]=="1":
        A=MaxHeap(INPUT_N)
        for i in range(INPUT_N):
            A.insert(list[i])
        a1=A.extractMax()
        a2=A.extractMax()
        while (a2!=0):
            A.insert(a1-a2)
            A.insert(0)
            a1=A.extractMax()
            a2=A.extractMax()
        return a1

    #priority queue using Python native list (faster, used for grading)
    if sys.argv[1]=="2" or sys.argv[1]=="0":
        A=list[:]
        for _ in range(INPUT_N-1):
            A.sort(reverse=True)
            diff=A[0]-A[1]
            A.remove(A[0])
            A.remove(A[0])
            A.append(diff)
        return A[0]

    #Python heapq package (fastest, not allowed for submission)
    if sys.argv[1]=="3":
        A=[]
        heapify(A)
        for i in range(INPUT_N):
            heappush(A,-1*list[i])
        a1=heappop(A)
        a2=heappop(A)
        while (a2!=0):
            heappush(A,a1-a2)
            heappush(A,0)
            a1=heappop(A)
            a2=heappop(A)
        return -a1

#generate random solution
def randSol():
    solution = []
    for i in range(INPUT_N):
        num = random.choice((-1,1))
        solution.append(num)
    return solution

#generate a random neighbor solution
def Neighbor(S):
    neighbor=list(S)
    i=random.randint(0, INPUT_N-1)
    neighbor[i]=(-1)*S[i]

    coin=random.choice((0,1))
    if (coin==0):
        return neighbor
    if (coin==1):
        j=random.randint(0, INPUT_N-1)
        r=[*range(0,i),*range(i+1,INPUT_N)] #create list of indexes excluding i
        j=random.choice(r)
        neighbor[j]=(-1)*S[j]
        return neighbor

def prePartition(A):
    P=[]
    A_prime=[0]*INPUT_N
    for i in range(INPUT_N):
        num=random.randint(0, INPUT_N-1)
        P.append(num)
    for i in range(INPUT_N):
        A_prime[P[i]]+=A[i]
    return A_prime, P

#generate neighbor solution for prepartioning
def preNeighbor(A,P):
    new_P=list(P)
    neighbor=[0]*INPUT_N
    i=random.randint(0, INPUT_N-1)
    r=[*range(0,P[i]),*range(P[i]+1,INPUT_N)] #create list of indexes excluding p_i
    j=random.choice(r)
    assert (P[i]!=j)
    new_P[i]=j
    for k in range(INPUT_N):
        neighbor[new_P[k]]+=A[k]
    return neighbor, new_P

#calculate residue
def residue(A, S):
    r=0
    for i in range(INPUT_N):
        r+=A[i]*S[i]
    return int(math.fabs(r))

# Repeated Random Algorithm for Standard Representation
def repRandom(A):
    S=randSol() #starting solution
    for i in range(MAX_ITER):
        S_prime=randSol() #random solution
        if (residue(A, S_prime) < residue(A, S)):
            S=list(S_prime)
    return residue(A,S)

# Repeated Random Algorithm for Prepartitioning
def preRepRandom(A):
    A_prime,P=prePartition(A) #starting solution
    for i in range(MAX_ITER):
        new_A_prime,new_P=prePartition(A) #random solution
        if (KarmarkarKarp(new_A_prime) < KarmarkarKarp(A_prime)):
            A_prime=list(new_A_prime)
    return KarmarkarKarp(A_prime)

### Hill Climbing Algorithm for Standard Representation
def hillClimb(A):
    S=randSol() #starting solution
    for i in range(MAX_ITER):
        S_prime=Neighbor(S) #random neighbor solution
        if (residue(A, S_prime) < residue(A, S)):
            S=list(S_prime)
    return residue(A,S)

### Hill Climbing Algorithm for Prepartitioning
def preHillClimb(A):
    A_prime,P=prePartition(A) #starting solution
    for i in range(MAX_ITER):
        new_A_prime,new_P=preNeighbor(A,P) #random neighbor solution
        if (KarmarkarKarp(new_A_prime) < KarmarkarKarp(A_prime)):
            A_prime=list(new_A_prime)
            P=list(new_P)
    return KarmarkarKarp(A_prime)

#cooling schedule
def T(iter):
    return pow(10,10)*pow(0.8, (iter/300))

### Simulated Annealing Algorithm for Standard Representation
def simAnneal(A):
    S=randSol() #starting solution
    S_prime_2=S
    for i in range(MAX_ITER):
        S_prime=Neighbor(S) #random neighbor solution
        if (residue(A, S_prime) < residue(A, S)):
            S=list(S_prime)
        else:
            prob=math.exp(-(residue(A,S_prime)-residue(A,S))/T(i))
            if (random.random() <= prob):
                S=list(S_prime)
        if (residue(A,S) < residue(A,S_prime_2)):
            S_prime_2=list(S)
    return residue(A,S_prime_2)

### Simulated Annealing Algorithm for Prepartitioning
def preSimAnneal(orig_A):
    A,P=prePartition(orig_A) #starting solution
    A_prime=list(A)
    P_prime=list(P)
    for i in range(MAX_ITER):
        new_A,new_P=preNeighbor(orig_A,P) #random neighbor solution
        if (KarmarkarKarp(new_A) < KarmarkarKarp(A)):
            A=list(new_A)
            P=list(new_P)
        else:
            prob=math.exp(-(KarmarkarKarp(new_A)-KarmarkarKarp(A))/T(i))
            if (random.random() <= prob):
                A=list(new_A)
                P=list(new_P)
        if (KarmarkarKarp(A) < KarmarkarKarp(A_prime)):
            A_prime=list(A)
            P_prime=list(P)
    return KarmarkarKarp(A_prime)

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Invalid arguments")
    else:
        # for submission
        if (sys.argv[1]=="0"):
            file=open(sys.argv[3],'r')
            A=[int(x) for x in file.readlines()]
            INPUT_N=len(A)
            MAX_ITER=25000
            if (sys.argv[2]=="0"):
                result=KarmarkarKarp(A)
            if (sys.argv[2]=="1"):
                result=repRandom(A)
            if (sys.argv[2]=="2"):
                result=hillClimb(A)
            if (sys.argv[2]=="3"):
                result=simAnneal(A)
            if (sys.argv[2]=="11"):
                result=preRepRandom(A)
            if (sys.argv[2]=="12"):
                result=preHillClimb(A)
            if (sys.argv[2]=="13"):
                result=preSimAnneal(A)
            print(result)

        #for testing
        else:
            # random.seed(8)
            MAX_ITER=25000
            NUM_INS=5
            output=[]

            for i in range(NUM_INS):
                A=[]
                for _ in range(100):
                    A.append(random.randint(1,10**12))
                INPUT_N=len(A)

                instance=i+1

                start=datetime.datetime.now()
                result=KarmarkarKarp(A)
                runtime=(datetime.datetime.now()-start).total_seconds()*1000
                output.append([instance,"KarmarkarKarp",result,runtime])

                start=datetime.datetime.now()
                result=repRandom(A)
                runtime=(datetime.datetime.now()-start).total_seconds()*1000
                output.append([instance,"Repeated Random",result,runtime])

                start=datetime.datetime.now()
                result=hillClimb(A)
                runtime=(datetime.datetime.now()-start).total_seconds()*1000
                output.append([instance,"Hill Climbing",result,runtime])

                start=datetime.datetime.now()
                result=simAnneal(A)
                runtime=(datetime.datetime.now()-start).total_seconds()*1000
                output.append([instance,"Simulated Annealing",result,runtime])

                start=datetime.datetime.now()
                result=preRepRandom(A)
                runtime=(datetime.datetime.now()-start).total_seconds()*1000
                output.append([instance,"Prepartitioned Repeated Random",result,runtime])

                start=datetime.datetime.now()
                result=preHillClimb(A)
                runtime=(datetime.datetime.now()-start).total_seconds()*1000
                output.append([instance,"Prepartitioned Hill Climbing",result,runtime])

                start=datetime.datetime.now()
                result=preSimAnneal(A)
                runtime=(datetime.datetime.now()-start).total_seconds()*1000
                output.append([instance,"Prepartitioned Simulated Annealing",result,runtime])

            df=pd.DataFrame(output, columns=['Instance','Algorithm','Residue','Runtime'])
            sum=df.iloc[:,1:].groupby('Algorithm', sort=False).describe()

            if sys.argv[1]=="1":
                df.to_excel('Partition_myheap.xlsx')
                sum.to_excel('Partition Summary_myheap.xlsx')

            if sys.argv[1]=="2":
                df.to_excel('Partition_list.xlsx')
                sum.to_excel('Partition Summary_list.xlsx')

            if sys.argv[1]=="3":
                df.to_excel('Partition_heapq.xlsx')
                sum.to_excel('Partition Summary_heapq.xlsx')
