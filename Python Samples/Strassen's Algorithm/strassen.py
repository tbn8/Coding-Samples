import sys
import numpy as np
import time
import random

def split(matrix,n):
    half=n//2
    return matrix[:half, :half], matrix[:half, half:], matrix[half:, :half], matrix[half:, half:]

def conventionalMultiply(x,y,n):
    x=x.tolist()
    y=y.tolist()
    z=[[0]*n for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                z[i][j]+=x[i][k]*y[k][j]
    z=np.array(z)
    return z


def variantStrassen(x,y,n,crossover):
    if n<=crossover:
        result=conventionalMultiply(x,y,n)
    else:
        a,b,c,d=split(x,n)
        e,f,g,h=split(y,n)
        half=n//2
        p1=variantStrassen(a, f-h, half, crossover)
        p2=variantStrassen(a+b, h, half, crossover)
        p3=variantStrassen(c+d ,e, half, crossover)
        p4=variantStrassen(d, g-e, half, crossover)
        p5=variantStrassen(a+d, e+h, half, crossover)
        p6=variantStrassen(b-d, g+h, half, crossover)
        p7=variantStrassen(c-a, e+f, half, crossover)

        q1=-p2+p4+p5+p6
        q2=p1+p2
        q3=p3+p4
        q4=p1-p3+p5+p7

        result=np.vstack((np.hstack((q1,q2)),np.hstack((q3,q4))))
    return result

def fastStrassen(x,y,n,crossover):
    if n<=crossover:
        result=conventionalMultiply(x,y,n)
    else:
        result=np.zeros(shape=(n,n), dtype=int) #final matrix product
        half=n//2

        #compute F-H and store in the F quadrant of Y
        y[:half,half:]-=y[half:,half:]
        #compute P1=A(F-H) and store in temp matrix
        temp=fastStrassen(x[:half,:half],y[:half,half:],half,crossover)
        #add P1 to Q2 and Q4 of final matrix
        result[:half,half:]+=temp
        result[half:,half:]+=temp
        #reset F
        y[:half,half:]+=y[half:,half:]

        #compute A+B and store in A quadrant of X
        x[:half,:half]+=x[:half,half:]
        #compute P2=(A+B)H and store in temp matrix
        temp=fastStrassen(x[:half,:half],y[half:,half:],half,crossover)
        #add P2 to Q1 and Q2 of final matrix
        result[:half,:half]-=temp
        result[:half,half:]+=temp
        #reset A
        x[:half,:half]-=x[:half,half:]

        #compute C+D and store in C quadrant of X
        x[half:,:half]+=x[half:,half:]
        #compute P3=(C+D)E and store in temp matrix
        temp=fastStrassen(x[half:,:half],y[:half,:half],half,crossover)
        #add P3 to Q3 and Q4 of final matrix
        result[half:,:half]+=temp
        result[half:,half:]-=temp
        #reset C
        x[half:,:half]-=x[half:,half:]

        #compute G-E and store in G quadrant of Y
        y[half:,:half]-=y[:half,:half]
        #compute P4=D(G-E) and store in temp matrix
        temp=fastStrassen(x[half:,half:],y[half:,:half],half,crossover)
        #add P4 to Q1 and Q3 of final matrix
        result[:half,:half]+=temp
        result[half:,:half]+=temp
        #reset G
        y[half:,:half]+=y[:half,:half]

        #compute A+D and store in A
        x[:half,:half]+=x[half:,half:]
        #compute E+H and store in E
        y[:half,:half]+=y[half:,half:]
        #compute P5
        temp=fastStrassen(x[:half,:half],y[:half,:half],half,crossover)
        #add P5 to Q1 and Q4
        result[:half,:half]+=temp
        result[half:,half:]+=temp
        #reset A and E
        x[:half,:half]-=x[half:,half:]
        y[:half,:half]-=y[half:,half:]

        #compute B-D and store in B
        x[:half,half:]-=x[half:,half:]
        #compute G+H and store in G
        y[half:,:half]+=y[half:,half:]
        #compute P6
        temp=fastStrassen(x[:half,half:],y[half:,:half],half,crossover)
        #add P6 to Q1
        result[:half,:half]+=temp
        #reset B and G
        x[:half,half:]+=x[half:,half:]
        y[half:,:half]-=y[half:,half:]

        #compute C-A and store in C
        x[half:,:half]-=x[:half,:half]
        #compute E+F and store in E
        y[:half,:half]+=y[:half,half:]
        #compute P7
        temp=fastStrassen(x[half:,:half],y[:half,:half],half,crossover)
        #add P7 to Q4
        result[half:,half:]+=temp
        #reset C and E
        x[half:,:half]+=x[:half,:half]
        y[:half,:half]-=y[:half,half:]

    return result

def pad(matrix,n,new_n):
    new_matrix=np.zeros(shape=(new_n,new_n), dtype=int)
    new_matrix[:n,:n]=matrix
    return new_matrix

def strassen(x,y,n,crossover):
    #if n is not a power of 2
    if (n & (n-1) != 0):
        n_0=2;
        while (n_0<n):
            n_0*=2;
        x_0=pad(x,n,n_0)
        y_0=pad(y,n,n_0)
        product=variantStrassen(x_0,y_0,n_0,crossover)
        # product=fastStrassen(x_0,y_0,n_0,crossover)
        result=product[:n,:n]
    else:
        result=variantStrassen(x,y,n,crossover)
        # result=fastStrassen(x,y,n,crossover)

    return result

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Invalid arguments")
    else:
        n=int(sys.argv[2])

        #for submission
        if (sys.argv[1]=="0"):
            file=sys.argv[3]
            arr=np.loadtxt(file, dtype='int')

            X=np.reshape(arr[:n**2],(n,n))
            Y=np.reshape(arr[n**2:],(n,n))

        #for testing
        else:
            np.random.seed(0)
            X=np.random.randint(low=0,high=3,size=(n,n))
            Y=np.random.randint(low=0,high=3,size=(n,n))
            # print(X)
            # print(Y)
        # start=time.time()
        result=strassen(X,Y,n,4)
        # print("--- %s seconds ---" % (time.time() - start))

        # output diagonal entries
        output=np.diag(result)
        print(*output, sep='\n')
