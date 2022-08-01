### Test maximum likelihood methods ###

import pandas as pd
import numpy as np
import numpy.random as nr
import statsmodels.api as sm
from scipy import stats
import seaborn as sns
from math import log, sqrt, sin
import matplotlib.pyplot as plt

#read data
housing = pd.read_csv('../CSCI-E-83/data/housing.csv')
housing = housing[pd.notnull(housing['medListPriceSqft'])]
housing['log_price']=np.log(housing['medListPriceSqft'])
sample=nr.choice(housing['log_price'],size=100)

#contour plot of log likelihood
def likelihood_contour_plot(sample,location,scale):
    L, S = np.meshgrid(location, scale)
    Z=np.empty(shape=(len(scale),len(location)))
    for l in range(len(location)):
        for s in range(len(scale)):
            Z[s,l]=stats.norm.logpdf(sample, loc=location[l], scale=scale[s]).sum()

    fig, ax = plt.subplots(figsize=(7,5))
    ax.contour(L,S,Z,levels=500, cmap='RdGy')
    ax.set_xlabel('location')
    ax.set_ylabel('scale')
    return ax

#compute gradient
def gradient(mu,s2,x):
    n=x.shape[0]
    grad_mu = 1/s2*(x-mu).sum()
    grad_s2 = -n/(2*s2) + ((x-mu)**2).sum()/(2*s2**2)
    return [grad_mu, grad_s2]

#gradient descent
def gradient_descent_update(x, mu_start=2.0, s2_start=1.0, learning_rate=0.00001) :
    mu = [mu_start]
    s2 = [s2_start]
    grad = gradient(mu_start,s2_start,x)
    norm = [np.linalg.norm(grad)]

    for i in range(1000) :
        while norm[-1] >= 0.1 :
            grad = gradient(mu[-1],s2[-1],x)
            mu_i, s2_i = mu[-1] + learning_rate*grad[0], s2[-1] + learning_rate*grad[1]

            mu.append(mu_i)
            s2.append(s2_i)
            norm.append(np.linalg.norm(grad))

    return mu, s2, norm

#stochastic gradient descent
def sgd_update(x, batchsize=16, mu_start=2, s2_start=1, learning_rate=0.0001):
    mu = [mu_start]
    s2 = [s2_start]

    for i in range(1000):
        nr.shuffle(x)
        for start in range(0,len(x),batchsize):
            stop = start + batchsize
            mini_batch = x[start:stop]
            if (i==0 & start==0):
                grad = gradient(mu_start,s2_start,mini_batch)
                norm = [np.linalg.norm(grad)]
            else:
                while norm[-1] >= 0.01:
                    grad = gradient(mu[-1], s2[-1], mini_batch)
                    mu_i, s2_i = mu[-1] + learning_rate*grad[0], s2[-1] + learning_rate*grad[1]
                    mu.append(mu_i)
                    s2.append(s2_i)
                    norm.append(np.linalg.norm(grad))
    return mu, s2, norm

#perform gradient descent update
[mu, s2, norm] = gradient_descent_update(sample, mu_start=2.0, s2_start=1.0, learning_rate=0.00001)
df = pd.DataFrame({'iteration' : range(len(mu)), 'location' : mu, 's2' : s2, 'norm of gradient' : norm})
df['scale'] = np.sqrt(df['s2'])

#plot parameter values and norm of gradient against iteration
fig, [ax1,ax2,ax3] = plt.subplots(1,3,figsize=(15,5))
sns.scatterplot(x='iteration',y='location', data=df, ax=ax1)
sns.scatterplot(x='iteration',y='scale', data=df, ax=ax2)
sns.scatterplot(x='iteration',y='norm of gradient', data=df, ax=ax3)
ax2.set_ylim(0.3,)

#plot parameter values on contour plot of the likelihood function
location = np.arange(0.5,6,0.05)
scale = np.arange(0.15,2,0.05)
ax=likelihood_contour_plot(sample,location,scale)
ax.scatter(df['location'], df['scale'])

#perform stochastic gradient descent update
[mu, s2, norm] = sgd_update(sample, batchsize=16, mu_start=2.0, s2_start=1.0, learning_rate=0.0001)
df2 = pd.DataFrame({'iteration' : range(len(mu)), 'location' : mu, 's2' : s2, 'norm of gradient' : norm})
df2['scale'] = np.sqrt(df2['s2'])

#plot parameter values and norm of gradient against iteration
fig, [ax1,ax2,ax3] = plt.subplots(1,3,figsize=(15,5))
sns.scatterplot(x='iteration',y='location', data=df2, ax=ax1)
sns.scatterplot(x='iteration',y='scale', data=df2, ax=ax2)
sns.scatterplot(x='iteration',y='norm of gradient', data=df2, ax=ax3)
ax2.set_ylim(0.3,)
