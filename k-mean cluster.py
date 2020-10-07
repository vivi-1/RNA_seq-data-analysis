#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 22 19:25:10 2020

@author: vivi
"""
import pandas
import os

current_dir='/Users/vivi/Desktop'
os.chdir(current_dir)
df1 = pandas.read_excel("FPKM values in C files_combined.xlsx")

df1.shape


## Niben101Scf08195g03020 or Niben101Ctg07414g00002 or Niben101Ctg1209t00020
#df1=df1.replace(regex=['Ctg'], value='999')
#df1=df1.replace(regex=['Scf'], value='888')
#df1=df1.replace(regex=['Niben101'], value='')
#df1=df1.replace(regex=['g'], value='1') ## change g to 1
#df1=df1.replace(regex=['t'], value='2') ## change g to 2

list1=df1['Gene_ID']
del df1['Gene_ID']

data=df1.to_numpy()
data=data.transpose()

import numpy as np

def kmeans_xufive(data, k):
    """k-means聚类算法 
    k       - 指定分簇数量
    data      - ndarray(m, n)，m个样本的数据集，每个样本n个属性值
    """
    
    m, n = data.shape # m：样本数量，n：每个样本的属性值个数
    result = np.empty(m, dtype=np.int) # m个样本的聚类结果
    cores = data[np.random.choice(np.arange(m), k, replace=False)] # 从m个数据样本中不重复地随机选择k个样本作为质心
    
    while True: # 迭代计算
        d = np.square(np.repeat(data, k, axis=0).reshape(m, k, n) - cores)
        distance = np.sqrt(np.sum(d, axis=2)) # ndarray(m, k)，每个样本距离k个质心的距离，共有m行
        index_min = np.argmin(distance, axis=1) # 每个样本距离最近的质心索引序号
        
        if (index_min == result).all(): # 如果样本聚类没有改变
            return result, cores # 则返回聚类结果和质心数据
        
        result[:] = index_min # 重新分类
        for i in range(k): # 遍历质心集
            items = data[result==i] # 找出对应当前质心的子样本集
            cores[i] = np.mean(items, axis=0) # 以子样本集的均值作为当前质心的位置


import matplotlib.pyplot as plt
k = 3
result, cores = kmeans_xufive(data, k)
plt.scatter(data[:,0], data[:,1], s=1, c=result.astype(np.int)+1)
#plt.scatter(data[:,2], data[:,3], s=1, c='g', label='rep3 and rep4')
plt.scatter(cores[:,0], cores[:,1], marker='x', c=np.arange(k))
plt.show()

k = 2
result, cores = kmeans_xufive(data, k)
plt.scatter(data[:,0], data[:,1], s=1, c=result.astype(np.int)+1)
#plt.scatter(data[:,2], data[:,3], s=1, c='g', label='rep3 and rep4')
plt.scatter(cores[:,0], cores[:,1], marker='x', c=np.arange(k))
plt.show()

k = 4
result, cores = kmeans_xufive(data, k)
plt.scatter(data[:,0], data[:,1], s=1, c=result.astype(np.int)+1)
#plt.scatter(data[:,2], data[:,3], s=1, c='g', label='rep3 and rep4')
plt.scatter(cores[:,0], cores[:,1], marker='x', c=np.arange(k))
plt.show()

k=1
result, cores = kmeans_xufive(data, k)
plt.scatter(data[:,0], data[:,1], s=1, c=result.astype(np.int)+1)
#plt.scatter(data[:,2], data[:,3], s=1, c='g', label='rep3 and rep4')
plt.scatter(cores[:,0], cores[:,1], marker='x', c=np.arange(k))
plt.show()



### sklearn cluster

import numpy as np
from sklearn import cluster
from sklearn.cluster import KMeans
estimator=KMeans(n_clusters=2)
estimator.fit(data)
label_pred=estimator.labels_
centroids=estimator.cluster_centers_
inertia=estimator.inertia_

import numpy as np
from sklearn import cluster
from sklearn.cluster import KMeans
estimator=KMeans(n_clusters=3)
estimator.fit(data)
label_pred=estimator.labels_
centroids=estimator.cluster_centers_
inertia=estimator.inertia_

import numpy as np
from sklearn import cluster
from sklearn.cluster import KMeans
estimator=KMeans(n_clusters=4)
estimator.fit(data)
label_pred=estimator.labels_
centroids=estimator.cluster_centers_
inertia=estimator.inertia_

## 
num_clusters=1
km_cluster=KMeans(n_clusters=num_clusters, max_iter=3000, init='k-means++', n_jobs=-1)
result=km_cluster.fit_predict(data)
km_cluster.labels_
km_cluster.predict(data)
print (result)
            
            
