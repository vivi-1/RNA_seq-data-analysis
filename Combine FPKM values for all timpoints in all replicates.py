#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 17:21:00 2020

@author: vivi
"""

import pandas
import os

current_dir='/Users/vivi/Desktop'
os.chdir(current_dir)
read_file= pandas.ExcelFile('C1-C7-four replicates.xlsx')

for i in (0,len(sheet_names))
df1 = pandas.read_excel("C1-C7-four replicates.xlsx")

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
