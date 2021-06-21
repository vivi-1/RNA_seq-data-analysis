#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 20 23:38:50 2021

@author: weiwang
"""

import xlrd
import pandas as pd
from pandas import DataFrame

def containKeyword(original, keyword):

    if original.find(keyword) == -1:
        print("not contain the keyword")
        return False
    else:
        print("Found the keyword in the string.")
        return True

sheet = []
results = []
annotations = []
wb = xlrd.open_workbook("/Volumes/WD1/Desktop/laboratory files/Results/Altria project/DEGlist/Wei_reDEGlist/Output/Eliminate outlier rep2/Class1_Gene_ID_Annotation.xlsx")

for i in range (0,7):
    sheet.append(wb.sheet_by_index(i))

for i in range(0,7):
    list = []
    #print(i, sheet[i].nrows)
    for j in range(1, sheet[i].nrows):
        list.append(sheet[i].cell_value(j,2))
    annotations.append(list)

keywords = ["sulf", "hypersensitive", "death", "cys", "glutathione", "GSL", "camalexin", "redox", "salicylic"]

for keyword in keywords:
    result=[]
    temp=[]
    for timepoint in range(0,len(annotations)):
        fileName = "/Volumes/WD1/Desktop/laboratory files/Results/Altria project/DEGlist/Wei_reDEGlist/Output/Eliminate outlier rep2/" + "T"+ str(timepoint+1) + "_up_Pnic_no rep2.csv"
        file = pd.read_csv(fileName, usecols= ['Gene_ID', 'log2FoldChange'])
        geneIDs = file['Gene_ID'].tolist()
        print(geneIDs)
        foldChange = file['log2FoldChange'].tolist()
        
        for j in range(0,len(annotations[timepoint])):            
            if containKeyword(annotations[timepoint][j], keyword):
                index = geneIDs.index((sheet[timepoint].cell_value(j,1)))                
                temp = [timepoint+1, sheet[timepoint].cell_value(j,1), foldChange[index], annotations[timepoint][j]]
                result.append(temp)
                print("hi")
    results.append(result)    


for i in range(0,len(results)):
    outputName="/Volumes/WD1/Desktop/laboratory files/Results/Altria project/DEGlist/Wei_reDEGlist/Output"+keywords[i]+'.xlsx'
    df = DataFrame(results[i])
    df = df.transpose()
    for j in range(0,len(results[i])):        
        df.to_excel(outputName)  
    
    
    
    
    
