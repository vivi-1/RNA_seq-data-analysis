# -*- coding: utf-8 -*-
"""
Spyder Editor

author@Wei Wang
"""
import pandas as pd
import os

os.chdir('/Users/weiwang/Downloads/Enrichment')

def GOdic():
    GoDictionary={}
    for i in range(2,32):
        print(i)
        fileName="/Users/weiwang/Downloads/Enrichment/GO/G1vsG"+str(i)+"_all_GOenrich.xlsx"
        sheetName="G1vsG"+str(i)+"_all_GOenrich"
        GOList = pd.read_excel(fileName, sheet_name = sheetName, usecols=['GOID'])
        geneIDList=pd.read_excel(fileName, sheet_name = sheetName,usecols=['geneID'])
        CategoryList=pd.read_excel(fileName, sheet_name = sheetName, usecols=['Category'])

        GOList=GOList.stack().tolist()
        geneIDList=geneIDList.stack().tolist()
        CategoryList=CategoryList.stack().tolist()
        geneID=[]
        for element in geneIDList:
            temp = element.split("/")
            geneID.append(temp)
        print(geneID)
        for i in range(0, len(geneID)):
            print(geneID[i])
            if GOList[i] in GoDictionary:
                for j in geneID[i]:
                    print(j, GoDictionary[GOList[i]])
                    GoDictionary[GOList[i]].append(j)
            if GOList[i] not in GoDictionary:
                for j in geneID[i]:
                    GoDictionary[GOList[i]]=[j]
    return GoDictionary

def createGeneGO(GoDictionary):
    geneGOList=[]
    for key, value in GoDictionary.items():
        for i in value:
            geneGOList.append([i, key])
    return geneGOList

GoDictionary=GOdic()
geneGoDic = createGeneGO(GoDictionary)

import xlsxwriter
workbook = xlsxwriter.Workbook('/Users/weiwang/Desktop/GO annotations.xlsx')
worksheet = workbook.add_worksheet("GO terms")
row=0

for i in geneGoDic:
    if row==0:
        worksheet.write(row, 0, "geneID")
        worksheet.write(row, 1, "GOID")
    if row>0:
        worksheet.write(row, 0, i[0])
        worksheet.write(row, 1, i[1])
    row+=1
workbook.close()
