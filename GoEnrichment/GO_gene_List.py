#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 12:50:22 2021

@author: weiwang
"""
import pandas as pd
import xlsxwriter
#create annotation of GO
def readInGOlist():
    fileName="/Users/weiwang/Desktop/GOterms annotations.xlsx"
    GOList = pd.read_excel(fileName, sheet_name = "GO terms", usecols=['GoID'])
    GOList=GOList.stack().tolist()
    return GOList

def getGenelistWithGo():
    GOs=readInGOlist()
    fileName="/Volumes/WD1/Desktop/laboratory files/Results/Altria project/DEGlist/Wei_reDEGlist/Input/Niben101_annotator.xlsx"
    GOList = pd.read_excel(fileName, sheet_name = "Gene_GOCluster", usecols=['GOID'])
    GOList=GOList.stack().tolist()
    GeneList = pd.read_excel(fileName, sheet_name = "Gene_GOCluster", usecols=['GeneID'])
    GeneList=GeneList.stack().tolist()
    if len(GOList)!= len(GeneList):
        print(len(GOList), len(GeneList))

    GoIDListWithGeneList=[]
    cnt=0
    lenGO=len(GOs)
    for GO in GOs:
        cnt+=1
        print(cnt,lenGO)
        for i in range(0, len(GOList)):
            if GO in GOList[i]:
                GoIDListWithGeneList.append([GO, GeneList[i]])
    return GoIDListWithGeneList

GoIDListWithGeneList=getGenelistWithGo()

def getGenelistWithGo():
    workbook = xlsxwriter.Workbook('/Users/weiwang/Desktop/GOterms_GeneList.xlsx')
    worksheet = workbook.add_worksheet("GO terms")
    row=0
    print("saving...")
    for i in GoIDListWithGeneList:
        print(row)
        if row==0:
            worksheet.write(row, 0, "GoID")
            worksheet.write(row, 1, "GeneID")
        if row>0:
            worksheet.write(row, 0, i[0])
            worksheet.write(row, 1, i[1])
        row+=1
    workbook.close()

getGenelistWithGo()
