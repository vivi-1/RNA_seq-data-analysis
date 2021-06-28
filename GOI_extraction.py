# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import pandas as pd
from pandas import *
import os
os.chdir('/Volumes/WD1/Desktop/laboratory files/Results/Altria project/DEGlist/Wei_reDEGlist/Output')
filePath="/Volumes/WD1/Desktop/laboratory files/Results/Altria project/DEGlist/Highlighted genes from P.palmivora.xlsx"
geneList = pd.read_excel(filePath, sheet_name="sequence and expression info", usecols=['Gene_ID'])
geneList=geneList.stack().tolist()
FClist=[]
for i in geneList:
    for j in range(1,8):
        fileName = 'T'+str(j) +' _all_Pnic_no rep2.csv'
        file=pd.read_csv(fileName)
        geneName=file["Row.names"].tolist()
        fcList=file["log2FoldChange"].tolist()
        plist = file["padj"].tolist()
        if i in geneName:
            index = geneName.index(i)
            result=fcList[index]
            pv = plist[index]
        if i not in geneName:
            result="NA"
            pv="NA"
        FClist.append([i,j,result,pv])

print(FClist)

import xlsxwriter
workbook = xlsxwriter.Workbook('/Users/weiwang/Desktop/GOI RNAseq.xlsx')
worksheet = workbook.add_worksheet("My sheet")
row=1
for i in FClist:
    worksheet.write(row, 0, i[0])
    worksheet.write(row, 1, i[1])
    worksheet.write(row, 2, i[2])
    worksheet.write(row, 3, i[3])
    row+=1
workbook.close()
