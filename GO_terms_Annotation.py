#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 12:50:22 2021

@author: weiwang
"""
import pandas as pd
import xlsxwriter
#create annotation of GO
def Go_annotation():
    fileName="/Volumes/WD1/Desktop/laboratory files/Results/Altria project/DEGlist/Wei_reDEGlist/Input/Niben101_annotator.xlsx"
    xl = pd.ExcelFile(fileName)
    print(xl.sheet_names)  # see all sheet names


    GOList = pd.read_excel(fileName, sheet_name = "GO_terms", usecols=['GoID'])
    GOList=GOList.stack().tolist()
    print(GOList)
    GoIDListWithAnnotation=[]
    for element in GOList:
        temp = element.split("), ")
        for i in temp:
            GoIDListWithAnnotation.append(i)
    GoIDs=[]

    for element in GoIDListWithAnnotation:
        temp = element.split(" (")
        GoIDs.append([temp[0].replace(' ',''), temp[1].replace(')','')])
    GoIDs.pop()
    print(GoIDs)
    return GoIDs

def writeGOAnnotation():
    GoIDs=Go_annotation()
    workbook = xlsxwriter.Workbook('/Users/weiwang/Desktop/GOterms annotations.xlsx')
    worksheet = workbook.add_worksheet("GO terms")
    row=0
    for i in GoIDs:
        if row==0:
            worksheet.write(row, 0, "GoID")
            worksheet.write(row, 1, "Annotation")
        if row>0:
            worksheet.write(row, 0, i[0])
            worksheet.write(row, 1, i[1])
        row+=1
    workbook.close()


writeGOAnnotation()
