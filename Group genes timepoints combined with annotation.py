#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 17:07:47 2020

@author: vivi
"""

import os
import xlrd
import re
import xlwt
import pandas
# from xlrd import load_workbook
currentDirectory = "/Volumes/WD1/DEGlist/Comparison files/Gene_ID for each group"
os.getcwd()
os.chdir(currentDirectory)
os.listdir()

#Example compare highlight files with group7 file 
group_name_variable='group9'
Name_groupfile='top200 by log2FoldChange in C_'+group_name_variable+'.xlsx'
read_groupfile=xlrd.open_workbook(Name_groupfile,'rb')
sheetnumber_groupfile=len(read_groupfile.sheet_names())
regex=re.compile(r'\d+')
k=int(regex.search(group_name_variable).group(0))
sheetname1=str(k)
list_temp=[]
row_number=0
annotation_path="/Volumes/WD1/DEGlist/Niben101_annotator.xlsx"
annotation_file=pandas.read_excel(annotation_path, sheet_name=0)
for i in range(0,sheetnumber_groupfile):
    read_sheet=read_groupfile.sheet_by_index(i)
    row_number_start=read_groupfile.sheet_by_index(0).nrows-1
    for j in range(1,read_sheet.nrows):
        if read_sheet.cell(j,0).value not in list_temp:
            list_temp.append(read_sheet.cell(j,0).value)
            row_number=row_number+1
            print(i,j)
            print(row_number)
        else:
            list_temp=list_temp
            row_number=row_number
            print(i,j)
            print(row_number)
if row_number==len(list_temp):
    print("total row number is correct, which is "+str(row_number))
    row_number=0
    print("row_number is cleared to 0")

## Find which timepoints does these genes showed up
wb=[]
ws=[]
wb=xlwt.Workbook('timepoints combined group'+sheetname1+'with annotation.xls','wrb')
ws=wb.add_sheet(sheetname1)
print("sheet"+ sheetname1 + "is created")
read_groupfile = pandas.ExcelFile(Name_groupfile)
for j in range(0,len(read_groupfile.sheet_names)):
    read_sheet=pandas.read_excel(read_groupfile, index_col=None, usecols="A", sheet_name=j)
    read_sheet_list=read_sheet.values.tolist()
    value_to_write=str(j+1)
    for i in range(0,len(list_temp)):
        ID_value=[list_temp[i]]
        annotation_to_write_row_number=annotation_file.Gene_ID[annotation_file.Gene_ID == list_temp[i]].index.tolist()
        if ID_value in read_sheet_list:           
            ws.write(i+1,j+1,value_to_write)
        else:
            ws.write(i+1,j+1,'N')
ws.write(0,0,'Gene_ID')
ws.write(0,8,'annotation')
for i in range(0,len(list_temp)):
    ws.write(i+1,0,list_temp[i])
    annotation_to_write_row_number=annotation_file.Gene_ID[annotation_file.Gene_ID == list_temp[i]].index.tolist()
    if annotation_to_write_row_number==[]:
        ws.write(i+1,8,"annotation doesn't exist")
        print(i,ID_value,"annotation doesn't exist")
    else:
        annotation_to_write_row=int(str(annotation_to_write_row_number).strip('[]'))
        annotation_value=annotation_file.loc[annotation_to_write_row].at['Unnamed: 4']
        ws.write(i+1,8,annotation_value)
        print(i,'annotation is',annotation_value)
for j in range(0,len(read_groupfile.sheet_names)):
    value="timepoint"+str(j+1)
    ws.write(0,j+1,value)

wb.save('timepoints combined group'+sheetname1+' with annotation.xls')
print("data is written in "+"sheet"+sheetname1)


