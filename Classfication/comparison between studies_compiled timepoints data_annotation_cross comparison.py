#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 21:03:13 2020

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

#E compare highlight files with group7 file 
group_name_variable='Group7'
Name_groupfile='Gene_ID_Compare_'+group_name_variable+'.xlsx'
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
    read_sheet=pandas.read_excel(read_groupfile, index_col=None, usecols="B", sheet_name=i)
    read_sheet_list=read_sheet.values.tolist()
    for j in range(0,len(read_sheet_list)):
        read_value=str(read_sheet_list[j]).strip('[]')
        read_value=str(read_value).strip("''")
        if read_value not in list_temp:
            list_temp.append(read_value)
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

## compile all genes in group7, saved in list_temp list
wb=[]
ws=[]
wb=xlwt.Workbook('all genes_timepoints combined group'+sheetname1+'with annotation.xls','wrb')
ws=wb.add_sheet(sheetname1)
print("sheet"+ sheetname1 + "is created")
read_groupfile = pandas.ExcelFile(Name_groupfile)
for j in range(0,len(read_groupfile.sheet_names)):
    read_sheet=pandas.read_excel(read_groupfile, index_col=None, usecols="B", sheet_name=j)
    read_sheet_list=read_sheet.values.tolist()
    value_to_write=str(j+1)
    for i in range(0,len(list_temp)):
        ID_value=[list_temp[i]]
        if ID_value in read_sheet_list:           
            ws.write(i+1,j+1,value_to_write)
            print(i,j,'is written')
        else:
            ws.write(i+1,j+1,'N')
            print(i,j,'N is written')
ws.write(0,0,'Gene_ID')
ws.write(0,8,'annotation')

## Add annotation to the compiled group7 genes
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
wb.save('all genes_timepoints combined group'+sheetname1+' with annotation.xls')
print("data is written in "+"sheet"+sheetname1)

## Read P.parasitica file, data downloaded from: https://link.springer.com/article/10.1007/s10725-016-0163-1#Sec17
searchfile_name='Nbenthi under P.parasitica_NBC to NB101.xls'
searchfile_dir="/Volumes/WD1/DEGlist/Comparison files/comparison to other studies/"
Name_searchfile=searchfile_dir+searchfile_name
Name_searchfile_workbook=xlrd.open_workbook(Name_searchfile,'rb')
search_sheet=Name_searchfile_workbook.sheet_by_index(0)
search_ID=[]
for i in range(0,search_sheet.nrows):
    for j in range(0,search_sheet.ncols): 
        j=1 ## Read the second column
    search_ID.append(search_sheet.cell(i,j).value)
search_ID=search_ID[1:]

## compare to p.parasitica file
wb=[]
ws1=[]
wb=xlwt.Workbook('all genes_comparison with P.parasitica_Group'+sheetname1+'.xls','wrb')
ws1=wb.add_sheet(sheetname1)
print("sheet"+ sheetname1 + "is created")
for i in range(0,len(list_temp)):
    ws1.write(i+1,0,list_temp[i])
    if list_temp[i] in search_ID:
        ws1.write(i+1,9,'Y')
    else:
        ws1.write(i+1,9,'N')
ws1.write(0,9,'If the gene showed up in Parasitica file')
wb.save('all genes_comparison with P.parasitica_Group'+sheetname1+'.xls')
    




