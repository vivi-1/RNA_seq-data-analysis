#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 18:45:18 2020

@author: vivi
"""


# Group7= [UC1-UA1-UB1]   
# Group8= [DC1-DA1-DB1] 
# Group9= [Group3-DA1-DB1]=[UB1,UC1]-[DA1]-[DB1]         (Class two)
# Group10= [UB1,DC1]-[UA1]-[DA1] (Class three)
# Group11=[DB1,UC1]-[UA1]-[DA1]  (Class four)

# sheetname7: C: 6,10,14,18,22,26,30,,,,6+4*(k-1)
# sheetname9: C: 6,10,14,18,22,26,30,,,,6+4*(k-1)
# sheetname9: B: 5,9,13,17,21,25,29,,,,,5+4*(k-1)


import os
from os import listdir
currentDirectory_up="/Volumes/WD1/DEGlist/up_regulate"
currentDirectory_comparison="/Volumes/WD1/DEGlist/Comparison files/Gene_ID for each group"
os.chdir(currentDirectory_comparison)
os.getcwd()
os.listdir()

import re
import pandas
from pandas import ExcelFile
from pandas import ExcelWriter
import xlwt

Name_comparefile=""
Name_comparefile="Gene_ID_Compare_Group9.xlsx" ## Change dynamicly
comparefile=pandas.ExcelFile(Name_comparefile)
total_sheet_number=len(comparefile.sheet_names) #7
regex=re.compile(r'\d+')
k=int(regex.search(Name_comparefile).group(0))
k_str=str(k)
Search_file=[] ##Orginal file that we need to read and grep the values
wb=xlwt.Workbook(currentDirectory_comparison+'/'+'Gene sorted by log2FoldChange in C_group'+k_str+'.xls','wrb') ## Change dynamicly

annotation_path="/Volumes/WD1/DEGlist/Niben101_annotator.xlsx"
annotation_file=pandas.read_excel(annotation_path, sheet_name=0)

for i in range (0,total_sheet_number):
    sheet_names=str(i)
    ws=wb.add_sheet(sheet_names)
    comparefile_sheet=pandas.read_excel(Name_comparefile, i)   
    del comparefile_sheet['Unnamed: 0']
    row_number=len(comparefile_sheet.index)
    print("timepoint",i,"row number is "+str(row_number))
    searchfile_value=6+4*(i)
    searchfile_value_str=str(searchfile_value)
    Searchfile_name="G1vsG" + searchfile_value_str +"_DEG_up.xlsx"
    print(i)
    Search_file=pandas.read_excel(currentDirectory_up+"/"+Searchfile_name)
    Search_file.dropna(inplace = True)  
    print(Searchfile_name+" is read")
    for j in range(0,row_number):
        ID_value=comparefile_sheet.at[j, 0]
        value_to_write_row_number=Search_file.Gene_ID[Search_file.Gene_ID == ID_value].index.tolist()
        value_to_write_row=int(str(value_to_write_row_number).strip('[]'))
        value_to_write=Search_file.loc[value_to_write_row].at['log2FoldChange']
        annotation_to_write_row_number=annotation_file.Gene_ID[annotation_file.Gene_ID == ID_value].index.tolist()
        if annotation_to_write_row_number==[]:
            ws.write(j+1,2,"annotation doesn't exist")
            print(i,j,ID_value,"annotation doesn't exist")
        else:
            annotation_to_write_row=int(str(annotation_to_write_row_number).strip('[]'))
            annotation_value=annotation_file.loc[annotation_to_write_row].at['Unnamed: 4']
            ws.write(j+1,1,value_to_write)
            ws.write(j+1,0,ID_value)
            ws.write(j+1,2,annotation_value)
            print(i,j,"file is written")
    ws.write(0,0,"Gene_ID")
    ws.write(0,1,"log2FoldChange")
    ws.write(0,2,"annotation")
    print("sheet",i,"header is written")
wb.save('Gene sorted by log2FoldChange in C_group'+k_str+'.xls')
print("data saved")

## Rank each sheet based on log2FoldChange value
File_to_sort= pandas.ExcelFile('Gene sorted by log2FoldChange in C_group'+k_str+'.xls')
file_to_save=pandas.ExcelWriter('top200 by log2FoldChange in C_group'+k_str+'.xlsx', engine='xlsxwriter')
for i in range (0,len(File_to_sort.sheet_names)):
    sheet_to_read=pandas.read_excel(File_to_sort, i)
    sorted_by_log2 = sheet_to_read.sort_values(['log2FoldChange'], ascending=True)
    if len(sorted_by_log2.index)>=200:
        top200=sorted_by_log2.head(200)
        top200.to_excel(file_to_save, index=False, sheet_name=str(i))
    else:
        sorted_by_log2.to_excel(file_to_save, index=False, sheet_name=str(i))
    file_to_save.bookworksheet = file_to_save.sheets[str(i)]
file_to_save.save()