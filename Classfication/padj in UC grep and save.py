#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 14:38:38 2020

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
import numpy
currentDirectory_up="/Volumes/WD1/DEGlist/up_regulate"
currentDirectory_comparison="/Volumes/WD1/DEGlist/Comparison files/comparison to other studies"
os.chdir(currentDirectory_comparison)
os.getcwd()
os.listdir()

import re
import pandas
import xlwt
Name_comparefile=""
Name_comparefile="comparison data_Group7.xls" ## Change dynamicly
read_comparefile=pandas.read_excel(Name_comparefile)
read_comparefile=read_comparefile.drop('Unnamed: 0', axis=1)
read_comparefile_col=len(read_comparefile.columns) #7 clumns
regex=re.compile(r'\d+')
k=int(regex.search(Name_comparefile).group(0))
sheetname1=str(k)
Search_file=[] ##Orginal file that we need to read and grep the values
wb=xlwt.Workbook(currentDirectory_comparison+'/'+'padj in C_ups_group7.xls','wrb') ## Change dynamicly
ws=wb.add_sheet(sheetname1)


## Gene_ID in other studies
import xlrd
Name_newfile=currentDirectory_comparison+'/'+'Gene_ID in other study.xlsx'
highlight_workbook=xlrd.open_workbook(Name_newfile,'rb')
highlight_sheet=highlight_workbook.sheet_by_index(0)
highlight_ID=[]
for i in range(0,highlight_sheet.nrows):
    for j in range(0,highlight_sheet.ncols): 
        j=0 ## Read the first column
    highlight_ID.append(highlight_sheet.cell(i,j).value)
highlight_ID.remove('Gene_ID')
## Create data file

ID_value=[]
searchfile_value_str=""
## Read the cell_value in compare file
for i in range(0,len(highlight_ID)):
    ID_value=highlight_ID[i]
    for j in range(0, read_comparefile_col):
        if read_comparefile.iat[i,j].isnumeric()==True:
            searchfile_value=6+4*(j)
            searchfile_value_str=str(searchfile_value)
            Searchfile_name="G1vsG" + searchfile_value_str +"_DEG_up.xlsx"
            Search_file=pandas.read_excel(currentDirectory_up+"/"+Searchfile_name)
            Search_file.dropna(inplace = True) 
            print(i,j)
            value_to_write_row_number=Search_file.Gene_ID[Search_file.Gene_ID == ID_value].index.tolist()
            value_to_write_row=int(str(value_to_write_row_number).strip('[]'))
            print(value_to_write_row)
            value_to_write=Search_file.loc[value_to_write_row].at['padj']
            print(value_to_write)
            ws.write(i+1,j+1,value_to_write)
        else:
            value_to_write='N'
            print(i,j)
            print(value_to_write)
            ws.write(i+1,j+1,value_to_write)
        print("saved data")              
wb.save('padj in C_ups_group7.xls') ## Change dynamicly
print("data saved")

##Check how long does it take to run the code
import time
start_time = time.time()
print("--- %s seconds ---" % (time.time() - start_time))      