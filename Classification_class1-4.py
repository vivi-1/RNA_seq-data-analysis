#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 16 19:16:31 2021

@author: weiwang
"""
import numpy
import pandas
from pandas import DataFrame
import os

def intersectGenes(A,B): # AnB
    return numpy.intersect1d(A,B)

def createExcelFile(currentDirectory,fileName): # creates excel file to manipulate
    newFileName = currentDirectory + "/" + fileName
    outputFile = pandas.ExcelWriter(newFileName, engine='xlsxwriter')
    return outputFile
#Class one: [UPnic-Uflag22]   
#Class two: [Uflag22,UPnic]-[Dflag22]
#Class three: [Uflag22,DPnic]
#Class four: [Dflag22,UPnic]
print(os.getcwd())
os.chdir('/Users/weiwang/Desktop/DEGlist/Output/Eliminate outlier rep2')
print(os.getcwd())

def class1(currentDirectory):
    excelFile = createExcelFile(currentDirectory,"Class1_Gene_ID.xlsx")
    for i in range(1, 8):
        print(i)
        UPnicName = "T" + str(i) + "_up_Pnic_no rep2.csv"
        Uflag22Name = "T" + str(i) + "_up_Flag22_no rep2.csv"
        print(UPnicName)
        print(Uflag22Name)
        UPnic = pandas.read_csv(UPnicName)
        UPnic = UPnic['Gene_ID']
        len1 = len(UPnic)
        print(UPnic)
        print(len1)
        Uflag22 = pandas.read_csv(Uflag22Name)
        Uflag22 = Uflag22['Gene_ID']
        print(Uflag22)
        result = []
        for j in range(0, len1):
            if UPnic[j] not in Uflag22:
                result.append(UPnic[j])        
        print(result)
        DataFrame({"Gene_ID":result}).to_excel(excelFile,sheet_name=str(i))
    excelFile.save()

def class2(currentDirectory):#Class two:[Uflag22,UPnic]-[Dflag22]
    excelFile = createExcelFile(currentDirectory,"Class2_Gene_ID.xlsx")
    for i in range(1, 8):
        print(i)
        UPnicName = "T" + str(i) + "_up_Pnic_no rep2.csv"
        Uflag22Name = "T" + str(i) + "_up_Flag22_no rep2.csv"
        Dflag22Name = "T" + str(i) + "_down_Flag22_no rep2.csv"
        print(UPnicName)
        print(Uflag22Name)
        UPnic = pandas.read_csv(UPnicName)
        UPnic = UPnic['Gene_ID']
        len1 = len(UPnic)
        print(UPnic)
        print(len1)
        Uflag22 = pandas.read_csv(Uflag22Name)
        Uflag22 = Uflag22['Gene_ID']
        print(Uflag22)
        Dflag22 = pandas.read_csv(Dflag22Name)
        Dflag22 = Dflag22['Gene_ID']
        print(Dflag22)
        result = []
        for j in range(0, len1):
            if (UPnic[j] in Uflag22) and (UPnic[j] not in Dflag22):
                result.append(UPnic[j])        
        print(result)
        
        DataFrame({"Gene_ID":result}).to_excel(excelFile,sheet_name=str(i))
    excelFile.save()

def class3(currentDirectory):#[Uflag22,DPnic]
    excelFile = createExcelFile(currentDirectory,"Class3_Gene_ID.xlsx")
    for i in range(1, 8):
        print(i)
        DPnicName = "T" + str(i) + "_down_Pnic_no rep2.csv"
        Uflag22Name = "T" + str(i) + "_up_Flag22_no rep2.csv"

        print(DPnicName)
        print(Uflag22Name)
        DPnic = pandas.read_csv(DPnicName)
        DPnic = DPnic['Gene_ID']
        len1 = len(DPnic)
        print(DPnic)
        print(len1)
        Uflag22 = pandas.read_csv(Uflag22Name)
        Uflag22 = Uflag22['Gene_ID']
        print(Uflag22)

        result = []
        for j in range(0, len1):
            print(DPnic[j])
            if DPnic[j] in Uflag22:
                print(DPnic[j])
                result.append(DPnic[j])        
        print(result)
        
        DataFrame({"Gene_ID":result}).to_excel(excelFile,sheet_name=str(i))
    excelFile.save()

def class4(currentDirectory):#[Dflag22,UPnic]
    excelFile = createExcelFile(currentDirectory,"Class4_Gene_ID.xlsx")
    for i in range(1, 8):
        print(i)
        UPnicName = "T" + str(i) + "_up_Pnic_no rep2.csv"
        Dflag22Name = "T" + str(i) + "_down_Flag22_no rep2.csv"
        print(UPnicName)
        print(Dflag22Name)
        
        UPnic = pandas.read_csv(UPnicName)
        UPnic = UPnic['Gene_ID']
        len1 = len(UPnic)
        print(UPnic)
        print(len1)
        Dflag22 = pandas.read_csv(Dflag22Name)
        Dflag22 = Dflag22['Gene_ID']
        print(Dflag22)

        result = []
        for j in range(0, len1):
            print(UPnic[j])
            if UPnic[j] in Dflag22:
                print(UPnic[j])
                result.append(UPnic[j])   
            else:
                result.append('N')
        print(result)
        
        DataFrame({"Gene_ID":result}).to_excel(excelFile,sheet_name=str(i))
    excelFile.save()

class1('/Users/weiwang/Desktop/DEGlist/Output/Eliminate outlier rep2')
class2('/Users/weiwang/Desktop/DEGlist/Output/Eliminate outlier rep2')
class3('/Users/weiwang/Desktop/DEGlist/Output/Eliminate outlier rep2')
class4('/Users/weiwang/Desktop/DEGlist/Output/Eliminate outlier rep2')


