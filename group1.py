#!/user/bin/python
#####Developed by Kevin Yu and Wei
################################## File manilupation ##########################
## Go to the directory containing all DEG files
import os
import glob
# import natsort
import pandas
import xlrd
import csv
import io

## FileExtensions should match up with folderNames and have the same size
currentDirectory = "/Volumes/WD1/DEGlist"
fileExtensions = ["*up.xls","*down.xls","*all.xls"]
folderNames = ['up_regulate','down_regulate','all']
nameOfFiles = []

os.getcwd()
os.chdir(currentDirectory)

## Finding up, down, & all gene files and put them in a new directory
## Obtaining the names of all current files (from fileExtension)
for i in range(0,len(fileExtensions)):
     tempFileNames = []
     for filename in glob.glob(fileExtensions[i]):
          tempFileNames.append(filename)
     nameOfFiles.append(tempFileNames)

## creating the directories (from folderNames)
for i in range(0,len(folderNames)):
     if not os.path.exists(folderNames[i]):
          os.makedirs(folderNames[i])

## moving all current files in currentDirecotry to currentDirectory+folderNames
for i in range(0,len(folderNames)):
     for j in range(0,len(nameOfFiles[i])):
          d1 = currentDirectory + "/" + nameOfFiles[i][j]
          d2 = currentDirectory + "/" + folderNames[i] + "/" + nameOfFiles[i][j]
          os.rename(d1,d2)

################################## Data manilupation #######################
## list files in upregulate:  os.listdir('/Volumes/WD1/DEGlist/up_regulate')
## sort files by names in Spyder

up_dir = '/Volumes/WD1/DEGlist/up_regulate'
down_dir = '/Volumes/WD1/DEGlist/down_regulate'
up_files = os.listdir(up_dir)
# print(up_files)
down_files = os.listdir(down_dir)
# natsort(up_files)
# os.chdir(currentDirectory + "/" + folderNames[i])

# data1 = []
# for i in range(0,len(folderNames)):
#     os.chdir(currentDirectory + "/" + folderNames[i])
#     tempFileNames = []
#     for filename in glob.glob(fileExtensions[i] + "x"):
#             tempFileNames.append(filename)
#     for j in range(0,len(tempFileNames)):
#         # d1 = currentDirectory + "/" + nameOfFiles[i][j]
#         d2 = currentDirectory + "/" + folderNames[i] + "/" + tempFileNames[j]
#         currentWorkbook = xlrd.open_workbook(d2)
#         test1 = pandas.read_excel(d2)
#         print(d2)
#         print(test1.shape)


data1 = []
for i in range(0,len(folderNames)):
    os.chdir(currentDirectory + "/" + folderNames[i])
    tempFileNames = []
    for filename in glob.glob(fileExtensions[i] + "x"):
            tempFileNames.append(filename)
    for j in range(0,len(tempFileNames)):
        # d1 = currentDirectory + "/" + nameOfFiles[i][j]
        d2 = currentDirectory + "/" + folderNames[i] + "/" + tempFileNames[j]
        currentWorkbook = xlrd.open_workbook(d2)
        test1 = pandas.read_excel(d2)
        print(d2)
        print(test1.shape)


## Organize the file matrix by the order of name
# import re
# numbers = re.compile(r'(\d+)')
# def numericalSort(FOI):`
    # parts = numbers.split(FOI)
    # parts[1::2] = map(int, parts[1::2])
    # return parts
# upfile_sorted=sorted(up_files, key = numericalSort)
# upfile_splice=upfile_sorted[2:]
# downfile_sorted = sorted(down_files, key = numericalSort)
# downfile_splice=downfile_sorted[2:]
## Group the files from g4-g7, G8-G11, G12-G15, G16-G19, G20-G23, G24-G27, G28-G31

    