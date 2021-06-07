#!/user/bin/python
################################## File manilupation ##########################
## Go to the directory containing all DEG files
import os
import glob
import pandas
import re
import numpy
import xlrd

def numericalSort(FOI): # obtains the VS number for each file
    numbers = re.findall(r"\d+", FOI)
    return int(numbers[1])
    
def getData(folder,letter,data): # folder = 'up_regulate' or 'down_regulate'; letter = 0:A, 1:B, 2:C, 3:D; data = dataOfDataStucturesNames
    output = []
    for i in range(0,len(data[folder])):
        tempData = data[folder][i][letter]
        output.append(tempData)
    return output

def subtractGenes(A,B): # A-B
    return set(A)-set(B)

def intersectGenes(A,B): # AnB
    return numpy.intersect1d(A,B)


def createExcelFile(currentDirectory,fileName): # creates excel file to manipulate
    newFileName = currentDirectory + "/" + fileName
    outputFile = pandas.ExcelWriter(newFileName, engine='xlsxwriter')
    return outputFile

################################## Global Variables for Data in File Explorer #######################
## FileExtensions should match up with folderNames and have the same size
currentDirectory = "/Volumes/WD1/DEGlist"
fileExtensions = ["*up.xls","*down.xls","*all.xls"]
folderNames = ['up_regulate','down_regulate','all']
nameOfFiles = []
os.getcwd()
os.chdir(currentDirectory)

################################## Creating Sub-folders and Moving Files ############################
## Finding up, down, & all gene files and put them in a new directory
## Obtaining the names of all current files (from fileExtension)
for i in range(0,len(fileExtensions)):
     tempFileNames = []
     for fileName in glob.glob(fileExtensions[i]):
          tempFileNames.append(fileName)
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

################################## Obtaining Files and Data #########################################
## list files in upregulate:  os.listdir('/Volumes/WD1/DEGlist/up_regulate')
## sort files by names
fullDataStructureIndices = []
fileNames = []
## loops through "folderNames". We do -1 because we don't care about the all folder
for i in range(0,len(folderNames)-1):
    os.chdir(currentDirectory + "/" + folderNames[i])
    tempFileNames = []
    for fileName in glob.glob(fileExtensions[i] + "x"): # gets all the file names
        tempFileNames.append(fileName)
    fileNames.append(tempFileNames)
    ## Group the files from G4-G7, G8-G11, G12-G15, G16-G19, G20-G23, G24-G27, G28-G31
    count = 4
    indexOfAllFileNames = []
    while count < 32:
        indexOfFileNames = []
        while len(indexOfFileNames) <4:
            for j in range(0,len(tempFileNames)):
                if numericalSort(tempFileNames[j]) == count:
                    indexOfFileNames.append(j)
                    count = count+1
                if len(indexOfFileNames) >3:
                    break
                if count > 31:
                    break
        indexOfAllFileNames.append(indexOfFileNames)
    fullDataStructureIndices.append(indexOfAllFileNames)

## gets the data from all of the files in a clean format
dataStructureNames = {}
dataOfDataStucturesNames = {}
for i in range(0,len(folderNames)-1):
    tempDictionaryNames = folderNames[i]
    tempData = {}
    tempData3 = {}
    for j in range(0,7):
        groupNumber = j
        tempData2 = []
        tempData4 = []
        for k in range(0,len(fullDataStructureIndices[i][j])):
            tempData2.append(currentDirectory + "/" + folderNames[i] + "/" + fileNames[i][fullDataStructureIndices[i][j][k]])
            fileData = pandas.read_excel(currentDirectory + "/" + folderNames[i] + "/" + fileNames[i][fullDataStructureIndices[i][j][k]]).fillna(0)
            geneNames=fileData.Gene_ID
            tempData4.append(geneNames)
        tempData[groupNumber] = tempData2
        tempData3[groupNumber] = tempData4
    dataStructureNames[tempDictionaryNames] = tempData
    dataOfDataStucturesNames[tempDictionaryNames] = tempData3


## Get file mitrixes as Up_regulate or down_regulate
#                      [A1,B1,C1,D1]   [G4,G5,G6,G7]
#                      [A2,B2,C2,D2]   [G8,G9,G10,G11]
#                      ...            ...
#                      [A7,B7,C7,D7]   [G28,G29,G30,G31]
# Group1=Group1_up=[UC1, UD1] = [U21,U31]
# Group2=Group1_down=[DC1,DD1]     
## Group3=Group2_up = [UB1,UC1]     
# Group4=Group2_down=[DB1,DC1]     
# Group5=Group3_original=[Group2_up, Group1_down]
# Group6=Group4_original=[Group1_up, Group2_down]
##Group7=[UC1-UA1-UB1]   
##Group8=[DC1-DA1-DB1] 
# Group9=[Group3-DA1-DB1]=[UB1,UC1]-[DA1]-[DB1] (Class two)
# Group10=[UB1,DC1]-[UA1]-[DA1] (Class three)
# Group11=[DB1,UC1]-[UA1]-[DA1]  (Class four)

#######################################################################################################
################################## Data manilupation ##################################################
## format for how to create all of the combinations of Gene's you would want    
#excelFile = createExcelFile(currentDirectory,"Gene_ID_Compare_Group7.xlsx")
#<place the data you want to get here>
#E.g. getData(<folder name>, <letter (0,1,2,3)>)
#for i in range(0,len(data1)):
#    <place the operations you want here (subtract, intersect)>
#                     ||
#    <place the operations you want here (subtract, intersect)>
#    data2File = pandas.DataFrame(<place the final data from above here>)
#    data2File.to_excel(excelFile,sheet_name=str(i))
#excelFile.save()

excelFile = createExcelFile(currentDirectory,"Gene_ID_Compare_Group1.xlsx")
data1 = getData('up_regulate',2,dataOfDataStucturesNames)
data2 = getData('up_regulate',3,dataOfDataStucturesNames)
for i in range(0,len(data1)):
    data12 = intersectGenes(data1[i],data2[i])
    data2File = pandas.DataFrame(data12)
    data2File.to_excel(excelFile,sheet_name=str(i))
excelFile.save()

excelFile = createExcelFile(currentDirectory,"Gene_ID_Compare_Group2.xlsx")
data1 = getData('down_regulate',2,dataOfDataStucturesNames)
data2 = getData('down_regulate',3,dataOfDataStucturesNames)
for i in range(0,len(data1)):
    data12 = intersectGenes(data1[i],data2[i])
    data2File = pandas.DataFrame(data12)
    data2File.to_excel(excelFile,sheet_name=str(i))
excelFile.save()

excelFile = createExcelFile(currentDirectory,"Gene_ID_Compare_Group3.xlsx")
data1 = getData('up_regulate',1,dataOfDataStucturesNames)
data2 = getData('up_regulate',2,dataOfDataStucturesNames)
for i in range(0,len(data1)):
    data12 = intersectGenes(data1[i],data2[i])
    data2File = pandas.DataFrame(data12)
    data2File.to_excel(excelFile,sheet_name=str(i))
excelFile.save()

#Group7= [UC1-UA1-UB1]           (Class one)
excelFile = createExcelFile(currentDirectory,"Gene_ID_Compare_Group7.xlsx")
data1 = getData('up_regulate',2,dataOfDataStucturesNames)
data2 = getData('up_regulate',0,dataOfDataStucturesNames)
data3 = getData('up_regulate',1,dataOfDataStucturesNames)
for i in range(0,len(data1)):
    data12 = subtractGenes(data1[i],data2[i])
    data123 = subtractGenes(data12,data3[i])
    data2File = pandas.DataFrame(data123)
    data2File.to_excel(excelFile,sheet_name=str(i))
excelFile.save()

excelFile = createExcelFile(currentDirectory,"Gene_ID_Compare_Group8.xlsx")
data1 = getData('down_regulate',2,dataOfDataStucturesNames)
data2 = getData('down_regulate',0,dataOfDataStucturesNames)
data3 = getData('down_regulate',1,dataOfDataStucturesNames)
for i in range(0,len(data1)):
    data12 = subtractGenes(data1[i],data2[i])
    data123 = subtractGenes(data12,data3[i])
    data2File = pandas.DataFrame(data123)
    data2File.to_excel(excelFile,sheet_name=str(i))
excelFile.save()

# Group3 = [UB1,UC1]     
# Group9 = [Group3-DA1-DB1]        (Class two)
excelFile = createExcelFile(currentDirectory,"Gene_ID_Compare_Group9.xlsx")
data1 = getData('up_regulate',1,dataOfDataStucturesNames)
data2 = getData('up_regulate',2,dataOfDataStucturesNames)
data3 = getData('down_regulate',0,dataOfDataStucturesNames)
data4 = getData('down_regulate',1,dataOfDataStucturesNames)
for i in range(0,len(data1)):
    data12 = intersectGenes(data1[i],data2[i])
    data123 = subtractGenes(data12,data3[i])
    data1234 = subtractGenes(data123,data4[i])
    data2File = pandas.DataFrame(data1234)
    data2File.to_excel(excelFile,sheet_name=str(i))
excelFile.save()

# Group10= [UB1,DC1]-[UA1]-[DA1] (Class three)
excelFile = createExcelFile(currentDirectory,"Gene_ID_Compare_Group10.xlsx")
data1 = getData('up_regulate',1,dataOfDataStucturesNames)
data2 = getData('down_regulate',2,dataOfDataStucturesNames)
data3 = getData('up_regulate',0,dataOfDataStucturesNames)
data4 = getData('down_regulate',0,dataOfDataStucturesNames)
for i in range(0,len(data1)):
    data12 = intersectGenes(data1[i],data2[i])
    data123 = subtractGenes(data12,data3[i])
    data1234 = subtractGenes(data123,data4[i])
    data2File = pandas.DataFrame(data1234)
    data2File.to_excel(excelFile,sheet_name=str(i))
excelFile.save()

# Group11=[DB1,UC1]-[UA1]-[DA1]  (Class four)
excelFile = createExcelFile(currentDirectory,"Gene_ID_Compare_Group11.xlsx")
data1 = getData('down_regulate',1,dataOfDataStucturesNames)
data2 = getData('up_regulate',2,dataOfDataStucturesNames)
data3 = getData('up_regulate',0,dataOfDataStucturesNames)
data4 = getData('down_regulate',0,dataOfDataStucturesNames)
for i in range(0,len(data1)):
    data12 = intersectGenes(data1[i],data2[i])
    data123 = subtractGenes(data12,data3[i])
    data1234 = subtractGenes(data123,data4[i])
    data2File = pandas.DataFrame(data1234)
    data2File.to_excel(excelFile,sheet_name=str(i))
excelFile.save()

##Summary of the workbook
##Comparison between worksheets

##############################################################################
###################Find highlighted rows in excel#############################

## COMPARING GROUP FILES WITH HIGHLIGHTED FILES
# TODO: change directory and read in highlighted file
currentDirectory = "/Volumes/WD1/DEGlist/Comparison files"
os.getcwd()
os.chdir(currentDirectory)
Name_newfile=[]
## Give new file name here
Name_newfile='Gene_ID in other study.xlsx'
highlight_workbook=xlrd.open_workbook(Name_newfile,'rb')
highlight_sheet=highlight_workbook.sheet_by_index(0)
highlight_ID=[]
for i in range(0,highlight_sheet.nrows):
    for j in range(0,highlight_sheet.ncols): 
        j=0 ## Read the first column
    highlight_ID.append(highlight_sheet.cell(i,j).value)

print(highlight_ID)

# read in group file
#comparisonData = [][]
import pandas
import re
import openpyxl
# from openpyxl import load_workbook


comparisonData = xlrd.open_workbook('comparison data.xlsx','rb')
# Name_groupfile=[]
new_groupfile=[]
#Example compare highlight files with group7 file 
Name_groupfile='Gene_ID_Compare_Group7.xlsx'
## Get number as k from the Name_groupfile
regex=re.compile(r'\d+')
k=int(regex.search(Name_groupfile).group(0))
read_groupfile = pandas.ExcelFile(Name_groupfile)
for j in range(0,len(read_groupfile.sheet_names)):
    read_sheet=pandas.read_excel(read_groupfile, index_col=None, usecols="A", sheet_name=j)
    read_sheet_list=read_sheet.values.tolist()
    for i in range(0,len(highlight_ID)):
        ID_value=highlight_ID[i]
        if ID_value in read_sheet:
        #if read_sheet.ismember(ID_value):
            #comparisionData[i][j]=j+1
            sheet_to_write=comparisonData.sheet_by_index(k-6)
            cell_to_write=sheet_to_write.cell[i][j]
            cell_to_write.value=j+1         
        else:
            sheet_to_write=comparisonData.sheet_by_index(k-6)
            # cell_to_write=sheet_to_write.cell[i][j]
            cell_to_write = 3
            cell_to_write.value='N'
comparisonData.save()



















