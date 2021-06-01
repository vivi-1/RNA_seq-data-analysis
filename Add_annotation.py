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

annotation_path="/Users/weiwang/Desktop/DEGlist/Input/Niben101_annotator.xlsx"
annotation_file=pandas.read_excel(annotation_path, sheet_name=0)
print(annotation_file)
print(annotation_file.Gene_ID)

def annotation(classfilePath, currentDirectory, classNumber):
    excelFile = createExcelFile(currentDirectory,"Class"+classNumber+"_Gene_ID_annotation.xlsx")

    for i in range(1,8):
        print("i")
        print(i)
        classFile = pandas.read_excel(classfilePath, sheet_name=str(i)).Gene_ID
        annotation = []
        for j in range(0,len(classFile)):
            annotation_to_write_row_number=annotation_file.Gene_ID[annotation_file.Gene_ID == classFile[j]].index.tolist()
            if annotation_to_write_row_number==[]:
                annotation.append("N")

            else:
                annotation_to_write_row=int(str(annotation_to_write_row_number).strip('[]'))
                annotation_value=annotation_file.loc[annotation_to_write_row].at['Unnamed: 3']
                annotation.append(annotation_value)
        DataFrame({"Annotation":annotation}).to_excel(excelFile,sheet_name=str(i))
    excelFile.save()

annotation('/Users/weiwang/Desktop/DEGlist/Output/Eliminate outlier rep2/Class2_2_Gene_ID.xlsx', '/Users/weiwang/Desktop/DEGlist/Output/Eliminate outlier rep2', "2_2")
