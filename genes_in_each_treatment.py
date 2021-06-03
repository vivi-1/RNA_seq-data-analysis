import pandas
import os

print(os.getcwd())
os.chdir('/Users/weiwang/Desktop/DEGlist/Output/Eliminate outlier rep2')
print(os.getcwd())

def find_gene(gene_ID):
    for i in range(1,8):
        print("timepoint ", i)
        down_or_up = ["down_", "up_"]
        treatment = ["Flag22_", "Pnic_", "Flag22_Pnic_"]
        for j in down_or_up:
            for n in treatment:
                fileName = "T" + str(i) + "_" + j + n +"no rep2.csv"
                file = pandas.read_csv(fileName)
                pool = file['Gene_ID']
                if gene_ID in pool.values:
                    print(gene_ID+ " timepoint"+str(i)+" "+n +" "+j+ "regulated")

find_gene("Niben101Scf00528g00001")
#timepoint1 Flag22_Pnic_ up_regulated

find_gene("Niben101Scf01738g02011")
find_gene("Niben101Scf01870g00003")
find_gene("Niben101Scf02348g10009")

find_gene("Niben101Scf03380g03005")
#timepoint2 Pnic_ up_regulated
#timepoint2 Flag22_Pnic_ up_regulated

find_gene("Niben101Scf06087g02007")
find_gene("Niben101Scf12330g02013")
