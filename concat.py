import pandas as pd
from os.path import join
import os
#RMSD 1.5 dan küçük olanlar alınacak geri kalanlar gruplanmak için bekleyecek
from sympy import false
def concater(prot, prot2, path, everexisted, opnum):
            if prot2.__contains__(".xlsx"):
                if prot != prot2 and ((prot.removeprefix(".xlsx").__contains__(prot2.removeprefix(".xlsx")) == false) and ((prot2.removeprefix(".xlsx").__contains__(prot.removeprefix(".xlsx")) == false))):
                    df1= pd.read_excel(join(path, prot))
                    df2= pd.read_excel(join(path, prot2))
                    pdbID1 = df1.iloc[:, 0].to_list()
                    pdbID2 = df2.iloc[:, 0].to_list()
                    exist = False
                    for element in pdbID1:
                        if element in pdbID2:
                            exist=True
                            break
                    if exist==True:
                        everexisted = True
                        df1=pd.concat([df1,df2])
                        newdf1 = df1.drop_duplicates(['PDB ID', 'Chain ID'])
                        filename = prot.removesuffix(".xlsx") + "_" + prot2
                        newdf1.to_excel(join(path,filename),index=False)
                        os.remove(join(path,prot))
                        os.remove(join(path,prot2))
                        opnum +=1
            return everexisted, opnum
def folderconcat(path_prefix):
    protein_list= os.listdir(path_prefix)
    #print(protein_list)
    opnum = 4
    while(opnum != 0):
        opnum = 0
        protein_list= os.listdir(path_prefix)
        for prot in protein_list:
            everexisted = False
            try:
                protein_list= os.listdir(path_prefix)
                if prot.__contains__(".xlsx"):
                    for prot2 in protein_list:
                        if prot2.__contains__(".xlsx"):
                            everexisted, opnum = concater(prot, prot2, path_prefix, everexisted, opnum)
                if everexisted:
                    print("modification is done")
                    protein_list= os.listdir(path_prefix)
            except:
                print(prot + " already merged!")
    print("Merge is done for, " + path_prefix)
    return path_prefix
"""
if __name__ == "__main__":
    folderconcat("/Users/ahmetolcum/Documents/Sabanci/ENS491/PoSSum/try")
"""