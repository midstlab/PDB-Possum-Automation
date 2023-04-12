from genericpath import exists
import xlsxwriter
import os
def txttoexcel(destination):
    #print(destination)
    arr = os.listdir(destination) #current directory (directory that txt files found)
    strtxt = ".txt"
    finaldest = destination + "/ExcelFiles"
    os.mkdir(finaldest)
    for txtfile in arr:
        if txtfile.__contains__(strtxt): # to take only txt files 
            name = txtfile[:4]
            fileObject = open(destination + "/" +txtfile, "r")
            if not exists (finaldest + "/"+ name + ".xlsx" ): #if it exist don't edit
                with open(finaldest + "/"+ name + ".xlsx", "a+") as c:
                    print("{}.txt file converted to xlsx file!".format(name))
                #print(name, ": ", txtfile) #debug purposes 
                workbook = xlsxwriter.Workbook(finaldest + "/"+ name + ".xlsx") #create a workbook
                worksheet = workbook.add_worksheet("structure") #create a worksheet named ‘structure’
                header ="structure1"
                worksheet.set_header(header)
                #data = open("/Users/ahmetolcum/Documents/Sabanci/ENS491/9/Report_PoSSuM (2).txt","r") #load data
                linelist = fileObject.readlines() #read each line in the txt
                count = len(linelist) #count the number of lines
                #print (count) #print the number of lines
                a=0
                linecount = 0 
                newrow=[]
                pdbid, Hetcode, uniprot = "","",""
                newline = ""
                for num in range(0, count): #create each line 
                    line = linelist[num] #load each line in variable line
                    if(line.find(":") != -1):
                        newline = line.split(":")
                    else:
                        newline = line.split(" ")
                    if(line[0] == "#"):
                        splitline = line.split("|") #split lines
                        splitline = splitline[1:]
                        worksheet.write_row(a, 0, splitline) #write each line in excel “{PDBid}.xlsx”
                        a+=1
                        if(a == 1):
                            for i in range(len(splitline)):
                                newrow.append("")
                            newrow[0] = pdbid
                            newrow[1] = Hetcode
                            newrow[9] = uniprot
                            worksheet.write_row(a,0,newrow)
                            a+=1
                        #print(splitline)
                    elif(newline[0]  == "PDB ID"):
                        pdbid = newline[1].strip("\n").strip(" ")
                    elif(newline[0] == "HET code"):
                        Hetcode = newline[1].strip("\n").strip(" ")
                    elif(newline[0] == "UniProt ID"):
                        uniprot = newline[1].strip("\n").strip(" ")
                    elif(line[0] == "-"):
                        linecount += 1 

                workbook.close() #<<<close workbook, important to complete the task
    print("Directory is completed") # feedback about last completed directory that have multiple txt files in it (named in grouplist)
    return finaldest
if __name__ == "__main__":
    destination = "/Users/ahmetolcum/Downloads/PDBSUM/COA"
    txttoexcel(destination)

