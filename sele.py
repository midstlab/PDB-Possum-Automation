##import libraries 
#from concurrent.futures import thread
from multiprocessing import Pool
import os
from selenium import webdriver
import time
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
##import files
from txtconverter import txttoexcel
from concat import folderconcat
from post_possum import post_possum
from grouping import grouping
from align import align
from query_pdb import caller

def matrixchecker(High_RMSDs, protein):
    exist = False
    for element in High_RMSDs:
        if protein == element[1]:
            exist = True
    return exist
def eliminatorhelper(prodict, protein, eliminates):
    for key in prodict:
        if protein in prodict[key]:
            prodict[key].remove(protein)
            eliminates.append(protein)
    return eliminates
def allelementsindictionary(prodict):
    numberofelements = 0
    for key in prodict:
        numberofelements += len(prodict[key])
    return numberofelements
def chromesetter(chromeDriver ,dest, headless):
    WINDOW_SIZE = "1920,1080"
    chromeOptions = Options()
    chromeOptions.add_experimental_option('excludeSwitches', ['enable-logging'])
    prefs = {"download.default_directory" : dest}
    if(headless):
        chromeOptions.add_argument("--headless")
        chromeOptions.add_argument("--window-size=%s" % WINDOW_SIZE)
    chromeOptions.add_experimental_option("prefs",prefs)
    browser = webdriver.Chrome(executable_path=chromeDriver, chrome_options=chromeOptions)
    return browser
def inputtaker():
    inpligands = []
    max_RMSD = float(0.0)
    dest1 = input("Enter a path to save required files (Leave empty for current directory): ")
    print("Your input may be a list of ligands or a single one! For multiple inputs put \", \" as a separator!")
    inpligands = input("Enter your ligands' letter code: ")
    while(inpligands == "" ):
        print("Wrong input format!")
        inpligands = input("Enter your ligand's three letter code: ")
    if(inpligands.find(",")):
        inpligands = inpligands.split(",")
        for element in inpligands:
            if element.find(" "):
                element = element.strip()
    for element in inpligands: 
        try:
            os.mkdir(dest1 + element)
        except:
            print(element+ " folder already existed and will be overwriten!")
    if dest1 != "":
        if dest1[len(dest1) - 1 != "/"]:
            dest1 = dest1 + "/"
    else:
        dest1 = os.getcwd()
    max_RMSD = input("Please enter maximum RMSD value (Leaving empty means no any limit): ")
    if(max_RMSD == ""):
        max_RMSD = "10"
    correctinp = 1
    while(correctinp <= 1):
        cleaninp = input("Do you want to have other ligands in your result? (Whether you want to see other items than yours)(Yes : \"y\", No : \"N\"): ")
        if cleaninp == "y":
            clean = False
            correctinp += 1
        elif(cleaninp == "N"):
            clean = True
            correctinp += 1
        else:
            print("Wrong Input!!!")
            correctinp = 0
    
    return inpligands,  dest1, clean, max_RMSD
def download_wait(path_to_downloads):
    seconds = 0
    dl_wait = True
    while dl_wait and seconds < 20:
        time.sleep(1)
        dl_wait = False
        for fname in os.listdir(path_to_downloads):
            if fname.endswith('.crdownload'):
                dl_wait = True
        seconds += 1
    return seconds
def PDBconnecter(dest, ligand, chromedriver, headless, prodict):
    liste = []
    numtotal =0
    for element in ligand:
        browser1 = chromesetter(chromedriver, dest, headless)
        browser1.get("https://www.rcsb.org")
        time.sleep(2)
        search = browser1.find_element(By.XPATH,"/html/body/div/div[2]/div/div/div[1]/div[2]/div/div[2]/table/tbody/tr[1]/td[2]/div/input")
        search.send_keys(element)
        time.sleep(1)
        find = browser1.find_element(By.XPATH,"/html/body/div/div[2]/div/div/div[1]/div[2]/div/div[2]/table/tbody/tr[1]/td[2]/div/div/div[2]")
        find.click()
        time.sleep(1)
        but = browser1.find_element(By.XPATH,"/html/body/div/div[3]/div[1]/div[3]/div[2]/input")
        but.click()
        time.sleep(1)
        down = browser1.find_element(By.XPATH,"/html/body/div/div[3]/div/div/div[3]/div[2]/div[3]/table/tbody/tr/td[2]/div/div[1]/div[3]/div/div[2]/span")
        down.click()
        time.sleep(1)
        try:
            browser1.switch_to.alert.dismiss()
        except:
            #fallback procedure if there is no any pop-up
            print("\n")
        time.sleep(2)
        listprot = browser1.find_element(By.TAG_NAME,"textarea")
        liste = []
        if len(liste) == 0:
            listetext = listprot.text
            if(listetext.find(",") == -1):
                liste.append(listetext)
            else:
                liste = listetext.split(",")
            for val in prodict.values():
                for item in val:
                    if item in liste:
                        liste.remove(item)
            numtotal += len(liste)
            prodict[element] = liste
    browser1.close()
    temp = []
    return prodict
def eliminator(destination, element, key, max_RMSD, prodict):
    fileObject = open("{}/{}.txt".format(destination, element), "r")
    lines = fileObject.readlines()
    newlines = []
    exist = False
    somethingtodelete = False
    High_RMSDs = []
    eliminates = []
    for line in lines:
        if (line[0] == "#"):
            row = line.split("|")
            if(len(row[1]) == 4):
                protein = row[1]
                eliminates = eliminatorhelper(prodict, protein, eliminates)
            if(len(row[8]) <= 4):
                rmsd_value = float(row[8])
                if(float(rmsd_value) > float(max_RMSD)):
                    exist = matrixchecker(High_RMSDs, row[1])
                    if (not exist):
                        High_RMSDs.append(row)
                    somethingtodelete = True
                else:
                    newlines.append(line)
            else:
                newlines.append(line)
        else:
            newlines.append(line)
    fileObject.close()
    if(somethingtodelete):
        fileObject2 = open("{}/{}.txt".format(destination, element), "w")
        for line in newlines:
            fileObject2.write(line)
    return High_RMSDs
def possumdownloader(ligand, prodict, destination, chromedriver, headless, max_RMSD):
    noresultlist = []
    SomeErrorList = []
    destroot = destination
    High_RMSDs, High_RMSDs2 = [],[]
    if prodict:
        time.sleep(1)
        done = 0 
        a = 0 
        toc = time.time()
        eliminates = []
        print("Estimation is calculating...")
        for key in prodict:
                arr = [f for f in os.listdir(destination + "/" + key) if not f.startswith('.')]#os.listdir(folderpath)
                for element in prodict[key]:
                    if(element+".txt" not in arr):
                        #print(element)
                        if len(element) == 4 :
                            try:
                                destination = destroot + "/" + key
                                driver = chromesetter(chromedriver, destination, headless)
                                driver.get("http://possum.cbrc.jp/PoSSuM/search_k.html")
                                time.sleep(1)
                                pdbid = driver.find_element(By.NAME,"params[0]")
                                ligand_loc = driver.find_element(By.NAME,"params[1]")
                                pdbid.send_keys(element)
                                #print(key)
                                ligand_loc.send_keys(key)
                                download = driver.find_element(By.XPATH,"/html/body/div/div[2]/div[2]/form/table/tbody/tr[11]/td[2]/input[2]")
                                download.click()
                                time.sleep(0.5)
                                submit = driver.find_element(By.NAME,"button")
                                submit.click()
                                time.sleep(0.5)
                                done += 1
                                namenotchanged = True
                                secs = download_wait(destination)
                                if namenotchanged:
                                    try:
                                        if (destination +"/Report_PoSSuM.txt"):
                                            filename = destination + "/" + element + ".txt"
                                            os.rename(destination +"/Report_PoSSuM.txt", filename)
                                            namenotchanged = False
                                    except:
                                        noresultlist.append(element)
                                        print("No any result for the " + element + "!")
                                if not namenotchanged and max_RMSD != None:
                                    High_RMSDs += eliminator(destination, element, key, max_RMSD, prodict)
                                #driver.close()
                                if (len(High_RMSDs) != 0):
                                    #print(High_RMSDs)
                                    for item in High_RMSDs:
                                        #print(item)
                                        destination = destroot + "/" + key
                                        driver = chromesetter(chromedriver, destination, headless)
                                        driver.get("http://possum.cbrc.jp/PoSSuM/search_k.html")
                                        time.sleep(1)
                                        pdbid = driver.find_element(By.NAME,"params[0]")
                                        ligand_loc = driver.find_element(By.NAME,"params[1]")
                                        pdbid.send_keys(item[1])
                                        ligand_loc.send_keys(item[2])
                                        download = driver.find_element(By.XPATH,"/html/body/div/div[2]/div[2]/form/table/tbody/tr[11]/td[2]/input[2]")
                                        download.click()
                                        time.sleep(0.5)
                                        submit = driver.find_element(By.NAME,"button")
                                        submit.click()
                                        time.sleep(0.5)
                                        namenotchanged = True
                                        secs = download_wait(destination)
                                        try:
                                            if (destination +"/Report_PoSSuM.txt"):
                                                #print(item[1])
                                                filename = destination + "/" + item[1] + ".txt"
                                                os.rename(destination +"/Report_PoSSuM.txt", filename)
                                                namenotchanged = False
                                        except:
                                            noresultlist.append(item)
                                            print("No any result for the " + item[1] + "!")
                                        if not namenotchanged and max_RMSD != None:
                                            High_RMSDs2 = eliminator(destination, item[1], key, max_RMSD, prodict)
                            except:
                                SomeErrorList.append(key+ ":" +element)
                    else:
                        done+=1    
                    High_RMSDs = []
                    eliminates += High_RMSDs2
                    numelement = allelementsindictionary(prodict)
                    #print(numelement)
                    percent = done * 100
                    percent = float(percent) / float(numelement)
                    print("{:.2f}".format(percent) ,"% completed! Download in progress...", )
                    if(percent > a):
                        a = a + 5
                        tic = time.time()
                        remtime = tic - toc
                        estimate = (100-percent) * remtime
                        estimate = estimate / percent
                        estimate = estimate / 60
                        print(f"Estimated time to complete {estimate:0.4f} minutes")
        print("Download completed!")
        print(SomeErrorList)         
    driver.close()
    return noresultlist, High_RMSDs2, destination
def main():
    #The desired directory needed to be changed. Destination can be found in inputtaker function 
    chromedriver = "/Users/ahmetolcum/Documents/DevSpace/Selenium/chromedriver" #Give the path of the chrome driver for the selenium 
    prodict = {} 
    resultless, lenlist = [], []
    clean = False
    ligand, destination, clean, max_RMSD = inputtaker()
    headless = True
    prodict = caller(ligand, prodict)
    ###For using selenium not PDB API 
    ###prodict = PDBconnecter(destination, ligand, chromedriver, headless, prodict)
    resultless, High_RMSDs, destination = possumdownloader(ligand, prodict, destination, chromedriver, headless, max_RMSD)
    destination = txttoexcel(destination)
    dest = destination+"/ResultFiles"
    os.mkdir(dest)
    arr = [f for f in os.listdir(destination) if ((not f.startswith('.')) & (f.endswith(".xlsx")))]
    lst = []
    for i in arr:
        lst.append((i,destination,ligand,clean))
    p = Pool(6)
    destination= p.starmap(post_possum,lst)
    destination = destination[0]
    dest = destination+"/AlignedResults"
    
    os.mkdir(dest)
    '''except:
        print("Folder exists will be overwritten: Aligned Result Files!")'''
    arr = [f for f in os.listdir(destination) if ((not f.startswith('.')) & (f.endswith(".xlsx")))]
    lst = []
    for i in arr:
        lst.append((i,destination))
    p = Pool(6)
    destination= p.starmap(align,lst)
if __name__ == "__main__":
    main()

############PROBLEMS TO FIX##############
"""

"""
