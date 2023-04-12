from cmath import nan
import pandas as pd
import numpy as np
from os.path import join
import requests
#%matplotlib inline
from urllib.request import Request, urlopen
from bs4 import BeautifulSoup
import lxml.html as lh
import ssl
import os
def findseq(id,chain):
  seq=""
  #https://www.uniprot.org/uniprot/P31572.fasta
  url = "https://www.uniprot.org/uniprot/"+id+".fasta"
  
  context = ssl._create_unverified_context()
  html = urlopen(url, context=context)
  soup = BeautifulSoup(html, 'lxml')
  type(soup)
  rows = soup.find('p')
  #print(rows.string)
  print(id) #protein uniprotta yoksa ne yapalım alignment konusunda?
  #uniprotid not found tarzı bir uyarı 
  if rows==None:
    seq="problem"
    return seq
  fasta = rows.string 
  fas_seq = fasta.split('\n')
  #print(len(fas_seq))
  #print(fas_seq)
  for k in range(1,len(fas_seq)-1):
    seq = seq + str(fas_seq[k])
  return seq

def align(item,folderpath):
    dest = folderpath+"/AlignedResults"
    #arr = [f for f in os.listdir(folderpath) if ((not f.startswith('.')) & (f.endswith(".xlsx")))]#os.listdir(folderpath)
    notfound = []
    '''try:
        os.mkdir(dest)
    except:
        print("Folder exists will be overwritten: Aligned Result Files!")'''
    
    #for item in arr:
    bind_sites=[]
    print(item)
    #try:
        #print(item)
    #item="1XVT.xlsx" # dosya adını değiştir
    df1= pd.read_excel(join(folderpath, item)) #açılacak dosya adı 
    uniprotID= df1['UniProt ID'].tolist()
    print(uniprotID[0])
    orig_protein = findseq(uniprotID[0],"A")
    if orig_protein=="problem":
        notfound.append(item)
    else:
        bind_sites.append(orig_protein)
        #print(orig_protein)
        pdbID=df1['PDB ID'].tolist()
        #print(len(pdbID))
        res = df1['Aligned residues (Ca atoms)'].tolist()
        for k in range(len(pdbID)):
            stri = ""
            if pd.isna(res[k]):
                bind_sites.append(stri)
            else:
                for i in range(len(orig_protein)):
                    stri = stri + "-"
                    lst1=res[k].split(",")
                for m in range(len(lst1)):
                    str1 = lst1[m]
                    idx = str1.rfind("_")
                    aa = str1[idx+1]
                    idx_ = str1.find("_")
                    idxline = str1.find("-")
                    loc = int(str1[idx_+2:idxline])
                    #print(loc)
                    stri = stri[:loc-1] + aa + stri[loc:]
                bind_sites.append(stri)
        #print(bind_sites)
        df1['Binding Site Alignment']=bind_sites[1:]
        w = open(dest + "/" + item.strip(".xlsx") + ".txt", "w",  encoding='utf-8')
        w.write(item.strip(".xlsx")+"\t"+orig_protein+"\n")
        for key in range(len(bind_sites[1:])):
            w.write(pdbID[key]+"\t"+bind_sites[key+1]+"\n")
        w.close()
    
        #except:
            #print("Encountered with a problem in alignment for: ", item)
    if len(notfound)>0:
        print("Uniprot ID could not be found for: \n")
        for id in notfound:
            print(id)
            print('\n')
    return dest
