import pdb
import Bio
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
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
#aynı proteinin farklı chainleri kıyaslamayı yapalım
def findsequni(id,chain):
  seq=""
  if id!="" or id!="none":
    try:
        #https://www.uniprot.org/uniprot/P31572.fasta
        url = "https://www.uniprot.org/uniprot/"+id+".fasta"
        
        context = ssl._create_unverified_context()
        html = urlopen(url, context=context)
        soup = BeautifulSoup(html, 'lxml')
        type(soup)
        rows = soup.find('p')
        #print(rows.string)
        fasta = rows.string 
        fas_seq = fasta.split('\n')
        #print(len(fas_seq))
        #print(fas_seq)
        for k in range(1,len(fas_seq)-1):
            seq = seq + str(fas_seq[k])
    except:
        return seq
    return seq
def findseq(id,chain):  
  seq=""
  url = "https://www.rcsb.org/fasta/entry/" + id + "/display"
  context = ssl._create_unverified_context()
  html = urlopen(url, context=context)
  soup = BeautifulSoup(html, 'lxml')
  type(soup)
  rows = soup.find('p')
  #print(rows.string)
  fasta = rows.string 
  fas_seq = fasta.split('\n')
  if len(fas_seq)>2:
    i=0
    k=0
    isfound = False
    while i < len(fas_seq) and isfound==False:
      header = fas_seq[i].split('|')
      if header!=['']:
        chain_info=header[1]
        #chain_info=chain_info.replace(',', '')
        chain_lst = chain_info.split(',')
        #chain_lst = chain_lst[1:]
        #print(chain_lst)
        if len(chain_lst)==1 and chain== "Chain " + chain:
          k=i
        else:
          for item in chain_lst:
            if item == chain_lst[0]:
              item = item.replace("Chains","")
              idx = item.find("auth")
              if idx!=-1:
                new_item = item.split('[')
                #print(new_item)
                for it in new_item:
                  if it == "auth " + str(chain) + "]":
                    k=i
                    #print("here")
                    isfound=True
                  elif it==chain:
                    k=i
                  elif item != chain_lst[0]:
                    idx = item.find("auth")
                    #print(idx)
                    if idx!=-1:
                      new_item = item.split('[')
                      #print(new_item)
                      for it in new_item:
                        if it== "auth " + str(chain) + "]":
                          k=i
                          isfound=True
                        elif it==chain:
                          k=i
                    elif item == " " + str(chain):
                      k=i
      i+=2
    seq= fas_seq[k+1]
  else:
    seq= fas_seq[1]
  return seq
def compare(i,k,df,pdbID,uniprotID,chains):
    to_keep = 0  
    chain_i = chains[i]
    chain_k = chains[k]
    ids_i= uniprotID[i].split(" ")
    ids_k= uniprotID[k].split(" ")
    seq_i=""
    seq_k=""
    if uniprotID[i]!="none" and uniprotID[k]!="none" and len(ids_i)==1 and len(ids_k)==1:
        seq_i = findsequni(uniprotID[i],chain_i)
        seq_k = findsequni(uniprotID[k],chain_k)
    else:
        seq_i = findseq(pdbID[i],chain_i)
        seq_k = findseq(pdbID[k],chain_k)
    if seq_i==seq_k:
        rmsd_i = df['RMSD(Ca)'][i + 1]
        rmsd_k = df['RMSD(Ca)'][k + 1]
        if rmsd_i==rmsd_k:
            al_i = df['Aligned length'][i + 1]
            al_k = df['Aligned length'][k + 1]
            if al_i==al_k:
                cos_i = df['Cosine value'][i + 1]
                cos_k = df['Cosine value'][k + 1]
                if cos_i==cos_k:
                    p_i = df['p value'][i + 1]
                    p_k = df['p value'][k + 1]
                    minp = min(p_i,p_k)
                    if minp==p_i:
                        to_keep = i + 1
                    else:
                        to_keep = k + 1
                else:
                    maxcos= max(cos_i,cos_k)
                    if maxcos == cos_i:
                        to_keep=i + 1
                    else:
                        to_keep=k + 1
            else:
                maxal= max(al_i,al_k)
                if maxal == al_i:
                    to_keep=i+1
                else:
                    to_keep=k + 1
        else:
            minrmsd = min(rmsd_i,rmsd_k)
            if minrmsd==rmsd_i:
                to_keep=i + 1
            else:
                to_keep=k + 1
    else:
        to_keep=-1
    
    return to_keep

def compareuni(i,k,df0):
    to_keep = 0
    het_i = df0['HET code'][i]
    het_k =df0['HET code'][k]
    if(het_i!=het_k):
        return -1
    rmsd_i = df0['RMSD(Ca)'][i]
    rmsd_k = df0['RMSD(Ca)'][k]
    if rmsd_i==rmsd_k:
        al_i = df0['Aligned length'][i]
        al_k = df0['Aligned length'][k]
        if al_i==al_k:
            cos_i = df0['Cosine value'][i]
            cos_k = df0['Cosine value'][k]
            if cos_i==cos_k:
                p_i = df0['p value'][i]
                p_k = df0['p value'][k]
                minp = min(p_i,p_k)
                if minp==p_i:
                    to_keep = i 
                else:
                    to_keep = k 
            else:
                maxcos= max(cos_i,cos_k)
                if maxcos == cos_i:
                    to_keep=i
                else:
                    to_keep=k
        else:
            maxal= max(al_i,al_k)
            if maxal == al_i:
                to_keep=i
            else:
                to_keep=k
    else:
        minrmsd = min(rmsd_i,rmsd_k)
        if minrmsd==rmsd_i:
            to_keep=i
        else:
            to_keep=k
    return to_keep

def find_occurence(element,pdbID):
    occurences=[]
    for i in range(len(pdbID)):
        if pdbID[i]==element:
            occurences.append(i)
    return occurences

def post_possum(filename,folderpath,ligand,clean):
    dest = folderpath+"/ResultFiles"
    arr = [f for f in os.listdir(folderpath) if not f.startswith('.')]#os.listdir(folderpath)
    #print(arr)
    folderpath = folderpath + "/"
    #os.mkdir(dest)
    nonputs = []
    #print(ligand)
    #for filename in arr:
    try:
        print("start " + filename + "\n")
        df1 = pd.read_excel(join(folderpath, filename))
        #df.head(32)
        names = ['PDB ID','HET code','Chain ID','Res. No.','Cosine value','p value','Aligned length','RMSD(Ca)','Protein Name','UniProt ID','UniRef50','EC No.','CATH code','SCOPe code','Aligned residues (Ca atoms)']
        df = df1.drop_duplicates(subset = ['PDB ID', 'Chain ID']).reset_index(drop = True)
        pdbID=df['PDB ID'].tolist()
        #print(pdbID)
        chains = df['Chain ID'].tolist()
        uniprotID= df['UniProt ID'].tolist()
        uniprotID.pop(0)
        chains.pop(0)
        chains_to_keep =  {}
        to_keep= []
        df_list= df.values.tolist()
        to_keep.append(df_list[0])
        print(df_list[0])
        seen = []
        pdbID.pop(0)
        for element in pdbID:
            if element not in seen:
                print(element)
                occurences = find_occurence(element,pdbID)
                print(len(occurences))
                copyocc = find_occurence(element,pdbID)
                if len(occurences)==1:
                    to_keep.append(df_list[occurences[0]+1])
                else:
                    for a in range(0,len(occurences)):
                        for b in range(a+1,len(occurences)):
                            idx=compare(occurences[a],occurences[b],df,pdbID,uniprotID,chains)
                            if idx==occurences[a]:
                                if occurences[b] in copyocc:
                                    copyocc.remove(occurences[b])
                            elif idx==occurences[b]:
                                if occurences[a] in copyocc:
                                    copyocc.remove(occurences[a])
                    for x in copyocc:
                        if df_list[x] not in to_keep:
                            to_keep.append(df_list[x])
                seen.append(element)
            else:
                print("already saw "+str(element))

        df0 = pd.DataFrame(to_keep)
        names = ['PDB ID','HET code','Chain ID','Res. No.','Cosine value','p value','Aligned length','RMSD(Ca)','Protein Name','UniProt ID','UniRef50','EC No.','CATH code','SCOPe code','Aligned residues (Ca atoms)']
        df0.columns=names
        to_keep2 = []
        df_list2= df0.values.tolist()
        #print(df_list2[0])
        #to_keep2.append(df_list2[0])
        seen2 = []
        uniprotID=df0['UniProt ID'].tolist()
        
        for element in uniprotID:
            #print(element)
            if element not in seen2:
                occurences = find_occurence(element,uniprotID)
                #print(len(occurences))
                copyocc = find_occurence(element,uniprotID)
                if len(occurences)==1:
                    to_keep2.append(df_list2[occurences[0]])
                else:
                    for a in range(0,len(occurences)):
                        for b in range(a+1,len(occurences)):
                            idx=compareuni(occurences[a],occurences[b],df0)
                            #print("idx" + str(idx))
                            if idx==occurences[a]:
                                if occurences[b] in copyocc:
                                    copyocc.remove(occurences[b])
                            elif idx==occurences[b]:
                                if occurences[a] in copyocc:
                                    copyocc.remove(occurences[a])
                    for x in copyocc:
                        if df_list2[x] not in to_keep2:
                            to_keep2.append(df_list2[x])
                seen2.append(element)
            else:
                print("already saw "+str(element))
        #print(to_keep2)
        
        df_final = pd.DataFrame(to_keep2)
        
        df_final.columns=names
        temp=[]
        if(clean):
            lst = df_final.values.tolist()
            for k in lst:
                lig = k[1].strip(' \n')
                #print(lig)
                if lig in ligand:
                    temp.append(k)
        if(len(temp)>1):
            df_final = pd.DataFrame(temp)
            df_final.columns=names 
        

        df_final.to_excel(join(dest, filename))
        print("end " + filename + "\n")
    except:
        nonputs.append(filename)
        print("Unexpected problem, while working with " + filename + ".xlsx!")
    return dest