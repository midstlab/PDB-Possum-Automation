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
def compare(i,k,df,pdbID,chains):
    to_keep = 0
    chain_i = chains[i]
    chain_k = chains[k]
    seq_i = findseq(pdbID[i],chain_i)
    seq_k = findseq(pdbID[k],chain_k)
    if seq_i==seq_k:
        rmsd_i = df['RMSD(Ca)'][i]
        rmsd_k = df['RMSD(Ca)'][k]
        if rmsd_i==rmsd_k:
            al_i = df['Aligned length'][i]
            al_k = df['Aligned length'][k]
            if al_i==al_k:
                cos_i = df['Cosine value'][i]
                cos_k = df['Cosine value'][k]
                if cos_i==cos_k:
                    p_i = df['p value'][i]
                    p_k = df['p value'][k]
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
    else:
        to_keep=-1
    
    return to_keep
def find_occurence(element,pdbID):
    occurences=[]
    for i in range(len(pdbID)):
        if pdbID[i]==element:
            occurences.append(i)
    return occurences           
def post_possum(folderpath):
  dest = folderpath+"/resultfiles"
  arr = [f for f in os.listdir(folderpath) if not f.startswith('.')]#os.listdir(folderpath)
  os.mkdir(dest)
  for filename in arr:
    print(filename)
    df = pd.read_excel(folderpath +"/" + filename)
    #df = df1.drop_duplicates(['PDB ID', 'Chain ID'])
    pdbID=df['PDB ID'].tolist()
    chains = df['Chain ID'].tolist()
    chains_to_keep =  {}
    to_keep= []
    queryname = filename.removesuffix('.xlsx')
    query = [queryname,'','','','','','','','','','','','','','' ]
    to_keep.append(query)
    df_list= df.values.tolist()
    processed=[]
    for element in pdbID:
        print(element)
        if element not in processed:
          occurences = find_occurence(element,pdbID)
          print(len(occurences))
          copyocc = find_occurence(element,pdbID)
          if len(occurences)==1:
              to_keep.append(df_list[occurences[0]])
          else:
              for a in range(0,len(occurences)):
                  for b in range(a+1,len(occurences)):
                      idx=compare(occurences[a],occurences[b],df,pdbID,chains)
                      if idx==occurences[a]:
                          if occurences[b] in copyocc:
                              copyocc.remove(occurences[b])
                      elif idx==occurences[b]:
                          if occurences[a] in copyocc:
                              copyocc.remove(occurences[a])
              print(len(copyocc))
              for x in copyocc:
                  if df_list[x] not in to_keep:
                      print(df_list[x])
                      to_keep.append(df_list[x])
          processed.append(element)
    df0 = pd.DataFrame(to_keep)
    names = ['PDB ID','HET code','Chain ID','Res. No.','Cosine value','p value','Aligned length','RMSD(Ca)','Protein Name','UniProt ID','UniRef50','EC No.','CATH code','SCOPe code' , 'Aligned residues (Ca atoms)']
    df0.columns=names
    pdbid=df0['PDB ID'].tolist()
    res = df0['Aligned residues (Ca atoms)'].tolist()
    bind_sites = []
    search_seq=""
    try:
      lst1=res[1].split(",")
      for k in range(len(lst1)):
        str1 = str(lst1[k])
        if(str1 != ""):
          idx = str1.find("_")
          search_seq = search_seq+str1[idx+1]
        else:
          print("site sequence resulted empty! Something is wrong!")
      if (search_seq != ""):
        bind_sites.append(search_seq)
      else:
        print(str1)
      for i in range(len(res)):
          lst1=res[i].split(",")
          site_seq=""
          for k in range(len(lst1)):
              str1 = str(lst1[k])
              if(str1 != ""):
                idx = str1.rfind("_")
                site_seq = site_seq+str1[idx+1]
              else:
                print("one of the other sequence resulted empty! Something is wrong!")
          if site_seq != "":
            bind_sites.append(site_seq)
          else:
            print(str1)
      print(df0)
      df0['Aligned residues (Ca atoms)']=bind_sites
      df0.to_excel(join(dest,filename))
    except:
      print("An error occured during the post possum!")
  return dest

