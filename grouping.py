#%pip install biopython
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

def grouping(folderpath):
    dest = folderpath+"/groupfiles"
    arr = [f for f in os.listdir(folderpath) if not f.startswith('.')] #os.listdir(folderpath)
    os.mkdir(dest)
    for filename in arr:
        df_groups = pd.read_excel(folderpath + "/" + filename)
        pdbID = df_groups.iloc[:, 1].to_list()
        pdbID.pop(0)
        print(pdbID)
        chains = df_groups.iloc[:, 3].to_list()
        chains.pop(0)
        print(chains)
        seq_dict = {}
        for id in pdbID:
            print(id)
            url = "https://www.rcsb.org/fasta/entry/" + id + "/display"
            context = ssl._create_unverified_context()
            html = urlopen(url, context=context)
            soup = BeautifulSoup(html, 'lxml')
            type(soup)
            rows = soup.find('p')
            #print(rows.string)
            fasta = rows.string 
            fas_seq = fasta.split('\n')
            #print(fas_seq[1])
            if len(fas_seq)>2:
                chain = chains[pdbID.index(id)]
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
                                    #print(idx)
                                    if idx!=-1:
                                        new_item = item.split('[')
                                        #print(new_item)
                                        for it in new_item:
                                            #print(it)
                                            #print("chain: " + chain)
                                            if it == "auth " + str(chain) + "]":
                                                k=i
                                                #print("here")
                                                isfound=True
                                            elif it==chain:
                                                k=i
                                    elif item== " " + str(chain):
                                        k=i
                                elif item != chain_lst[0]:
                                    idx = item.find("auth")
                                    #print(idx)
                                    if idx!=-1:
                                        new_item = item.split('[')
                                        #print(new_item)
                                        for it in new_item:
                                            #print(it)
                                            #print("chain: " + chain)
                                            #print("auth " + str(chain) + "]")
                                            if it== "auth " + str(chain) + "]":
                                                #print("here")
                                                k=i
                                                isfound=True
                                            elif it==chain:
                                                k=i
                                    elif item == " " + str(chain):
                                        k=i
                    i+=2
    
        
            seq_dict[id]= fas_seq[k+1]
        else:
            seq_dict[id]= fas_seq[1]
    #print(seq_dict)

        seq_names = list(seq_dict.keys())

        final_list = list(seq_dict.keys())
        groups = {}
        #print(final_list)
        for i in range(len(seq_names)):
            for k in range(i+1,len(seq_names)):
                #print(i)
                #print(k)
                seq1_name = seq_names[i]
                seq2_name = seq_names[k]
                seq1 = seq_dict[seq1_name]
                seq2 = seq_dict[seq2_name]
                rem = seq2_name
                if len(seq1) >= len(seq2):
                    min = len(seq2)
                else:
                    min = len(seq1)
        
        
                score = pairwise2.align.globalxx(seq1, seq2,score_only=True)
                #print("score:",score)
                #print ("min: ", ((min*90)/100))
                if score >= ((min*90)/100) and rem in final_list:
                    final_list.remove(rem)
                    if seq1_name not in groups.keys():
                        groups[seq1_name]=seq2_name
                    elif seq1_name in groups.keys():
                        groups[seq1_name]=groups[seq1_name] + " " + seq2_name
        file = filename.removesuffix(".xlsx")+".txt"
        w = open(dest + "/" + file, "w",  encoding='utf-8')
        i=1
        for key in groups.keys():
            w.write("Group" + " " + str(i)+":")
            w.write("\n")
            w.write(key)
            w.write("\n")
            seq_names.remove(key)
            elements = groups[key].split()
            for element in elements:
                w.write(element)
                w.write("\n")
                seq_names.remove(element)
            i = i+1
        if len(seq_names)>0:
            for seq in seq_names:
                w.write("Group" + " " + str(i)+":")
                w.write("\n")
                w.write(seq)
                w.write("\n")
                i=i+1
        w.close()
