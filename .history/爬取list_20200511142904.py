"""
尝试通过biopython 下载accession list
"""

# -*- coding: UTF-8 -*-
import re
from io import StringIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import datetime
from Bio import Entrez

term = 'txid11118 [Organism]'
filename = 'C:\\Users\\ly\\Desktop\\txid11118_acc_list'
Entrez.email = "sunxiuqiang15@mails.ucas.ac.cn"
search_handle = Entrez.esearch(db="nucleotide",term=term,usehistory="y")
search_results = Entrez.read(search_handle)
search_handle.close()
count = int(search_results["Count"])
webenv = search_results["WebEnv"]
query_key = search_results["QueryKey"]
batch_size = 500
out_handle = open(filename, "w")
for start in range(0,count,batch_size):
    end = min(count, start+batch_size)
#    print ("Going to download record %i to %i" % (start+1, end))report=accnlist
    fetch_handle = Entrez.efetch(db="nucleotide", report="accnlist", retmode="text", retstart=start, retmax=batch_size, webenv=webenv, query_key=query_key)
    #fetch_handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", retstart=start, retmax=batch_size, webenv=webenv, query_key=query_key)
    data = fetch_handle.read()
    fetch_handle.close()
    out_handle.write(data)
out_handle.close()

print ("over")