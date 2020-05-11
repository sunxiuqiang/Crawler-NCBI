#!/usr/bin/python
#coding=utf-8

#LOCUS DEFINITION  submitDateFormat 都会有，其他字段需要判断
import re
import os
import argparse
import sys
from io import StringIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Entrez
import datetime
import socket
import random
import time,requests

timeout = 20
socket.setdefaulttimeout(timeout)
sleep_download_time = random.randint(0,6)

date = datetime.datetime.now().strftime('%Y_%m_%d')

Entrez.email = "sunxiuqiang15@mails.ucas.ac.cn"
'''
用于 nucleotide 数据库的下载，下载新增的fasta  和 gbk文件，以及对应genome 且长度在15000以上的locus的gff文件。
'''
#============================== 参数设置 ==============================
parser = argparse.ArgumentParser(prog='pipline',usage='%(prog)s [opthions] [value]',description='根据输入文件中需要更新的version，生成相应的文件')

parser.add_argument('--in',help='需要更新的version文件')
parser.add_argument('--dir',help='输出文件的路径')
parser.add_argument('--pro_num',help='Protein表的条数')
parser.add_argument('--paper_IN_num',help='ID2paper表的条数')
parser.add_argument('--Gene_num',help='Gene表的条数')
parser.add_argument('--fasta_cds_aa',help='fasta_cds_aa的文件路径')
parser.add_argument('--nucleotide_file',help='nucleotide_file的文件路径')
parser.add_argument('--pro_GI_ID_file',help='pro_GI_ID_file的文件路径')
parser.add_argument('--uniprot_pdb_file',help='uniprot_pdb_file的文件路径')
parser.add_argument('--fasta_cds_na',help='fasta_cds_na的文件路径')

argv = vars(parser.parse_args())
#==============================判断参数并设置默认参数==============================
if argv['in'] == None:
    sys.exit('请填写输入文件参数\n')
else:
    in_dir = argv['in'].strip()

if argv['dir'] == None:
    sys.exit('请填写输出文件路径\n')
else:
    dir = argv['dir'].strip()
if argv['pro_num'] == None:
    Protein_id_num = 1
else:
    Protein_id_num = argv['pro_num'].strip()
if argv['paper_IN_num'] == None:
    sys.exit('paper_IN_num ERROR!\n')
else:
    paper_id_num = argv['paper_IN_num'].strip()
if argv['Gene_num'] == None:
    Gene_id_num = 1
else:
    Gene_id_num = argv['Gene_num'].strip()

if argv['fasta_cds_aa'] == None:
    sys.exit('请填写fasta_cds_aa文件路径\n')
else:
    fasta_cds_aa = argv['fasta_cds_aa'].strip()
if argv['nucleotide_file'] == None:
    sys.exit('请填写nucleotide_file文件路径\n')
else:
    nucleotide_file = argv['nucleotide_file'].strip()
if argv['pro_GI_ID_file'] == None:
    sys.exit('请填写pro_GI_ID_file文件路径\n')
else:
    pro_GI_ID_file = argv['pro_GI_ID_file'].strip()
if argv['uniprot_pdb_file'] == None:
    sys.exit('请填写uniprot_pdb_file文件路径\n')
else:
    uniprot_pdb_file = argv['uniprot_pdb_file'].strip()

if argv['fasta_cds_na'] == None:
    sys.exit('请填写fasta_cds_na文件路径\n')
else:
    fasta_cds_na = argv['fasta_cds_na'].strip()



input_file = "%s/update_GI_%s"%(in_dir,date)
file1 = open("%s"%(input_file))


if os.path.isfile("%s/Coronaviridae_protein_%s/temp/add_protein_info_%s.txt"%(dir,date,date)):
    os.remove("%s/Coronaviridae_protein_%s/temp/add_protein_info_%s.txt"%(dir,date,date))
if os.path.isfile("%s/temp/add_Gene_info_%s.txt"%(dir,date)):
    os.remove("%s/temp/add_Gene_info_%s.txt"%(dir,date))
if os.path.isfile("/home/sxq/NCBI/data/download_data/mysql/nucleotide/paper/add_paper_P_info_%s.txt"%(date)):
    os.remove("/home/sxq/NCBI/data/download_data/mysql/nucleotide/paper/add_paper_P_info_%s.txt"%(date))



file_mysql_pro = open ("%s/temp/add_protein_info_%s.txt"%(dir,date),"a")
file_mysql_pro.write("Protein_ID\tNucleotide_ID\tNucleotideLocus\tSpciesname\tStrain\tGeneID\tProteinLocus\tGI\tProteinDefinition\tProteinID\tProteinAccession\tProteinClassification\tProteinLength\tProteinSequence\tUniProtID\tPDB_ID\tSource\n")
file_mysql_Gene = open("%s/temp/add_Gene_info_%s.txt"%(dir,date),"a")
file_mysql_Gene.write("Gene_ID\tNucleotide_ID\tProtein_ID\tGeneID\tgene\tProteinClassification\tGeneSymbol\tGeneDescription\tLocusTag\tOrient\tStart\tEnd\tGeneLength\tNucleotideSequence\n")
file_mysql_pmid = open ("/home/sxq/NCBI/data/download_data/mysql/nucleotide/paper/temp/add_paper_P_info_%s.txt"%(date),"a")
file_mysql_pmid.write("ID2Paper_id\tEntry_ID\tEntryType\tPaper_ID\tPaperType\tSource\tTitle\tJOURNAL\tAUTHORS\n")


def get_Gene_info(fasta_cds_aa):
    """
    从fasta_cds_aa 文件中提取相关的字符
    """
    list_protein_id = []
    dic_pro_id_gene = {}
##=============================对fasta_cds_aa文件中的信息进行提取，包括 核酸version protein id gene 在基因组的位置信息等===========================
    record_iterator  = SeqIO.parse(fasta_cds_aa, "fasta")
    for seq_record in record_iterator:
        id_line = seq_record.id
        #Nucle_version = id_line.split("_prot_")[0].lstrip("lcl|")
        protein_ID = "_".join(id_line.split("_prot_")[1].split("_")[:-1])
        result = seq_record.description
        #print ("descrip")
        #result = "#".join(re.findall(r'\[(.*?)\]', descrip))
        #print ("result",result)
        if "gene=" in result:
            gene = re.match(r'.*gene=(.*?)].*', result).group(1)
            #print ("gene",gene)
        else:
            gene = " "
        if "protein=" in result:
            protein = re.match(r'.*protein=(.*?)].*', result).group(1)
        else:
            protein = " "       
        if "location=" in result:
            loca = re.findall(r'.*location=(.*?)].*', result)[0]
            #print ("loca:",loca)
            if "join" in loca:
                orient = "+"
                location = loca.lstrip("join(").rstrip(")").split("..")[0]+".."+loca.lstrip("join(").rstrip(")").split("..")[-1]
                if "<" in location and ">" not in location:
                    start = location.split("..")[0].lstrip("<")
                    end = location.split("..")[1]
                elif "<" not in location and ">" in location:
                    start =  location.split("..")[0]
                    end = location.split("..")[1].lstrip(">")
                elif "<" in location and ">" in location:
                    start =  location.split("..")[0].lstrip("<")
                    end = location.split("..")[1].lstrip(">")
                else:
                    start = location.split("..")[0]
                    end = location.split("..")[1]
            elif "complement(" in loca:
                orient = "-"
                location = loca.lstrip("complement(").rstrip(")")
                if "<" in location and ">" not in location:
                    start = location.split("..")[0].lstrip("<")
                    end = location.split("..")[1]
                elif "<" not in location and ">" in location:
                    start = location.split("..")[0]
                    end = location.split("..")[1].lstrip(">")
                elif "<" in location and ">" in location:
                    start = location.split("..")[0].lstrip("<")
                    end = location.split("..")[1].lstrip(">")
                else:
                    start = location.split("..")[0]
                    end = location.split("..")[1]
            else:
                orient = "+"
                if "<" in loca and ">" not in loca:
                    start = loca.split("..")[0].lstrip("<")
                    end = loca.split("..")[1]
                elif "<" not in loca and ">" in loca:
                    start = loca.split("..")[0]
                    end = loca.split("..")[1].lstrip(">")
                elif "<" in loca and ">" in loca:
                    start = loca.split("..")[0].lstrip("<")
                    end = loca.split("..")[1].lstrip(">")
                else:
                    start = loca.split("..")[0]
                    end = loca.split("..")[1]
        ProteinLength =str(len(seq_record.seq))
        ProteinSequence = str(seq_record.seq)
        dic_pro_id_gene[protein_ID] = [gene,protein,orient,start,end,ProteinLength,ProteinSequence]
    return dic_pro_id_gene

def nucleotide_info(nucleotide_file):
    """
    从/mysql/nucleotide/nucleotide_info 中获得核酸表所需要的相关字段，包括Nucleotide_ID NucleotideLocus Spciesname Strain
    """
    dic_Nucle_version_info = {} 
    #dic_Nucle_version_locus = {}
    with open(nucleotide_file) as nucleotide_file:
        title = nucleotide_file.readline()
        for eachline in nucleotide_file:
            line = eachline.strip().split("\t")
            Nucleotide_locus = line[1]
            Nucleotide_ID = line[0]
            Nucleotide_Version = line[3]
            Spciesname = line[5]
            Strain = line[10]
            dic_Nucle_version_info[Nucleotide_Version] = [Nucleotide_ID,Nucleotide_locus,Nucleotide_Version,Spciesname,Strain]
            #dic_Nucle_version_locus[Nucleotide_Version] = Nucleotide_locus
    return dic_Nucle_version_info

def Uniprot_info(pro_GI_ID_file,uniprot_pdb_file):
    """
    根据蛋白的GI去uniprot 中获得 uniprot id 和PDB id的信息  蛋白GI信息home/sxq/NCBI/scripts/test/all_ncbi/get_GI/protein_GI/locus_pro_GI_2020_03_17" 这个需要考虑到数据更新的问题  "/home/sxq//NCBI/scripts/test/all_ncbi/get_GI/idmapping_selected.tab" 这个文件 来自uniprot 数据库  8周更新一次
    """
    dic_GI_Uniprot = {}
    dic_GI_EntrezGene = {}
    dic_GI_PDB = {}
    dic_protein_id_GI = {}
    dic_protein_id_Nucle_version = {}

    with open(pro_GI_ID_file) as GI_file:
        title = GI_file.readline()
        for eachline in GI_file:
            line = eachline.strip().split(",")
            protein_id = line[1]      
            dic_protein_id_GI[protein_id] = line[2]
            dic_protein_id_Nucle_version[protein_id] = line[0]   #通过protein_id 得到对应的核酸version号
            
    with open(uniprot_pdb_file) as uniprot_pdb_file:
        for eachline in uniprot_pdb_file:
            line = eachline.strip().split("\t")
            GI_list = line[4]
            gi = GI_list.split(";")
            for each in gi:
                GI = each.strip()
                dic_GI_Uniprot[GI] = line[0]
                dic_GI_EntrezGene[GI] = line[2]  # 即Gene库里的Gene ID
                dic_GI_PDB[GI] = line[5]
    return [dic_GI_PDB,dic_GI_Uniprot,dic_GI_EntrezGene,dic_protein_id_GI,dic_protein_id_Nucle_version] #得到protein_id 和核酸locus，protein GI的对应关系,以及GI 和uniprot pdb GeneID的对应关系 
def protein_id_list(input_file):
    """
    获得蛋白库中对应物种的protein id list，因为有些protein id 在核酸库中不存在
    """
    pro_id_list = []
    with open(input_file) as protein_id_list_file:
        for eachline in protein_id_list_file:
            line = eachline.strip().split(",")
            protein_id = line[1]
            pro_id_list.append(protein_id)
    return pro_id_list

def pro_gp_file(filename,paper_IN_num,Protein_ID_num):
    """
    从gp文件中提取ProteinDefinition ProteinAccession ProteinClassification
    """
    dic_pro_id_accssion_definition = {}
    record = SeqIO.read(filename, "genbank")
    ProteinAccession = record.name
    ProteinDefinition = record.description
    pro_id = record.id
    dic_pro_id_accssion_definition[pro_id] = [ProteinAccession,ProteinDefinition]
    #paper_IN_num =int(paper_IN_num)
    for each in record.annotations["references"]:
        for each_id in each.pubmed_id:
            if each_id != "":
                paper_IN_num +=1
                EntryType = "P"
                Paper_ID = each.pubmed_id
                PaperType = "PMID"
                Source = "PubMed"
                Title = each.title
                JOURNAL = each.journal
                AUTHORS = each.authors
                file_mysql_pmid.write(str(paper_IN_num)+"\t"+str(Protein_ID_num)+"\t"+EntryType+"\t"+Paper_ID+"\t"+PaperType+"\t"+Source+"\t"+Title+"\t"+JOURNAL+"\t"+AUTHORS+"\n")
    #print ("dic_pro_id_accssion_definition",dic_pro_id_accssion_definition)
    return dic_pro_id_accssion_definition
def fasta_cds_na_info(fasta_cds_na):   
    """
    核酸库下载的fasta_cds_na文件，提取protein id 对应的核酸序列及序列长度
    """
    dic_pro_id_nucle_length = {}
    dic_pro_id_nucle_sequence = {}
    record_iterator  = SeqIO.parse(fasta_cds_na, "fasta")
    for seq_record in record_iterator:
        protein_id = "_".join(seq_record.name.split("_cds_")[1].split("_")[:-1])
        dic_pro_id_nucle_length[protein_id] = str(len(seq_record.seq))
        dic_pro_id_nucle_sequence[protein_id] = str(seq_record.seq)
    return [dic_pro_id_nucle_length,dic_pro_id_nucle_sequence]


def main():
    pro_id_not_nucle = open("/home/sxq/NCBI/data/download_data/mysql/protein/pro_id_not_nucle_%s.txt"%(date),"a") #记录那些尽在蛋白库中才有的protein id，核酸库里没有此protein id
    pro_id_GI = open("/home/sxq/NCBI/data/update_ncbi/protein/Nucle_locus_proID_proGI/protein_id_pro_GI","a")  #记录Nucle_locus_proID_proGI_2020_04_26中没有的的GI ，两列： protein ID 和GI
    pro_id_not_nucle.write("protein_id\tGI\n")
    count_number = 0
    Protein_ID_num = int(Protein_id_num)
    Gene_ID_num = int(Gene_id_num)
    paper_IN_num = int(paper_id_num)
    f = get_Gene_info(fasta_cds_aa)
    #nucleotide_file = "/home/sxq/NCBI/data/download_data/mysql/nucleotide/file/nucleotide_info.txt"  
    f2 = nucleotide_info(nucleotide_file)
    #print ("nucleotide_info:",f2)
    #pro_GI_ID_file = "/home/sxq/NCBI/data/update_ncbi/protein/Nucle_locus_proID_proGI/Nucle_locus_proID_proGI_2020_04_26"  ##需要把此文件 换成每日更新的file  3列：核酸locus protein_id  protein_id的GI号
    #uniprot_pdb_file = "/home/sxq/NCBI/data/update_ncbi/protein/UniProt/idmapping_selected.tab"    ##此文件8周更新一次
    f3 = Uniprot_info(pro_GI_ID_file,uniprot_pdb_file)
    #print ("Uniprot_info",f3)
    #protein_id_list_file = "/home/sxq/NCBI/data/update_ncbi/protein/update_GI_2020_04_20_test"    ##需要把此文件 换成每日更新的file  两列信息 蛋白locus 蛋白的version protein_id
    
    #print ("pro_gp_file",f4)
    #fasta_cds_na = "/home/sxq/NCBI/data/update_ncbi/nucleotide/fasta_cds_na/txid11118_fasta_cds_na_2020_04_26.txt"  ##需要把此文件 换成每日更新的file
    f5 = fasta_cds_na_info(fasta_cds_na)
    #print ("fasta_cds_na_info",f5)
    #print ("GeneID:",f3[2].get(f3[3].get(pro_id)))
    for eachline in file1:  #输入file1  每日更新的蛋白库数据，两列信息，pro_locus和pro_id
        count_number +=1
        line = eachline.strip().split(",")
        pro_locus = line[0]
        pro_id = line[1]
        print ("正在下载第%s个数据，pro_id:%s"%(count_number,pro_id))
        filename = "%s/Coronaviridae_protein/%s.gp"%(dir,pro_locus)
        
        if os.path.isfile(filename):
            if pro_id in f3[3]:
                GI = f3[3].get(pro_id)
            else:
                time.sleep(sleep_download_time)
                gi_handle = Entrez.esummary(db="protein", id=pro_id)
                record = Entrez.read(gi_handle)
                GI = record[0]["Id"]
                pro_id_GI.write(pro_id+","+GI+"\n")
            #print ("GI:",GI)
            #fasta_cds_aa = fasta_cds_aa ##需要把此文件 换成每日更新的file
            if pro_id in f:
                Protein_ID_num +=1
                Gene_ID_num +=1 
                f4 = pro_gp_file(filename,paper_IN_num,Protein_ID_num)
                Protein_ID = str(Protein_ID_num)
                Nucleotide_ID = f2.get(f3[4].get(pro_id))[0]
                NucleotideLocus = f2.get(f3[4].get(pro_id))[1]
                #NucleotideVersion = f2.get(f3[4].get(pro_id))[2]
                Spciesname = f2.get(f3[4].get(pro_id))[3]
                Strain = f2.get(f3[4].get(pro_id))[4]
                GI = f3[3].get(pro_id)
                ProteinClassification = " "
                ProteinLength = f.get(pro_id)[5]
                ProteinSequence = f.get(pro_id)[6]
                Source = "NCBI"
                if f3[3].get(pro_id) in f3[1]: #dic_GI_Uniprot
                    UniProtID = f3[1].get(f3[3].get(pro_id))
                    if f3[3].get(pro_id) in f3[0]:  #dic_GI_PDB
                        PDB_ID = f3[0].get(f3[3].get(pro_id))  
                        if f3[3].get(pro_id) in f3[2]: #dic_GI_EntrezGene
                            GeneID = f3[2].get(f3[3].get(pro_id))
                        else:
                            GeneID = " "
                    else:
                        if f3[3].get(pro_id) in f3[2]: #dic_GI_EntrezGene
                            GeneID = f3[2].get(f3[3].get(pro_id))
                            PDB_ID = " "    
                        else:
                            PDB_ID = " "   
                            GeneID = " "            
                else:
                    if f3[3].get(pro_id) in f3[0]:  #dic_GI_PDB
                        PDB_ID = f3[0].get(f3[3].get(pro_id))  
                        UniProtID = " "
                        if f3[3].get(pro_id) in f3[2]: #dic_GI_EntrezGene
                            GeneID = f3[2].get(f3[3].get(pro_id))
                        else:
                            GeneID = " "
                    else:
                        PDB_ID = " "
                        GeneID = " "
                        UniProtID = " "

                #print (Protein_ID+"\t"+Nucleotide_ID+"\t"+NucleotideLocus+"\t"+Spciesname+"\t"+GeneID+"\t"+pro_locus+"\t"+GI+"\t"+f4.get(pro_id)[1]+"\t"+pro_id+"\t"+f4.get(pro_id)[0]+"\t"+ProteinClassification+"\t"+ProteinLength+"\t"+UniProtID+"\t"+PDB_ID+"\t"+Source+"\n") 
                file_mysql_pro.write(Protein_ID+"\t"+Nucleotide_ID+"\t"+NucleotideLocus+"\t"+Spciesname+"\t"+Strain+"\t"+GeneID+"\t"+pro_locus+"\t"+GI+"\t"+f4.get(pro_id)[1]+"\t"+pro_id+"\t"+f4.get(pro_id)[0]+"\t"+ProteinClassification+"\t"+ProteinLength+"\t"+ProteinSequence+"\t"+UniProtID+"\t"+PDB_ID+"\t"+Source+"\n")            
                file_mysql_Gene.write(str(Gene_ID_num)+"\t"+Nucleotide_ID+"\t"+Protein_ID+"\t"+GeneID+"\t"+f.get(pro_id)[0]+"\t \t \t \t \t"+f.get(pro_id)[2]+"\t"+f.get(pro_id)[3]+"\t"+f.get(pro_id)[4]+"\t"+f5[0].get(pro_id)+"\t"+f5[1].get(pro_id)+"\n")
            else:
                print ("这些protein id 不在核酸库中：",pro_id)
                pro_id_not_nucle.write(pro_id+"\t"+GI+"\n")            
        else:
            time.sleep(sleep_download_time)
            gi_handle = Entrez.esummary(db="protein", id=pro_id)
            record = Entrez.read(gi_handle)
            GI = record[0]["Id"]
            pro_id_GI.write(pro_id+","+GI+"\n")
            print ("在下载gp文件呀。。。\n")
            net_handle = Entrez.efetch(db="protein",id=GI,rettype="gp", retmode="text")
            out_handle = open(filename, "w")
            out_handle.write(net_handle.read())
            out_handle.close()
            net_handle.close()
            if pro_id in f:
                Protein_ID_num +=1
                Gene_ID_num +=1 
                f4 = pro_gp_file(filename,paper_IN_num,Protein_ID_num)
                Protein_ID = str(Protein_ID_num)
                Nucleotide_ID = f2.get(f3[4].get(pro_id))[0]
                NucleotideLocus = f2.get(f3[4].get(pro_id))[1]
                Spciesname = f2.get(f3[4].get(pro_id))[3]
                Strain = f2.get(f3[4].get(pro_id))[4]
                GI = f3[3].get(pro_id)
                ProteinClassification = " "
                ProteinLength = f.get(pro_id)[5]
                ProteinSequence = f.get(pro_id)[6]
                Source = "NCBI"
                if f3[3].get(pro_id) in f3[1]: #dic_GI_Uniprot
                    UniProtID = f3[1].get(f3[3].get(pro_id))
                    if f3[3].get(pro_id) in f3[0]:  #dic_GI_PDB
                        PDB_ID = f3[0].get(f3[3].get(pro_id))  
                        if f3[3].get(pro_id) in f3[2]: #dic_GI_EntrezGene
                            GeneID = f3[2].get(f3[3].get(pro_id))
                        else:
                            GeneID = " "
                    else:
                        if f3[3].get(pro_id) in f3[2]: #dic_GI_EntrezGene
                            GeneID = f3[2].get(f3[3].get(pro_id))
                            PDB_ID = " "    
                        else:
                            PDB_ID = " "   
                            GeneID = " "            
                else:
                    if f3[3].get(pro_id) in f3[0]:  #dic_GI_PDB
                        PDB_ID = f3[0].get(f3[3].get(pro_id))  
                        UniProtID = " "
                        if f3[3].get(pro_id) in f3[2]: #dic_GI_EntrezGene
                            GeneID = f3[2].get(f3[3].get(pro_id))
                        else:
                            GeneID = " "
                    else:
                        PDB_ID = " "
                        GeneID = " "
                        UniProtID = " "

                #print (Protein_ID+"\t"+Nucleotide_ID+"\t"+NucleotideLocus+"\t"+Spciesname+"\t"+GeneID+"\t"+pro_locus+"\t"+GI+"\t"+f4.get(pro_id)[1]+"\t"+pro_id+"\t"+f4.get(pro_id)[0]+"\t"+ProteinClassification+"\t"+ProteinLength+"\t"+UniProtID+"\t"+PDB_ID+"\t"+Source+"\n") 
                file_mysql_pro.write(Protein_ID+"\t"+Nucleotide_ID+"\t"+NucleotideLocus+"\t"+Spciesname+"\t"+Strain+"\t"+GeneID+"\t"+pro_locus+"\t"+GI+"\t"+f4.get(pro_id)[1]+"\t"+pro_id+"\t"+f4.get(pro_id)[0]+"\t"+ProteinClassification+"\t"+ProteinLength+"\t"+ProteinSequence+"\t"+UniProtID+"\t"+PDB_ID+"\t"+Source+"\n")            
                file_mysql_Gene.write(str(Gene_ID_num)+"\t"+Nucleotide_ID+"\t"+Protein_ID+"\t"+GeneID+"\t"+f.get(pro_id)[0]+"\t \t \t \t \t"+f.get(pro_id)[2]+"\t"+f.get(pro_id)[3]+"\t"+f.get(pro_id)[4]+"\t"+f5[0].get(pro_id)+"\t"+f5[1].get(pro_id)+"\n")
            else:
                print ("这些protein id 不在核酸库中：",pro_id)
                pro_id_not_nucle.write(pro_id+"\t"+GI+"\n")                
    print("数据更新现在结束,一共更新%s条数据=============================================================================\n"%(count_number))

    
if __name__ == '__main__': 
    main()
    file_mysql_pro.close()
    file_mysql_Gene.close()
    file_mysql_pmid.close()


os.system("tail -1  /home/sxq/NCBI/data/download_data/mysql/protein/temp/add_protein_info_%s.txt |awk '{print $1}'>> /home/sxq/NCBI/data/download_data/mysql/protein/config/number_protein"%(date))
os.system("tail -1 /home/sxq/NCBI/data/download_data/mysql/protein/temp/add_Gene_info_%s.txt|awk '{print $1}' >>/home/sxq/NCBI/data/download_data/mysql/protein/config/number_Gene"%(date))
#os.system("tail -1 /home/sxq/NCBI/data/download_data/mysql/nucleotide/paper/add_paper_P_info_%s.txt|awk '{print $1}' >>/home/sxq/NCBI/data/download_data/mysql/nucleotide/config/number_paper"%(date))

