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
sleep_download_time = random.randint(1,5)

date = datetime.datetime.now().strftime('%Y_%m_%d')

Entrez.email = "sunxiuqiang15@mails.ucas.ac.cn"
'''
用于 nucleotide 数据库的下载，下载新增的fasta  和 gbk文件，以及对应genome 且长度在15000以上的locus的gff文件。
'''
#============================== 参数设置 ==============================
parser = argparse.ArgumentParser(prog='pipline',usage='%(prog)s [opthions] [value]',description='根据输入文件中需要更新的version，生成相应的文件')

parser.add_argument('--in',help='需要更新的version文件')
parser.add_argument('--dir',help='输出文件的路径')
parser.add_argument('--num',help='核酸的数目')
parser.add_argument('--paper_IN_num',help='')
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

if argv['num'] == None:
    num = 1
else:
    num = argv['num'].strip()
if argv['paper_IN_num'] == None:
    paper_IN_num = 1
else:
    paper_IN_num = argv['paper_IN_num'].strip()

input_file = "%s/update_GI_%s"%(in_dir,date)
file1 = open("%s"%(input_file))

if not os.path.exists("%s/Coronaviridae_genbank_%s/gffs"%(dir,date)):
    os.makedirs("%s/Coronaviridae_genbank_%s/gffs"%(dir,date))
if not os.path.exists("%s/Coronaviridae_genbank_%s/fasta"%(dir,date)):
    os.makedirs("%s/Coronaviridae_genbank_%s/fasta/"%(dir,date))
if os.path.isfile("%s/Coronaviridae_genbank_%s/add_info_%s.txt"%(dir,date,date)):
    os.remove("%s/Coronaviridae_genbank_%s/add_info_%s.txt"%(dir,date,date))
if not os.path.exists("%s/gyc_%s/"%(dir,date)):
    os.makedirs("%s/gyc_%s/"%(dir,date))


file_out = open ("%s/Coronaviridae_genbank_%s/add_info_%s.txt"%(dir,date,date),"a")
file_out.write("locus\tspciesname\tdefinition\tisolate\tstrain\thost\tisolation_source\tglength\tcollection_date\tyear\tcountry\tsubmitdate\tcollectionDateFormat\tsubmitDateFormat\tGenome\tSARS-CoV-2\n")
file_mysql = open ("%s/mysql/nucleotide/add_nucleotide_info_%s.txt"%(dir,date),"a")
file_mysql.write("NucleotideID\tNucleotideLocus\tNucleotidAccession\tNucleotidVersion\tisRefSeq\tSpciesName\tTaxonID\tDefinition\tType\tIsolate\tStrain\tHost\tHostFormat\tIsolationSource\tGlength\tCollectionDate\tYear\tCountry\tCountryFormat\tDistrict\tSubmitDate\tSubmitLocation\tCollectionDateFormat\tSubmitDateFormat\tGenome\tSARS-CoV-2\tXref_ID\tSource\n")
file_mysql_pmid = open ("%s/mysql/nucleotide/paper/temp/add_paper_Nucle_info_%s.txt"%(dir,date),"a")
file_mysql_pmid.write("ID2Paper_id\tEntry_ID\tEntryType\tPaper_ID\tPaperType\tSource\tTitle\tJOURNAL\tAUTHORS\n")

def  extract_info(filename,version,locus,Nucleotide_ID_num,paper_IN_num,GI):
    genome_gff = []
    paper_IN_num = int(paper_IN_num)
    record = SeqIO.read(filename, "genbank")
    for each in record.annotations["references"]:
        if each.title == "Direct Submission":
            SubmitLocation = each.journal.split(",")[-1].strip()
        else:
            SubmitLocation = " "
        for each_id in each.pubmed_id:
            if each_id != "":
                paper_IN_num +=1
                Entry_ID = locus
                EntryType = "N"
                Paper_ID = each.pubmed_id
                PaperType = "PMID"
                Source = "PubMed"
                Title = each.title
                JOURNAL = each.journal
                AUTHORS = each.authors
                file_mysql_pmid.write(str(paper_IN_num)+"\t"+str(Nucleotide_ID_num)+"\t"+EntryType+"\t"+Paper_ID+"\t"+PaperType+"\t"+Source+"\t"+Title+"\t"+JOURNAL+"\t"+AUTHORS+"\n")
    with open(filename) as input_file:
        info = "".join(input_file.readlines()).replace("\n", "")
        #print (info)
        locus_length = re.findall("LOCUS\s+\w+\s+\w+\sbp",info)[0].split()
        length = locus_length[2]
        accession = re.findall(r'.*ACCESSION\s+(.*)?VERSION',info)[0].replace("  ","")
        keywords = re.findall(r'.*KEYWORDS\s+(.*)?SOURCE',info)[0]
        if "RefSeq." in keywords:
            isRefseq = "1"
        else:
            isRefseq = "0"
        definition = re.findall(r'.*DEFINITION\s+(.*)?ACCESSION',info)[0].rstrip('.').replace("  ","")
        definite = definition.split(",")[-1].strip()
        if  definite == "complete genome" or definite == "completegenome":
            Type = "complete genome"
        elif definite == "partial genome" or definite == "partialgenome":
            Type = "partial genome"
        elif definite == "partial cds" or definite == "partialcds":
            Type = "partial cds"
        elif definite == "complete cds" or definite == "completecds":
            Type = "complete cds"
        else:
            Type = " "

        if "genome" in definition and int(length) >15000 or "genomic" in definition and int(length) >15000:
            genome = "Yes"
            genome_b = "1"
            genome_gff.append(locus)
            gff_file = "%s/Coronaviridae_genbank_%s/gffs/%s.gff3"%(dir,date,locus)
            if not os.path.isfile(gff_file):
##======================写出gff文件===================================================================
                gff_url = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=gff3&id=%s"%(GI)
                time.sleep(sleep_download_time)
                user_agents=['Mozilla/5.0 (Windows NT 6.1; rv:2.0.1) Gecko/20100101 Firefox/4.0.1','Mozilla/5.0 (Windows; U; Windows NT 6.1; en-us) AppleWebKit/534.50 (KHTML, like Gecko) Version/5.1 Safari/534.50','Opera/9.80 (Windows NT 6.1; U; en) Presto/2.8.131 Version/11.11']
                proxy_list = [ '74.59.132.126:49073','176.117.255.182:53100',
                    '39.137.69.7:80','112.195.81.243:8118']
                proxy = random.choice(proxy_list)
                pxs = {'http': 'http://' + proxy}
                headers={'User-Agent':random.choice(user_agents)}

                r = requests.request('GET',gff_url,proxies = pxs,headers=headers)
                with open ("%s//Coronaviridae_genbank_%s/gffs/%s.gff3"%(dir,date,locus),"wb") as gff_file:
                    gff_file.write(r.content)
                with open ("%s//Coronaviridae_genbank_%s/gffs/%s.gff3"%(dir,date,locus)) as f1:
                    alllines = f1.readlines()
                with open("%s//Coronaviridae_genbank_%s/gffs/%s.gff3"%(dir,date,locus),"w+") as final_gff:
                    for eachline in alllines:
                        line =re.sub(version,locus,eachline)
                        final_gff.writelines(line)
        else:
            genome = " "
            genome_b = "0"
        if "/isolate" in info:
            isolate = re.findall(r'.*/isolate="(.*)?"',info)[0].split('"')[0]
        #isolate = re.findall(r'.*/isolate="(.*)"\s+',info)
        else:
            isolate = " "
        if "/strain=" in info:
            strain = re.findall(r'.*/strain="(.*)?"',info)[0].split('"')[0]
        else:
            strain = " "
        if "/host=" in info:
            host = re.findall(r'.*/host="(.*)?"',info)[0].split('"')[0].split(";")[0]
            if host == "Homo sapiens":
                HostFormat = "Human"
            else:  
                HostFormat = " "
 
        else:
            host = " "
            HostFormat = " "
        if "/collection_date=" in info:
            collection_date = re.findall(r'.*/collection_date="(.*)?"',info)[0].split('"')[0]
            print (collection_date)
            collectionDateFormat = collection_date
        else:
            collection_date = " "
            collection_Date_Format = "0001-01-01 "
        if "/country=" in info:
            country = re.findall(r'.*/country="(.*)?"',info)[0].split('"')[0]
            if country == "Taiwan":
                country = "China: Taiwan"
            elif country == "Hong Kong":
                country = "China: Hong Kong"
            con = country.split(":")
            if len(con) ==2:
                CountryFormat = con[0]
                LocationFormat = con[1].strip()
            else:
                CountryFormat = con[0]
                LocationFormat = " "
        else:
            country = " "
            CountryFormat = " "
            LocationFormat = " "
        if "/isolation_source=" in info:
            isolationSource = re.findall(r'.*/isolation_source="(.*)?"',info)[0].split('"')[0].replace("  ","")
        else:
            isolationSource = " "
        spciesname = re.findall(r'.*SOURCE\s+(.*)\s+ORGANISM',info)[0].split('"')[0]
        if "Severe acute respiratory syndrome coronavirus 2" in spciesname or "SARS-CoV-2" in definition:
            SARS_CoV_2 = "Yes"
            SARS_CoV_2_b ="1"
        else:
            SARS_CoV_2 = " "
            SARS_CoV_2_b = "0"
        if "JOURNAL   Submitted (" in info:
            p1 = re.compile(r'.*JOURNAL\s+Submitted\s+[(](.*?)[)]', re.S)
            submitDateFormat = re.findall(p1,info)[0]
            submit_date = re.findall(p1,info)[0]
        else: 
            submit_date = " "
            year = "0001"
            submitDateFormat = " "
            submit_Date_Format = "0001-01-01"
        fasta = ' '.join(re.findall(r'.*ORIGIN\s+(.*)//',info)[0].split())
        fasta_sequence = re.sub(r'\d*\s*','',fasta)
##==========================写出fasta序列信息===================================================
        with open ("%s/Coronaviridae_genbank_%s/%s.fasta"%(dir,date,locus),"w") as fa_file:
            fa_file.write(">%s\n"%(locus))
            for eachbase in range(0,len(fasta_sequence),70):
                fa_file.write(''.join(fasta_sequence[eachbase:eachbase+70])+"\n")

        #print (len(fasta_sequence))
        if "/db_xref=" in info:
            TaxonID = re.findall(r'.*/db_xref="taxon:(.*)?"',info)[0].split('"')[0]
        else:
            TaxonID = ' '
    ##=========================调整时间格式=========================================================
        #month = ''
        if submitDateFormat != " ":    
            if "Jan" in submitDateFormat or "JAN"in submitDateFormat:
                month = '01'
            elif "Feb" in  submitDateFormat or "FEB" in  submitDateFormat:
                month = '02'
            elif "Mar" in  submitDateFormat or "MAR"in  submitDateFormat:
                month = '03'
            elif "Apr" in  submitDateFormat or "APR" in  submitDateFormat:
                month = '04'
            elif "May" in  submitDateFormat or "MAY" in  submitDateFormat:
                month = '05'
            elif "Jun" in  submitDateFormat or "JUN" in  submitDateFormat:
                month = '06'
            elif "Jul" in  submitDateFormat or "JUL"in  submitDateFormat:
                month = '07'
            elif "Aug" in  submitDateFormat or "AUG"in  submitDateFormat:
                month = '08'
            elif "Sept" in  submitDateFormat or "SEP" in submitDateFormat:
                month = '09'
            elif "Oct" in  submitDateFormat or "OCT"in  submitDateFormat:
                month = '10'
            elif "Nov" in  submitDateFormat or "NOV" in  submitDateFormat:
                month = '11'
            elif "Dec" in  submitDateFormat or "DEC" in submitDateFormat:
                month = '12'
            else:
                month = '01'
            submitDate = submitDateFormat.split("-")
            #if int(submitDate[0]) < 10:
            #    submit_Date_Format  = submitDate[2]+"-"+month+"-0"+submitDate[0]
            submit_Date_Format = submitDate[2]+"-"+month+"-"+submitDate[0]
            year = submitDate[2]
            #print (submit_Date_Format)
        if collection_date != " ":
            if len(collection_date.split("/")) ==1:       
                collectionDate = collectionDateFormat.split("-")            

                if len (collectionDate) == 3:
                    M = collectionDate[1]
                    if re.match(r'\D+',collectionDate[1]) != "None":
                        if "Jan" in collectionDateFormat or "JAN"in collectionDateFormat:
                            M = "01"
                        elif "Feb"  in collectionDateFormat or "FEB"in collectionDateFormat:
                            M = "02"
                        elif "Mar" in  collectionDateFormat or "MAR"in  collectionDateFormat:
                            M = '03'
                        elif "Apr" in  collectionDateFormat or "APR" in  collectionDateFormat:
                            M = '04'
                        elif "May" in  collectionDateFormat or "MAY" in  collectionDateFormat:
                            M = '05'
                        elif "Jun" in  collectionDateFormat or "JUN" in  collectionDateFormat:
                            M = '06'
                        elif "Jul" in  collectionDateFormat or "JUL"in  collectionDateFormat:
                            M = '07'
                        elif "Aug" in  collectionDateFormat or "AUG"in  collectionDateFormat:
                            M = '08'
                        elif "Sep" in  collectionDateFormat or "SEP" in collectionDateFormat:
                            M = '09'
                        elif "Oct" in  collectionDateFormat or "OCT"in  collectionDateFormat:
                            M = '10'
                        elif "Nov" in  collectionDateFormat or "NOV" in  collectionDateFormat:
                            M = '11'
                        elif "Dec" in  collectionDateFormat or "DEC" in collectionDateFormat:
                            M = '12'
                    else:
                        if collectionDate[1].startswith("0"):
                            M = collectionDate[1]
                        else:
                            if len(collectionDate[1]) == 1:
                                M = "0"+ collectionDate[1]
                            else:
                                M = collectionDate[1]
                        
                    if len(collectionDate[2]) != 4:
                        Day = collectionDate[2]
                        Year = collectionDate[0]
                        if len (collectionDate[2]) != 1:
                                Day = collectionDate[2]
                        else:
                            Day = "0"+ collectionDate[2]
                    else:
                        Year = collectionDate[2]
                        Day = collectionDate[0]
                        if len (collectionDate[2]) != 1:
                            Day = collectionDate[0]
                        else:
                            Day = "0"+ collectionDate[0]
                                
                elif len (collectionDate) == 2:
                    Day = "01"
                    if len(collectionDate[1]) != 4:
                        M = collectionDate[1]
                        Year =  collectionDate[0]
                        if "Jan" in collectionDateFormat or "JAN"in collectionDateFormat:
                            M = "01"
                        elif "Feb"  in collectionDateFormat or "FEB"in collectionDateFormat:
                            M = "02"
                        elif "Mar" in  collectionDateFormat or "MAR"in  collectionDateFormat:
                            M = '03'
                        elif "Apr" in  collectionDateFormat or "APR" in  collectionDateFormat:
                            M = '04'
                        elif "May" in  collectionDateFormat or "MAY" in  collectionDateFormat:
                            M = '05'
                        elif "Jun" in  collectionDateFormat or "JUN" in  collectionDateFormat:
                            M = '06'
                        elif "Jul" in  collectionDateFormat or "JUL"in  collectionDateFormat:
                            M = '07'
                        elif "Aug" in  collectionDateFormat or "AUG"in  collectionDateFormat:
                            M = '08'
                        elif "Sep" in  collectionDateFormat or "SEP" in collectionDateFormat:
                            M = '09'
                        elif "Oct" in  collectionDateFormat or "OCT"in  collectionDateFormat:
                            M = '10'
                        elif "Nov" in  collectionDateFormat or "NOV" in  collectionDateFormat:
                            M = '11'
                        elif "Dec" in  collectionDateFormat or "DEC" in collectionDateFormat:
                            M = '12'
                        else:
                            if collectionDate[1].startswith("0"):
                                M = collectionDate[1]
                            else:
                                if len(collectionDate[1]) ==1:
                                    M = "0"+ collectionDate[1]
                                else:
                                    M = collectionDate[1]
                    else:
                        M =  collectionDate[0]
                        Year =  collectionDate[1]
                        if "Jan" in collectionDateFormat or "JAN"in collectionDateFormat:
                            M = "01"
                        elif "Feb"  in collectionDateFormat or "FEB"in collectionDateFormat:
                            M = "02"
                        elif "Mar" in  collectionDateFormat or "MAR"in  collectionDateFormat:
                            M = '03'
                        elif "Apr" in  collectionDateFormat or "APR" in  collectionDateFormat:
                            M = '04'
                        elif "May" in  collectionDateFormat or "MAY" in  collectionDateFormat:
                            M = '05'
                        elif "Jun" in  collectionDateFormat or "JUN" in  collectionDateFormat:
                            M = '06'
                        elif "Jul" in  collectionDateFormat or "JUL"in  collectionDateFormat:
                            M = '07'
                        elif "Aug" in  collectionDateFormat or "AUG"in  collectionDateFormat:
                            M = '08'
                        elif "Sep" in  collectionDateFormat or "SEP" in collectionDateFormat:
                            M = '09'
                        elif "Oct" in  collectionDateFormat or "OCT"in  collectionDateFormat:
                            M = '10'
                        elif "Nov" in  collectionDateFormat or "NOV" in  collectionDateFormat:
                            M = '11'
                        elif "Dec" in  collectionDateFormat or "DEC" in collectionDateFormat:
                            M = '12'
                        else:
                            if collectionDate[0].startswith("0"):
                                M = collectionDate[0]
                            else:
                                if len(collectionDate[0]) ==1:
                                    M = "0"+ collectionDate[0]
                                elif collectionDate[0] == "0":
                                    M = "01"
                                else:
                                    M = collectionDate[0]
                else:
                    Year = collectionDate[0]
                    M = "01"
                    Day = "01"        
            else:
                Year = collection_date.split("/")[-1]
                M = "01"
                Day = "01"        

            collection_Date_Format = Year + "-" +M +"-" + Day       

        file_mysql.write(str(Nucleotide_ID_num)+"\t"+locus+"\t"+accession+"\t"+version+"\t"+isRefseq+"\t"+spciesname+"\t"+TaxonID+"\t"+definition+"\t"+Type+"\t"+isolate+"\t"+strain+"\t"+host+"\t"+HostFormat+"\t"+isolationSource+"\t"+length+"\t"+collection_date+"\t"+year+"\t"+country+"\t"+CountryFormat+"\t"+LocationFormat+"\t"+submit_date+"\t"+SubmitLocation+"\t"+collection_Date_Format+"\t"+submit_Date_Format+"\t"+genome_b+"\t"+SARS_CoV_2_b+"\t \tNCBI\n")
        file_out.write(locus+"\t"+spciesname+"\t"+definition+"\t"+isolate+"\t"+strain+"\t"+host+"\t"+isolationSource+"\t"+length+"\t"+collection_date+"\t"+year+"\t"+country+"\t"+submit_date+"\t"+collection_Date_Format+"\t"+submit_Date_Format+"\t"+genome+"\t"+SARS_CoV_2+"\n")
    for each_locus in genome_gff:
        os.system("cp %s/Coronaviridae_genbank_%s/%s.fasta  %s/Coronaviridae_genbank_%s/fasta/"%(dir,date,each_locus,dir,date))


def main():
    count_number = 0
    Nucleotide_ID_num = int(num)
    for eachline in file1:
        count_number +=1
        Nucleotide_ID_num +=1
        line = eachline.strip().split(",")
        locus = line[0]
        version = line[1]
        print ("正在下载第%s个数据，locus:%s"%(count_number,locus))
        filename = "%s/Coronaviridae_genbank_%s/%s.gbk"%(dir,date,locus)
        #filename = "/home/sxq/NCBI/data/download_data/mysql/nucleotide/Coronaviridae_genbank/%s.gbk"%(locus)
        if os.path.isfile(filename):
            gi_handle = Entrez.esummary(db="nucleotide", id=version)
            record = Entrez.read(gi_handle)
            GI = record[0]["Id"]
            extract_info(filename,version,locus,Nucleotide_ID_num,paper_IN_num,GI)
        else:
            gi_handle = Entrez.esummary(db="nucleotide", id=version)
            record = Entrez.read(gi_handle)
            GI = record[0]["Id"]
            net_handle = Entrez.efetch(db="nucleotide",id=GI,rettype="gb", retmode="text")
            out_handle = open(filename, "w")
            out_handle.write(net_handle.read())
            out_handle.close()
            net_handle.close()
            extract_info(filename,version,locus,Nucleotide_ID_num,paper_IN_num,GI)
    print("数据更新现在结束,一共更新%s条数据=============================================================================\n"%(count_number))    

if __name__ == '__main__':
    main()
    file_mysql.close()
    file_out.close()
    file_mysql_pmid.close()



os.chdir("%s/Coronaviridae_genbank_%s/"%(dir,date))
os.system("mkdir fasta_file gbk_file")
os.system("cp *.gbk /home/sxq/NCBI/data/download_data/mysql/nucleotide/Coronaviridae_genbank")
#os.system("tar -zhcvf Coronaviridae_genbank_fasta_%s.tar.gz  fasta/"%(date))

#os.system("tar -zhcvPf Coronaviridae_genbank_gff_%s.tar.gz  gffs/"%(date))
#os.system("rm -r %s/Coronaviridae_genbank_%s/fasta/"%(dir,date))
#os.system("rm -r %s/Coronaviridae_genbank_%s/gffs/"%(dir,date))
#os.system("mv %s/Coronaviridae_genbank_%s/Coronaviridae_genbank_gff_%s.tar.gz %s/gyc_%s"%(dir,date,date,dir,date))
#os.system("mv %s/Coronaviridae_genbank_%s/Coronaviridae_genbank_fasta_%s.tar.gz %s/gyc_%s"%(dir,date,date,dir,date))
os.system("mv %s/Coronaviridae_genbank_%s/fasta/ %s/gyc_%s"%(dir,date,dir,date))
os.system("mv %s/Coronaviridae_genbank_%s/gffs/ %s/gyc_%s"%(dir,date,dir,date))

os.system("mv *.fasta fasta_file")
os.system("mv *.gbk gbk_file")
os.system("tail -1 /home/sxq/NCBI/data/download_data/mysql/nucleotide/add_nucleotide_info_%s.txt|awk '{print $1}' >>/home/sxq/NCBI/data/download_data/mysql/nucleotide/config/number_nucleotide"%(date))
#os.system("tail -1 /home/sxq/NCBI/data/download_data/mysql/nucleotide/paper/add_paper_Nucle_info_%s.txt|awk '{print $1}'  >>/home/sxq/NCBI/data/download_data/mysql/nucleotide/config/number_paper"%(date))
#os.system("cat /home/sxq/NCBI/data/download_data/mysql/nucleotide/add_nucleotide_info_%s.txt >>/home/sxq/NCBI/data/download_data/mysql/nucleotide/file/nucleotide_info.txt"%(date))
#os.system("cat /home/sxq/NCBI/data/download_data/mysql/nucleotide/paper/add_paper_Nucle_info_%s.txt >>/home/sxq/NCBI/data/download_data/mysql/nucleotide/file/Paper_info.txt"%(date))
