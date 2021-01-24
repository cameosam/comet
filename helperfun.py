import shlex
import subprocess
import time
import mmap
import re
from Bio import Entrez, SeqIO
import json
import os
import urllib.request
import requests
from io import StringIO

Entrez.email = "casamesh@lakeheadu.ca"

def subprocess_cmd(command):
    process = subprocess.Popen(command,stdout=subprocess.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()
    return (proc_stdout)

def parsefile(filename, chrom):
    rslist = []
    for line in open('uploads/'+filename):
        if '#' not in line and 'i' not in line:
            temp = re.split(r'\t', line.rstrip('\n'))
            if temp[1] == chrom and temp[3] != "--":
                rslist.append(temp)
    return rslist

def findpdb(gene):
    params = {
        "query": 
        {
            "type": "group", 
            "logical_operator": "and", 
            "nodes": [
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "operator": "exact_match",
                "value": "Homo sapiens",
                "attribute": "rcsb_entity_source_organism.taxonomy_lineage.name"
            }
        },
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "operator": "exact_match",
                "value": gene,
                "attribute":"rcsb_entity_source_organism.rcsb_gene_name.value"
            }
        }
    ]
    },
        "request_options": {
            "pager": {
                "start": 0,
                "rows": 100
            }
    },
    "return_type": "entry"
    }

    params_json = json.dumps(params)
    url = "https://search.rcsb.org/rcsbsearch/v1/query?json=%s"
    response = requests.get(url % params_json)

    if response.status_code == 200:
        pdbs = []
        count = response.json().get("total_count")
        for i in range(count):
            pdbs.append((response.json().get("result_set"))[i].get("identifier"))
        return (pdbs)
    else:
        return "N/A"

# def pdbfilter(pdbs, gene):
#     filteredpdb = []
#     print(pdbs)
#     for pdb in pdbs:
#         print(pdb)
#         urllib.request.urlretrieve('http://files.rcsb.org/download/'+pdb+'.pdb', 'prediction/tmp/'+pdb+'.pdb')
#         for line in open('prediction/tmp/'+pdb+'.pdb'):
#             listval = line.split()
#             if listval[0] == 'DBREF':
#                 if gene in listval[7] or gene in listval[6]:
#                     filteredpdb.append(pdb)
#         os.remove('prediction/tmp/'+pdb+'.pdb')
#     return filteredpdb

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() == 'txt'

def getsnpinfo(rsid):
    handle = Entrez.esummary(db="snp", term="[snp_id]", id=rsid[2:])
    record = Entrez.read(handle)['DocumentSummarySet']['DocumentSummary'][0]

    # gene, chromosome, amino acid substitution, clinical significance
    gene_name = record['GENES'][0]['NAME'] if (record['GENES']) else "N/A"
    chromosome = record['CHR']
    docsum = record['DOCSUM']
    clinical = record['CLINICAL_SIGNIFICANCE'] if record['CLINICAL_SIGNIFICANCE'] !=  '' else "N/A"

    # substitutions
    nuclist = []
    aalist = []
    first_aa = []
    entries = docsum.split(",")
    for i in entries:
        if ">" in i[-3:] and i[-3:] not in nuclist:
            nuclist.append(i[-3:])
        if "p." in i and (i.partition("p.")[2] not in aalist and len(i.partition("p.")[2]) < 10):
            aalist.append(i.partition("p.")[2])
    nuccount = []
    aacount = []
    for i in nuclist:
        nuccount.append(docsum.count(i))
    for i in aalist:
        aacount.append(docsum.count(i))
    if aacount:
        sorted_nuclist = sorted(list(zip(nuclist, nuccount)),key=lambda x: (x[1]), reverse=True)
        sorted_aalist = sorted(list(zip(aalist, aacount)),key=lambda x: (x[1]), reverse=True)
        subs = [sorted_nuclist, sorted_aalist]
    else:
        sorted_nuclist = sorted(list(zip(nuclist, nuccount)),key=lambda x: (x[1]), reverse=True)
        subs = [sorted_nuclist]
    #subs = sorted(list(zip(nuclist, nuccount, aalist, aacount)),key=lambda x: (x[1]), reverse=True) if aacount else sorted(list(zip(nuclist, nuccount)),key=lambda x: (x[1]), reverse=True)
    first_aa = sorted_aalist[0][0] if aacount else 'N/A'
    print(subs)
    # minor allele frequency
    maf = record['GLOBAL_MAFS']
    freq_kg = 'N/A'
    freq_hm = 'N/A'
    for study in maf:
        if '1000Genomes' in study.values():
            freq_kg = study['FREQ'].partition("/")[0]
        if 'HapMap' in study.values():
            freq_hm = study['FREQ'].partition("/")[0]

    return gene_name, chromosome, freq_kg, freq_hm, clinical.split(","), subs, first_aa

def getsequence(uniprotcode):
    baseUrl = "http://www.uniprot.org/uniprot/"
    currentUrl = baseUrl+uniprotcode+".fasta"
    response = requests.post(currentUrl)
    cData = ''.join(response.text)
    Seq = StringIO(cData)
    pSeq = list(SeqIO.parse(Seq, 'fasta'))
    return str(pSeq[0].seq)

def get_uniprot_code(gene):
    baseUrl = "http://www.uniprot.org/uniprot/"
    query = "?query=reviewed:yes+AND+gene_exact:"+gene+"+organism:9606&format=list"
    currentUrl = baseUrl+query
    response = requests.post(currentUrl)
    return response.text.split("\n")[0]

def muteffect(ddg, posdes):
    if posdes == True:
        return "Decrease in stability" if float(ddg) > 0.0 else "Increase in stability"
    else:
        return "Decrease in stability" if float(ddg) < 0.0 else "Increase in stability"
