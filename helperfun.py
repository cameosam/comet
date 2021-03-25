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
import pandas as pd

Entrez.email = "casamesh@lakeheadu.ca"

def subprocess_cmd(command):
    process = subprocess.Popen(command,stdout=subprocess.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()
    return (proc_stdout)

def parse_df(filename, chrom):
    rslist = []
    rs_df = pd.read_pickle("./uploads/"+filename+".pkl")
    chr_df = pd.read_pickle("./miss_chr1.pkl")
    for i in range(len(rs_df)):
        if rs_df.at[i,'rsid'] in chr_df[1].values:
            index_val = chr_df[chr_df[1]==rs_df.at[i,'rsid']].index.values
            rslist.append(list(rs_df.iloc[i,:].append(chr_df[2][index_val])))
    return rslist

def is_missense(snp_info):
    if snp_info[1] in ["1","2","3","4","5","6","7","8"]:
        if snp_info[3] != "--":
            chr_df = pd.read_pickle("./missense-db/miss_chr"+snp_info[1]+".pkl")
            if snp_info[0] in chr_df[1].values:
                return (chr_df.loc[chr_df[1] == snp_info[0], 2].item())
            else:
                return False
        else:
            return False       
    else:
        return False
    
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

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() == 'txt' and filename != ''

def get_clinlitvar(rsid):
    handle = Entrez.esearch(db="clinvar", term=rsid)
    record = Entrez.read(handle)
    conditions = []
    if len(record["IdList"]) > 0:
        handle = Entrez.esummary(db="clinvar", id=record["IdList"][0])
        record = Entrez.read(handle)
        for i in range(0,len(record['DocumentSummarySet']['DocumentSummary'][0]['trait_set'])):
            trait = record['DocumentSummarySet']['DocumentSummary'][0]['trait_set'][i]['trait_name']
            if "not" not in trait and trait not in conditions:
                conditions.append(trait)
    
    response = requests.get("https://www.ncbi.nlm.nih.gov/research/bionlp/litvar/api/v1/entity/litvar/"+rsid+"%23%23")
    if response.text != "":
        diseases = response.json()['diseases']
        disease_name = []
        disease_freq = []
        for key, value in diseases.items():
            disease_name.append(key.replace(",",";"))
            disease_freq.append(value)
    else:
        disease_name = disease_freq = "N/A"
    if disease_name == []:
        disease_name = disease_freq = "N/A"
    return [conditions, disease_name, disease_freq]

def getsnpinfo(rsid):
    handle = Entrez.esummary(db="snp", term="[snp_id]", id=rsid[2:])
    record = Entrez.read(handle)['DocumentSummarySet']['DocumentSummary'][0]

    # gene, chromosome, amino acid substitution, clinical significance
    gene_name = record['GENES'][0]['NAME'] if (record['GENES']) else "N/A"
    chromosome = record['CHR']
    docsum = record['DOCSUM']
    clinical = record['CLINICAL_SIGNIFICANCE'] if record['CLINICAL_SIGNIFICANCE'] !=  '' else "N/A"
    clinical = clinical.replace("-"," ").split(",")
    print(record['CLINICAL_SIGNIFICANCE'], clinical)
    # substitutions
    nuclist = []
    aalist = []
    first_aa = []
    entries = docsum.split(",")
    for i in entries:
        if any(x in i for x in ["AC","NC","NG","NT","NW","NZ"]) and ">" in i[-3:] and i[-3:] not in nuclist:
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
        zipped_nuclist = sorted(list(zip(nuclist, nuccount)),key=lambda x: (x[1]), reverse=True)
        zipped_aalist = sorted(list(zip(aalist, aacount)),key=lambda x: (x[1]), reverse=True)
        first_aa = zipped_aalist[0][0]
        sorted_nuclist = [[i[0] for i in zipped_nuclist],[i[1] for i in zipped_nuclist]]
        sorted_aalist = [[i[0] for i in zipped_aalist],[i[1] for i in zipped_aalist]]
    else:
        zipped_nuclist = sorted(list(zip(nuclist, nuccount)),key=lambda x: (x[1]), reverse=True)
        sorted_nuclist = [[i[0] for i in zipped_nuclist],[i[1] for i in zipped_nuclist]]
        sorted_aalist = first_aa = 'N/A'


    # minor allele frequency
    maf = record['GLOBAL_MAFS']
    freq_kg = 'N/A'
    freq_hm = 'N/A'
    for study in maf:
        if '1000Genomes' in study.values():
            freq_kg = study['FREQ'].partition("/")[0]
        if 'HapMap' in study.values():
            freq_hm = study['FREQ'].partition("/")[0]

    return gene_name, chromosome, freq_kg, freq_hm, clinical, sorted_nuclist, sorted_aalist, first_aa

def getsequence(uniprotcode):
    baseUrl = "http://www.uniprot.org/uniprot/"
    currentUrl = baseUrl+uniprotcode+".fasta"
    response = requests.post(currentUrl)
    cData = ''.join(response.text)
    Seq = StringIO(cData)
    pSeq = list(SeqIO.parse(Seq, 'fasta'))
    return str(pSeq[0].seq)

def get_uniprot_code_ints(gene):
    baseUrl = "http://www.uniprot.org/uniprot/"
    # query = "?query=reviewed:yes+AND+gene_exact:"+gene+"+organism:9606&format=list"
    query = "?query=reviewed:yes+AND+gene_exact:"+gene+"+organism:9606&columns=id,interactor&format=tab"  
    currentUrl = baseUrl+query
    response = requests.post(currentUrl)
    response_vals = response.text.split("\n")[1]
    code = response_vals.split("\t")[0]
    interactors = response_vals.split("\t")[1]
    interactors_list = []
    if ";" in interactors:
        for i in interactors.split("; "):
            if i not in interactors_list:
                interactors_list.append(i)
    else:
        interactors_list.append(interactors)
    return code, interactors_list

def get_uniprot_names(codes, protselect):
    names = []
    uniprotcode = "N/A"
    for code in codes:
        baseUrl = "http://www.uniprot.org/uniprot/"
        query = "?query=reviewed:yes+id:"+code+"+organism:9606&columns=entry name&format=tab"  
        currentUrl = baseUrl+query
        response = requests.post(currentUrl)
        if "Entry name" in response.text:
            if response.text.split("\n")[1] == protselect:
                uniprotcode = code
            names.append(response.text.split("\n")[1])
    if uniprotcode == "N/A":
        uniprotcode = code
    return names, uniprotcode

def muteffect(ddg, posdes):
    if posdes == True:
        return "Decrease in stability" if float(ddg) > 0.0 else "Increase in stability"
    else:
        return "Decrease in stability" if float(ddg) < 0.0 else "Increase in stability"

def deletefiles(files):
    now = time.time()
    for filename in files:
            # automatic delete of upladed files after 24 hours
            if (now - os.stat(filename).st_mtime) > (24 * 60 * 60):
                command = "rm {0}".format(filename)
                subprocess.call(command, shell=True)

def comparegeno(genotype, sorted_nuclist, freq_kg, freq_hm):
    mutlist = []
    count = 0
    if len(genotype) > 1:
        for i in range(len(sorted_nuclist[0])):
            mutlist.append(sorted_nuclist[0][i][2])
        if genotype[0] in mutlist:
            count += 1
        if genotype[1] in mutlist:
             count += 1
    else:
        for i in range(len(sorted_nuclist[0])):
            mutlist.append(sorted_nuclist[0][i][2])
        if genotype in mutlist:
            count -= 1

    if freq_kg != 'N/A' and freq_hm != 'N/A':
        avg_freq =  round((float(freq_kg[3:]) + float(freq_hm[3:]))/2*100,2)
    elif freq_kg != 'N/A':
        avg_freq = round(float(freq_kg[3:])*100,2)
    elif freq_hm != 'N/A':
        avg_freq = round(float(freq_hm[3:])*100,2)
    else:
        avg_freq = -1
    # same = False
    # if avg_freq != -1:
    #     if freq_kg[0] in genotype or freq_hm[0] in genotype:
    #         same = True

    return [count, avg_freq]
