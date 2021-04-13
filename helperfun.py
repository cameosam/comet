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
            try:
                pdbs.append((response.json().get("result_set"))[i].get("identifier"))
            except IndexError:
                break
        return (pdbs)
    else:
        return "N/A"

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() == 'txt' and filename != ''

def get_clinlitvar(rsid):
    # clinvar conditions
    handle = Entrez.esearch(db="clinvar", term=rsid)
    record = Entrez.read(handle, validate=False)
    cond_freq = {}
    if len(record["IdList"]) > 0:
        handle = Entrez.esummary(db="clinvar", id=record["IdList"][0])
        record = Entrez.read(handle, validate=False)
        for i in range(0,len(record['DocumentSummarySet']['DocumentSummary'][0]['trait_set'])):
            trait = record['DocumentSummarySet']['DocumentSummary'][0]['trait_set'][i]['trait_name']
            trait = trait.replace(",",";")
            if (trait in cond_freq):
                cond_freq[trait] += 1
            elif ("not" not in trait):
                cond_freq[trait] = 1
    # litvar diseases
    response = requests.get("https://www.ncbi.nlm.nih.gov/research/bionlp/litvar/api/v1/entity/litvar/"+rsid+"%23%23")
    dis_freq = response.json()['diseases'] if response.text != "" else {}
    dis_freq = {key.replace(",",";") : value for key, value in dis_freq.items()}
    return [list(cond_freq.keys()), list(cond_freq.values()),  list(dis_freq.keys()), list(dis_freq.values())]

def getsnpinfo(rsid):
    handle = Entrez.esummary(db="snp", term="[snp_id]", id=rsid[2:])
    record = Entrez.read(handle, validate=False)['DocumentSummarySet']['DocumentSummary'][0]

    # gene, chromosome, amino acid substitution, clinical significance
    gene_name = record['GENES'][0]['NAME'] if (record['GENES']) else "N/A"
    chromosome = record['CHR']
    docsum = record['DOCSUM']
    clinical = record['CLINICAL_SIGNIFICANCE'] if record['CLINICAL_SIGNIFICANCE'] !=  '' else "N/A"
    clinical = clinical.replace("-"," ").split(",")
    # substitutions
    nuc_freq = {}
    aa_freq = {}
    entries = docsum.split(",")
    for i in entries:
        if any(x in i for x in ["AC","NC","NG","NT","NW","NZ"]) and ">" in i[-3:]:
            n_sub = i[-3:]
            if (n_sub in nuc_freq):
                nuc_freq[n_sub] += 1
            else:
                nuc_freq[n_sub] = 1
        if ("p." in i) and (len(i.partition("p.")[2]) < 10):
            aa_sub = i.partition("p.")[2]
            if (aa_sub in aa_freq):
                aa_freq[aa_sub] += 1
            else:
                aa_freq[aa_sub] = 1
    nuc_freq = dict(sorted(nuc_freq.items(), key=lambda item: item[1], reverse=True))
    aa_freq = dict(sorted(aa_freq.items(), key=lambda item: item[1], reverse=True))
    sorted_nuclist = [list(nuc_freq.keys()),list(nuc_freq.values())]
    sorted_aalist = [list(aa_freq.keys()),list(aa_freq.values())]
    try:
        first_aa = sorted_aalist[0][0]
    except IndexError:
        first_aa = "N/A"

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
    try:
        return str(pSeq[0].seq)
    except IndexError:
        return "N/A"

def get_uniprot_code_ints(gene):
    baseUrl = "http://www.uniprot.org/uniprot/"
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
    # separate genotype
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
    # calcaulte average nucleotide frequency
    if 'N/A' not in freq_kg and 'N/A' not in freq_hm:
        try:
            avg_freq = round((float(freq_kg[3:]) + float(freq_hm[3:]))/2*100,2)
        except ValueError:
            avg_freq = -1
    elif 'N/A' not in freq_kg:
        try:
            avg_freq = round(float(freq_kg[3:])*100,2)
        except ValueError:
            avg_freq = -1
    elif 'N/A' not in freq_hm:
        try:
            avg_freq = round(float(freq_hm[3:])*100,2)
        except ValueError:
            avg_freq = -1
    else:
        avg_freq = -1
    return [count, avg_freq]
