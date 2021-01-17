import urllib.request
from helperfun import subprocess_cmd, getsequence
import pandas as pd
import os
import requests as r
import subprocess
import sys

aminoacids = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Asx': 'B', 'Cys': 'C', 'Glu': 'E', 'Gln': 'Q', 'Glx': 'Z', 'Gly': 'G', 'His': 'H',
              'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V', }

def ddgcalcs(pdb, aasub, gene):
    # parse amino acid substitution -> wild_position_mutant
    wild = aminoacids.get(aasub[:3])
    mutant = aminoacids.get(aasub[-3:])
    position = int(aasub[3:-3])
    # retrieve pdb file
    urllib.request.urlretrieve('http://files.rcsb.org/download/'+pdb+'.pdb', 'prediction/tmp/'+pdb+'.pdb')

    # read pdb file to determine where mutation is on the pdb (chain and nucleotide number)
    dbref = []
    chain = 'A'
    shift = 0
    for line in open('prediction/tmp/'+pdb+'.pdb'):
        listval = line.split()
        if listval[0] == 'DBREF':
            if len(listval) == 10:
                # chain, pos1, uniprot, protein, pos2
                dbref.append(
                    [listval[2], listval[3], listval[6], listval[7], listval[8]])
            else:
                # pos1, uniprot, protein, pos2
                dbref.append([listval[2], listval[5], listval[6], listval[7]])
    for mol in dbref:
        if len(mol) == 5:
            if gene in mol[3]:
                shift = int(mol[4]) - int(mol[1])
                chain = mol[0]
                uniprotcode = mol[2]
            else:
                secondchain = mol[0]
                secounduniprotcode = mol[2]
        elif len(mol) < 5 and gene in mol[2]:
            shift = int(mol[3]) - int(mol[0])
            uniprotcode = mol[1]

    # SAAMBE-3D calculation
    saambe_out = subprocess_cmd('source /Users/cameosameshima/opt/anaconda3/etc/profile.d/conda.sh; conda activate py2;\
        cd prediction/saambe;\
        python Mutation_pred.py -i ../tmp/'+pdb+'.pdb -c '+chain+' -r '+str(position - shift)+' -w '+wild+' -m '+mutant+' -d 1').decode("utf-8")
    #saambe_eff = subprocess_cmd('source /Users/cameosameshima/opt/anaconda3/etc/profile.d/conda.sh; conda activate py2;\
    #    python Mutation_pred.py -i tmp/'+pdb+'.pdb -c '+chain+' -r '+str(position - shift)+' -w '+wild+' -m '+mutant+' -d 0').decode("utf-8")
    saambe_val = saambe_out.split("\n")[1]
    if saambe_val:
        saambe_eff = "Decrease in stability" if float(saambe_val) > 0.0 else "Increase in stability"
    else:
        saambe_val = saambe_eff = 'N/A'

    # imut2.0 struc calculation
    subprocess_cmd('source /Users/cameosameshima/opt/anaconda3/etc/profile.d/conda.sh; conda activate py2;\
        cd prediction/imutant;\
        ./mkdssp -i ../tmp/'+pdb+'.pdb -o ../tmp/'+pdb+'.dssp')
    imut2_out1 = subprocess_cmd('source /Users/cameosameshima/opt/anaconda3/etc/profile.d/conda.sh; conda activate py2;\
        cd prediction/imutant;\
        python -O  I-Mutant2.0.py -pdbv ../tmp/'+pdb+'.pdb ../tmp/'+pdb+'.dssp '+chain+' '+str(position - shift)+' '+mutant)
    imut2_out2 = subprocess_cmd('source /Users/cameosameshima/opt/anaconda3/etc/profile.d/conda.sh; conda activate py2;\
         cd prediction/imutant;\
        python -O  I-Mutant2.0.py -pdb ../tmp/'+pdb+'.pdb ../tmp/'+pdb+'.dssp '+chain+' '+str(position - shift)+' '+mutant)
    if "I-Mutant" in imut2_out1.decode("utf-8"):
        imut2_val = (imut2_out1.decode("utf-8").split("RSA")
                     [1].split("WT")[0]).split()[3]
        imut2_eff = ((imut2_out2.decode("utf-8").split("RSA")
                      [1].split("WT")[0]).split()[3])+' in stability'
    else:
        imut2_val = imut2_eff = 'N/A'

    # imut2.0 seq calculation
    sequence = getsequence(uniprotcode)
    with open("prediction/tmp/sequence.seq", "w") as text_file:
        text_file.write(sequence)

    imut2_seq_out1 = subprocess_cmd('source /Users/cameosameshima/opt/anaconda3/etc/profile.d/conda.sh; conda activate py2;\
        cd prediction/imutant;\
        python -O  I-Mutant2.0.py -seqv ../tmp/sequence.seq '+str(position)+' '+mutant)
    imut2_seq_out2 = subprocess_cmd('source /Users/cameosameshima/opt/anaconda3/etc/profile.d/conda.sh; conda activate py2;\
        cd prediction/imutant;\
        python -O  I-Mutant2.0.py -seq ../tmp/sequence.seq '+str(position)+' '+mutant)
    if "I-Mutant" in imut2_seq_out1.decode("utf-8"):
        imut2_seq_val = (imut2_seq_out1.decode("utf-8").split("pH    T")
                         [1].split("WT")[0]).split()[3]
        imut2_seq_eff = ((imut2_seq_out2.decode("utf-8").split("pH    T")
                          [1].split("WT")[0]).split()[3])+' in stability'
    else:
        imut2_seq_val = imut2_seq_eff = 'N/A'

    # UEP calculation
    os.system('cd prediction/uep; python3 UEP.py --pdb=../tmp/'+pdb +'.pdb --interface='+chain+','+secondchain)
    uep_out = pd.read_csv('prediction/tmp/'+pdb+'_UEP_'+chain+'_'+secondchain+'.csv')
    location = uep_out.loc[uep_out['Unnamed: 0'] == chain+'_' + str(position - shift)+'_'+aasub[:3].upper(), aasub[-3:].upper()]
    uep_val = "N/A" if location.empty else location.values[0]
    if uep_val != "N/A":
        uep_eff = "Improved binding affinity" if float(uep_val) < 0 else "Improved binding affinity" 
    else:
        uep_val = uep_eff = "N/A"

    # panda calculation
    secondseq = "'" + getsequence(secounduniprotcode) + "'"
    mutseq = "'" + sequence[:position-1] + mutant + sequence[position:] + "'"
    sequence = "'" + sequence + "'"

    panda_output = subprocess_cmd('cd prediction/panda;\
        (echo "from panda import *"; echo "print(predict_affinity('+secondseq+','+sequence+','+secondseq+','+mutseq+'))") | python')
    panda_val = panda_output.decode("utf-8")[1:-1]
    if panda_val:
        panda_eff = "Increasing affinity" if float(panda_val) > 0.0 else "Decreasing affinity"
    else:
        panda_val = panda_eff = 'N/A'

    # get rid of files
    os.remove('prediction/tmp/'+pdb+'.pdb')
    os.remove('prediction/tmp/'+pdb+'.dssp')
    os.remove('prediction/tmp/sequence.seq')
    os.remove('prediction/tmp/'+pdb+'_UEP_'+chain+'_'+secondchain+'.csv')

    return [["SAAMBE-3D",saambe_val, saambe_eff], ["I-Mutant2.0 Structure", imut2_val, imut2_eff], ["I-Mutant2.0 Sequence",imut2_seq_val, imut2_seq_eff], ["PANDA", panda_val, panda_eff], ["UEP",uep_val, uep_eff]], chain