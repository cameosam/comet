import urllib.request
from helperfun import *
import pandas as pd
import os
import requests as r
import subprocess
import sys
from os import path
from prediction.panda.panda import *

condapath = 'C:\\Users\\CameoSameshima\\Desktop\\thesis\\comet'
aminoacids = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Asx': 'B', 'Cys': 'C', 'Glu': 'E', 'Gln': 'Q', 'Glx': 'Z', 'Gly': 'G', 'His': 'H',
              'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V', }

def ddgcalcs(pdb, aasub, gene, protselect):
    # parse amino acid substitution -> wild_position_mutant
    wild = aminoacids.get(aasub[:3]) 
    mutant = aminoacids.get(aasub[-3:]) 
    position = int(aasub[3:-3])

    # retrieve sequence
    uniprotcode, interactors = get_uniprot_code_ints(gene)
    chain = "*"
    secondchain = "*"
    currentprot = "N/A"
    otherprot = "N/A"
    seconduniprotcode = "N/A"
    if pdb != "N/A":
        # retrieve pdb file
        urllib.request.urlretrieve('http://files.rcsb.org/download/'+pdb+'.pdb', 'prediction/tmp/'+pdb+'.pdb')

        # read pdb file to determine where mutation is on the pdb (chain and nucleotide number)
        dbref = []
        otherprot = []
        for line in open('prediction/tmp/'+pdb+'.pdb'):
            listval = line.split()
            if listval[0] == 'DBREF':
                if len(listval) == 10 and listval[5] != 'PDB':
                    # chain, pos1, uniprotcode, protein, pos2
                    dbref.append([listval[2], listval[3], listval[6], listval[7], listval[8]])
                else:
                    # pos1, uniprotcode, protein, pos2
                    dbref.append([listval[2], listval[5], listval[6], listval[7]])
        for mol in dbref:
            if len(mol) == 5:
                if gene in mol[3] or uniprotcode in mol[2]:
                    shift = int(mol[4]) - int(mol[1])
                    chain = mol[0]
                    uniprotcode = mol[2]
                else:
                    if mol[3] not in otherprot:
                        otherprot.append(mol[3])
                    if protselect != "N/A":
                        if mol[3] == protselect:
                            currentprot = mol[3]
                            secondchain = mol[0]
                            seconduniprotcode = mol[2]
                        else:
                            currentprot_temp = mol[3]
                            secondchain_temp = mol[0]
                            seconduniprotcode_temp = mol[2]
                    else:
                        currentprot = mol[3]
                        secondchain = mol[0]
                        seconduniprotcode = mol[2]
            elif len(mol) < 5 and (gene in mol[2] or uniprotcode in mol[1]):
                shift = int(mol[3]) - int(mol[0])
                uniprotcode = mol[1]
        if protselect != "N/A" and seconduniprotcode == "N/A":
            currentprot = currentprot_temp
            secondchain = secondchain_temp
            seconduniprotcode = seconduniprotcode_temp
        if interactors != [''] and otherprot == []:
            otherprot, seconduniprotcode = get_uniprot_names(interactors, protselect)
    else:
        seconduniprotcode = "N/A"

    sequence = getsequence(uniprotcode)

    if pdb != "N/A" and chain != "*":
        # SAAMBE-3D calculation
        saambe_out = subprocess_cmd('conda activate py2 && cd prediction/saambe && python Mutation_pred.py -i ../tmp/'+pdb+'.pdb -c '+chain+' -r '+str(position - shift)+' -w '+wild+' -m '+mutant+' -d 1').decode("utf-8")
        saambe_val = saambe_out.split("\n")[1] if "\n" in saambe_out else 'N/A'
        saambe_eff = muteffect(saambe_val,True) if saambe_val != 'N/A' and muteffect(saambe_val,True) else 'N/A'

        # imut2.0 struc calculation 
        # DSSP dependency not working
        # subprocess_cmd('cd prediction/imutant && mkdssp.exe -i ../tmp/'+pdb+'.pdb -o ../tmp/'+pdb+'.dssp')
        # imut2_out = subprocess_cmd('conda activate py2 && cd prediction/imutant && python -O  I-Mutant2.0.py -pdbv ../tmp/'+pdb+'.pdb ../tmp/'+pdb+'.dssp '+chain+' '+str(position - shift)+' '+mutant)
        # imut2_val = (imut2_out.decode("utf-8").split("RSA")[1].split("WT")[0]).split()[3] if "I-Mutant" in imut2_out.decode("utf-8") else 'N/A'
        # imut2_eff = muteffect(imut2_val,False) if imut2_val != 'N/A' and muteffect(imut2_val,False) else 'N/A'
        # get rid of files
        # os.remove('prediction/tmp/'+pdb+'.dssp')
        imut2_val = "Under construction"
        imut2_eff = "Under construction"
        
        if secondchain != "*":
            # UEP calculation
            os.system('cd prediction/uep && python UEP.py --pdb=../tmp/'+pdb +'.pdb --interface='+chain+','+secondchain)
            uep_file_path = 'prediction/tmp/'+pdb+'_UEP_'+chain+'_'+secondchain+'.csv'
            if path.exists(uep_file_path):
                uep_out = pd.read_csv(uep_file_path)
                location = uep_out.loc[uep_out['Unnamed: 0'] == chain+'_' + str(position - shift)+'_'+aasub[:3].upper(), aasub[-3:].upper()]
                uep_val = "N/A" if location.empty else location.values[0]
                uep_eff = muteffect(uep_val,True) if uep_val != "N/A" and muteffect(uep_val,True) else 'N/A'
                os.remove('prediction/tmp/'+pdb+'_UEP_'+chain+'_'+secondchain+'.csv')
            else:
                uep_val = uep_eff = 'N/A'
        else:
            uep_val = uep_eff = 'N/A'
    
        # get rid of files
        os.remove('prediction/tmp/'+pdb+'.pdb')
    else:
        saambe_val = saambe_eff = imut2_val = imut2_eff = uep_val = uep_eff = 'N/A'

    if seconduniprotcode != "N/A" and sequence != "N/A":
        # panda calculation
        mutseq = "'" + sequence[:position-1] + mutant + sequence[position:] + "'"
        sequence = "'" + sequence[:position-1] + wild + sequence[position:] + "'"
        secondseq = "'" + getsequence(seconduniprotcode) + "'"
        
        os.chdir(condapath+"\\prediction\\panda")
        
        panda_val = predict_affinity(secondseq,sequence,secondseq,mutseq)
        panda_val = panda_val[0]
        panda_eff = muteffect(panda_val,False) if panda_val != 'N/A' and muteffect(panda_val,False) else 'N/A'
        os.chdir(condapath)
        # print(os.getcwd())
    else:
        panda_val = panda_eff = 'N/A'

    # imut2.0 seq calculation
    if sequence != "N/A":
        with open(condapath+"\\prediction\\tmp\\sequence.seq", "w") as text_file:
            text_file.write(sequence)
        imut2_seq_out = subprocess_cmd('conda activate py2 && cd prediction/imutant && python -O  I-Mutant2.0.py -seqv ../tmp/sequence.seq '+str(position)+' '+mutant)
        if "I-Mutant" in imut2_seq_out.decode("utf-8"):
            imut2_seq_val = (imut2_seq_out.decode("utf-8").split("pH    T")
                                [1].split("WT")[0]).split()[3]
        else:
            imut2_seq_val ='N/A'
        imut2_seq_eff = muteffect(imut2_seq_val,False) if  imut2_seq_val != 'N/A' and muteffect(imut2_seq_val,False) else 'N/A'
        os.remove(condapath+'\\prediction\\tmp\\sequence.seq')
    else:
        imut2_seq_val = imut2_seq_eff = 'N/A'

    # Consensus
    increase = len(list(filter(lambda x: x == "Increase in stability", [saambe_eff, imut2_eff, imut2_seq_eff, panda_eff, uep_eff])))
    decrease = len(list(filter(lambda x: x == "Decrease in stability", [saambe_eff, imut2_eff, imut2_seq_eff, panda_eff, uep_eff])))
    if increase > decrease:
        consensus = "Increase in stability"
    elif decrease > increase:
        consensus = "Decrease in stability"
    else:
        consensus = "No change in stability"

    return [["SAAMBE-3D",saambe_val, saambe_eff], ["I-Mutant2.0 Structure", imut2_val, imut2_eff], ["I-Mutant2.0 Sequence",imut2_seq_val, imut2_seq_eff], ["PANDA", panda_val, panda_eff], ["UEP",uep_val, uep_eff],["Stability Consensus","",consensus]], [chain, secondchain, otherprot, currentprot]