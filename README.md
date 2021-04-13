# COMET
Honours Thesis Project. <br>
COMET (Comparative mutation effect) is a web application that analyzes genotype data in terms of protein-protein interactions (PPIs) by calculating the change in binding affinity (∆∆G) 

# Dependencies and Installation
This application uses four ∆∆G predictors: <br>
SAAMBE 3D : Pahari, S., Li, G., Murthy, A. K., Liang, S., Fragoza, R., Yu, H., & Alexov, E. (2020). SAAMBE-3D: Predicting Effect of Mutations on Protein–Protein Interactions. International journal of molecular sciences, 21(7), 2563.
<a href="http://compbio.clemson.edu/saambe_webserver/">SAAMBE Website Link</a>

I-Mutant2.0: Capriotti, E., Fariselli, P., & Casadio, R. (2005). I-Mutant2.0: predicting stability changes upon mutation from the protein sequence or structure.
Nucleic acids research, 33, W306 - W310.
<a href="https://folding.biofold.org/i-mutant/i-mutant2.0.html">I-Mutant2.0 Website Link</a>

PANDA: Abbasi, W. A.; Abbas, S. A.; Andleeb, S. PANDA: Predicting the Change in Proteins Binding Affinity upon Mutations Using Sequence Information. ArXiv200908869 Cs Q-Bio Stat 2020.
<a href="https://github.com/wajidarshad/panda">PANDA GitHub Link</a>

UEP: Amengual-Rigo, P., Fernández-Recio, J., & Guallar, V. (2020). UEP: an open-source and fast classifier for predicting the impact of mutations in protein–protein complexes. Bioinformatics.
<a href="https://github.com/pepamengual/UEP">UEP GitHub Link</a>  
      
To use multiple Python versions, conda environments were created: <br>
py3: conda create -n py3 python=3.6 scikit-learn=0.20.3 biopython=1.77 flask=1.1.2 pandas=1.1.3 requests=2.25.1. <br>
Prody (1.10.8) and compress_pickle (1.2.0) (for UEP) were installed using pip. <br>
py2: conda create -n py2 python=2.7 numpy=1.16.6 py-xgboost=0.90 pandas=0.22.0.<br>

Two folders must be created: uploads and prediction. The uploads folder stores the uploaded genotype text and pickle files. The prediction folder has the ∆∆G predictors folders (imutant, panda, saambe, and uep) and a temporary folder (tmp) for the pdb and sequences files. The ∆∆G predictors' program files can be downloaded from the above links. 