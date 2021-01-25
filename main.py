from flask import Flask, redirect, url_for, render_template, request, session, jsonify, flash
from werkzeug.utils import secure_filename
from Bio import Entrez
import itertools
import os
import json
import pandas as pd
from helperfun import *
from ddgcalc import ddgcalcs

app = Flask(__name__)
app.secret_key = "secret_key"
app.config['UPLOAD_FOLDER'] = './uploads/'

@app.route("/")
@app.route("/home")
def index():
    return render_template("index.html")

@app.route("/about")
def about():
    return render_template("about.html")

@app.route("/endsession")
def endsession():
    if session.get("file") != None:
        os.remove('uploads/'+session["file"])
        [session.pop(key) for key in list(session.keys())]
        flash('Thanks for using COMET. All session data has been deleted.')
    else:
        flash('Nothing to be deleted!')
    return redirect(url_for('index'))

@app.route("/input", methods=["POST", "GET"])
def input():
    if request.method == "POST":
        if 'file' in request.files:
            file = request.files['file']
            if allowed_file(file.filename):
                filename = secure_filename(file.filename)
                file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
                session["file"] = filename
                return redirect(url_for("select"))
    else:
        return render_template("input.html")

@app.route("/select", methods=["POST", "GET"])
def select():
    if request.method == "POST" and 'chrom' in request.form:
        session["chrom"] = request.form['chrom']
    return render_template("select.html", list=list)

@app.route('/showtable')
def index_get_data():
    chrom = session["chrom"] if session.get("chrom") != None else str(1)
    filename = session["file"]
    rslist = parsefile(filename, chrom)
    cols = ['rsid', 'chromosome', 'position', 'genotype']
    df = pd.DataFrame(rslist, columns=cols)
    datatable = df.to_json(orient="table")
    jsontable = json.loads(datatable)
    return jsonify(jsontable)

@app.route('/nextsnp')
def next_snp():
    curr_snp = session["rsid"]
    filename = session["file"]
    chrom = session["chrom"] if session.get("chrom") != None else str(1)
    rslist = parsefile(filename, chrom)
    curr_index = [index for index, row in enumerate(rslist) if curr_snp in row]
    session["rsid"] = rslist[curr_index[0]+1][0]
    session["genotype"] = rslist[curr_index[0]+1][3]
    return redirect(url_for("output"))

@app.route('/prevsnp')
def prev_snp():
    curr_snp = session["rsid"]
    filename = session["file"]
    chrom = session["chrom"]
    rslist = parsefile(filename, chrom)
    curr_index = [index for index, row in enumerate(rslist) if curr_snp in row]
    session["rsid"] = rslist[curr_index[0]-1][0]
    session["genotype"] = rslist[curr_index[0]-1][3]
    return redirect(url_for("output"))

@app.route("/getchoice", methods=["POST"])
def testid():
    if request.method == "POST":
        session["rsid"] = request.args.get('rsid', '')
        session["genotype"] = request.args.get('genotype', '')
        return redirect(url_for("output"))

@app.route("/output", methods=["POST", "GET"])
def output():
    if request.method == "POST":
        rsid, genotype, gene_name, chromosome, pdbs, first_pdb, first_aa, freq_kg, freq_hm, clinical, sorted_nuclist, sorted_aalist = session[
                "output"]
        if 'ddgcalc' in request.form:
            if request.form['ddgcalc'] == 'mutated':
                if first_aa != "N/A":
                    ddgresults, chain = ddgcalcs(first_pdb, first_aa, gene_name, True)
                else:
                    chain = "*"
                    ddgresults = [["N/A" for i in range(2)] for j in range(4)]
                return render_template("output.html", snp=rsid, genotype = genotype, gene=gene_name, chr=chromosome, pdb=pdbs, pdbselect=first_pdb, aa1=first_aa, ddgresults=ddgresults, freq1000g=freq_kg, freqhapmap=freq_hm, clin=clinical, sorted_nuclist=sorted_nuclist, sorted_aalist=sorted_aalist, chain=chain, zip=zip, len=len)
            else:
                if first_aa != "N/A":
                    ddgresults, chain = ddgcalcs(first_pdb, first_aa, gene_name, False)    
                else:
                    chain = "*"
                    ddgresults = [["N/A" for i in range(2)] for j in range(4)]
                return render_template("output.html", snp=rsid, genotype = genotype, gene=gene_name, chr=chromosome, pdb=pdbs, pdbselect=first_pdb, aa1=first_aa, ddgresults=ddgresults, freq1000g=freq_kg, freqhapmap=freq_hm, clin=clinical, sorted_nuclist=sorted_nuclist, sorted_aalist=sorted_aalist, chain=chain, zip=zip, len=len) 
        elif 'pdbselect' in request.form:
            pdbselect = request.form['pdbselect']
            if first_aa != "N/A":
                ddgresults, chain = ddgcalcs(pdbselect, first_aa, gene_name, True)
            return render_template("output.html", snp=rsid, genotype = genotype, gene=gene_name, chr=chromosome, pdb=pdbs, pdbselect=pdbselect, aa1=first_aa, ddgresults=ddgresults, freq1000g=freq_kg, freqhapmap=freq_hm, clin=clinical, sorted_nuclist=sorted_nuclist, sorted_aalist=sorted_aalist, chain=chain, zip=zip, len=len)
        
    if "rsid" in session:
        rsid = session["rsid"]
        genotype = session["genotype"]
        gene_name, chromosome, freq_kg, freq_hm, clinical, sorted_nuclist, sorted_aalist, first_aa = getsnpinfo(rsid) 
        pdbs = findpdb(gene_name)

        if pdbs != "N/A":
            first_pdb = pdbs[0] if len(pdbs) > 1 else pdbs
        else:
            first_pdb = "N/A"

        session["output"] = [rsid, genotype, gene_name, chromosome, pdbs, first_pdb, first_aa, freq_kg, freq_hm, clinical, sorted_nuclist, sorted_aalist]
 
        return render_template("output.html", snp=rsid, genotype= genotype,gene=gene_name,
                               chr=chromosome, pdb=pdbs, pdbselect=first_pdb,
                               aa1=first_aa, freq1000g=freq_kg,
                               freqhapmap=freq_hm, clin=clinical,
                               sorted_nuclist=sorted_nuclist, sorted_aalist=sorted_aalist, chain = "*",
                               zip=zip, len=len)
    else:
        return redirect(url_for('index'))

if __name__ == "__main__":
    app.run(debug=True)
