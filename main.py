from flask import Flask, redirect, url_for, render_template, request, session, jsonify, flash
from werkzeug.utils import secure_filename
from Bio import Entrez
import itertools
import os
import json
import pandas as pd
import secrets
import time

from helperfun import *
# change ddgcalc_os to ddgcalc_win when using windows
# from ddgcalc_win import ddgcalcs
from ddgcalc_os import ddgcalcs


app = Flask(__name__)
# MAX 30 mb
app.config['MAX_CONTENT_LENGTH'] = 1024 * 1024 * 30 
app.config['UPLOAD_FOLDER'] = 'uploads/'
app.config['SECRET_KEY'] = "temporarysecretkey"
app.config['SESSION_PERMANENT'] = True
database = "snp-db/snpDB_v5.pkl"

@app.route("/")
@app.route("/home")
def index():
    uploadfiles = [os.path.join(app.config['UPLOAD_FOLDER'], filename) for filename in os.listdir( app.config['UPLOAD_FOLDER'])]
    tmpfiles = [os.path.join('prediction/tmp/', filename) for filename in os.listdir('prediction/tmp/')]
    deletefiles(uploadfiles)
    deletefiles(tmpfiles)
    return render_template("index.html")

@app.route("/about")
def about():
    return render_template("about.html")

@app.route("/endsession")
def endsession():
    if session.get("file") != None:
        if session["file"] != "N/A":
            os.remove(app.config['UPLOAD_FOLDER']+session["file"]+'.pkl')
            os.remove(app.config['UPLOAD_FOLDER']+session["file"]+'unfiltered.pkl')
            [session.pop(key) for key in list(session.keys())]
            flash('Thanks for using COMET. All session data and files have been deleted.')
        else:
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
                filename = filename[:-4] + secrets.token_urlsafe(16)  
                file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
                session["file"] = filename
                rs_list = []
                for line in open(app.config['UPLOAD_FOLDER']+filename):
                    if '#' not in line and 'i' not in line and "rsid" not in line:
                        # 23andMe
                        if 'rs' in line and len(re.split(r'\t', line.rstrip('\n'))) == 4:
                            temp = re.split(r'\t', line.rstrip('\n'))
                            if temp[3] != "--":
                                rs_list.append(temp)
                        # Ancestry
                        elif 'rs' in line and len(re.split(r'\t', line.rstrip('\n'))) == 5:
                            temp = re.split(r'\t', line.rstrip('\n'))
                            if temp[3] != "0" and temp[3] != "I" and temp[3] != "D":
                                temp[3] = temp[3] + temp[4]
                                del temp[-1]
                                rs_list.append(temp)
                        # MyHeritage
                        elif 'rs' in line and len(re.split(',', line.rstrip('\n'))) == 4:
                            temp = re.split('","', line.rstrip('\n'))
                            if temp[3][:-1] != "--":
                                temp[0] = temp[0][1:]
                                temp[3] = temp[3][:-1]
                                rs_list.append(temp)
                if len(rs_list) > 0:
                    rs_df = pd.DataFrame(rs_list,columns=['rsid','chr','pos','geno'])
                    rs_df.replace({'chr': {"23": "X", "24": "Y", "25": "Y", "26" : "MT"}})
                    db = pd.read_pickle(database)
                    db = db.drop(['chr', 'pos'], axis=1)
                    rs_df = rs_df.join(db.set_index('rsid'), on='rsid')
                    rs_df.to_pickle(app.config['UPLOAD_FOLDER']+filename+"unfiltered.pkl")
                    rs_df = rs_df.dropna()
                    rs_df = rs_df.reset_index(drop=True)
                    os.remove(app.config['UPLOAD_FOLDER']+filename)
                    rs_df.to_pickle(app.config['UPLOAD_FOLDER']+filename+".pkl")
                    return redirect(url_for("select"))
                else:
                    os.remove(app.config['UPLOAD_FOLDER']+filename)
                    return render_template("input.html", error = "Error: File could not be parsed. Please check that this is a raw genotype file.")
            else:
                return render_template("input.html", error = "Error: Incorrect file type. Please upload a .txt file.")
    else:
        return render_template("input.html")

@app.route("/select", methods=["POST", "GET"])
def select():
    if request.method == "POST": 
        if 'turnoff' in request.form:
            session["turnoffdb"] = True
        if 'turnon' in request.form:
            session.pop("turnoffdb")
        if 'chrom' in request.form:
            session["chrom"] = request.form['chrom']
    if "file" in session and session["file"] != "N/A":
        if "turnoffdb" in session:
            return render_template("select.html", list=list, geno=True, dboff = 1)
        return render_template("select.html", list=list, geno=True, dboff = 0)
    else:
        session["file"] = "N/A"
        return render_template("select.html", list=list, geno=False, dboff = -1)

@app.route('/showtable')
def index_get_data():
    chrom = session["chrom"] if session.get("chrom") != None else str(1)
    if session["file"] != "N/A":
        filename = session["file"]
        if "turnoffdb" in session:
            rs_df = pd.read_pickle(app.config['UPLOAD_FOLDER']+filename+"unfiltered.pkl")
            cols = ['rsid', 'chromosome', 'position', 'genotype', 'substitution', 'gene','freq']
        else:
            rs_df = pd.read_pickle(app.config['UPLOAD_FOLDER']+filename+".pkl")
            cols = ['rsid', 'chromosome', 'position', 'genotype', 'substitution', 'gene','freq']
        rs_df_chr = rs_df.loc[rs_df['chr'] == chrom].values.tolist() if chrom != "All" else rs_df.values.tolist()
    else:
        rs_df = pd.read_pickle(database)
        rs_df_chr = rs_df.loc[rs_df['chr'] == chrom].values.tolist() if chrom != "All" else rs_df.values.tolist()
        cols = ['rsid', 'chromosome', 'position', 'substitution', 'gene','freq']
    df = pd.DataFrame(rs_df_chr, columns=cols)    
    datatable = df.to_json(orient="table")
    jsontable = json.loads(datatable)
    return jsonify(jsontable)

@app.route('/nextsnp')
def next_snp():
    filename = session["file"]
    if "turnoffdb" in session:
        filename = filename+"unfiltered"
    filename = app.config['UPLOAD_FOLDER']+filename+".pkl" if filename != "N/A" else database
    curr_snp = session["rsid"]
    rs_df = pd.read_pickle(filename)
    index = rs_df[rs_df['rsid']==curr_snp].index.values.astype(int)[0]
    if index == len(rs_df)-1:
        index = 0
    session["rsid"] = rs_df.loc[rs_df.index[index+1], 'rsid']
    if 'geno' in rs_df.columns:
        session["genotype"] = rs_df.loc[rs_df.index[index+1], 'geno']
    return redirect(url_for("output"))

@app.route('/prevsnp')
def prev_snp():
    filename = session["file"]
    if "turnoffdb" in session:
        filename = filename+"unfiltered"
    filename = app.config['UPLOAD_FOLDER']+filename+".pkl" if filename != "N/A" else database
    curr_snp = session["rsid"]
    rs_df = pd.read_pickle(filename)
    index = rs_df[rs_df['rsid']==curr_snp].index.values.astype(int)[0]
    if index == 0:
        index = len(rs_df)-1
    session["rsid"] = rs_df.loc[rs_df.index[index-1], 'rsid']
    if 'geno' in rs_df.columns:
        session["genotype"] = rs_df.loc[rs_df.index[index-1], 'geno']
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
        rsid, genotype, gene_name, chromosome, pdbs, first_pdb, first_aa, freq_kg, freq_hm, clinical, sorted_nuclist, sorted_aalist, condition, summary= session[
                "output"]
        if 'ddgcalc' in request.form:
            # calculate the ddg values for the first time
            if first_aa != "N/A":
                protselect = "N/A"
                ddgresults, chain = ddgcalcs(first_pdb, first_aa, gene_name, protselect)
            else:
                chain = "*"
                ddgresults = [["N/A" for i in range(2)] for j in range(4)]
            return render_template("output.html", snp=rsid, genotype = genotype, gene=gene_name, chr=chromosome, pdb=pdbs, pdbselect=first_pdb, aa1=first_aa, ddgresults=ddgresults, freq1000g=freq_kg, freqhapmap=freq_hm, clin=clinical, sorted_nuclist=sorted_nuclist, sorted_aalist=sorted_aalist, condition=condition, summary=summary, chain=chain, zip=zip, len=len,float=float, round=round)
        elif 'pdbselect' in request.form:
            # recalculate the ddg values with a different pdb/protein
            pdbselect = request.form['pdbselect']
            protselect = request.form['protselect'] if 'protselect' in request.form else "N/A"
            if first_aa != "N/A":
                ddgresults, chain = ddgcalcs(pdbselect, first_aa, gene_name, protselect)
            return render_template("output.html", snp=rsid, genotype = genotype, gene=gene_name, chr=chromosome, pdb=pdbs, pdbselect=pdbselect, aa1=first_aa, ddgresults=ddgresults, freq1000g=freq_kg, freqhapmap=freq_hm, clin=clinical, sorted_nuclist=sorted_nuclist, sorted_aalist=sorted_aalist, condition=condition, summary=summary, chain=chain, zip=zip, len=len,float=float, round=round)
    if "rsid" in session:
        # results page (no ddg values)
        rsid = session["rsid"]
        genotype = session["genotype"]
        gene_name, chromosome, freq_kg, freq_hm, clinical, sorted_nuclist, sorted_aalist, first_aa = getsnpinfo(rsid) 
        summary = comparegeno(genotype, sorted_nuclist, freq_kg, freq_hm)
        condition = get_clinlitvar(rsid)
        pdbs = findpdb(gene_name)
        if pdbs != "N/A" and len(pdbs) > 0:
            first_pdb = pdbs[0]
        else:
            first_pdb = "N/A"
        session["output"] = [rsid, genotype, gene_name, chromosome, pdbs, first_pdb, first_aa, freq_kg, freq_hm, clinical, sorted_nuclist, sorted_aalist, condition, summary]
        return render_template("output.html", snp=rsid, genotype= genotype,gene=gene_name,
                               chr=chromosome, pdb=pdbs, pdbselect=first_pdb,
                               aa1=first_aa, freq1000g=freq_kg,
                               freqhapmap=freq_hm, clin=clinical,
                               sorted_nuclist=sorted_nuclist, sorted_aalist=sorted_aalist, condition=condition, summary=summary, chain = "*",
                               zip=zip, len=len, float=float, round=round)
    else:
        return redirect(url_for('index'))

if __name__ == "__main__":
    app.run(debug=True)