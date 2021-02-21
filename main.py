from flask import Flask, redirect, url_for, render_template, request, session, jsonify, flash
from werkzeug.utils import secure_filename
from Bio import Entrez
import itertools
import os
import json
import pandas as pd
from helperfun import *
from ddgcalc import ddgcalcs
import secrets

# from flask_session import Session
# from datetime import timedelta
# from flask_sqlalchemy import SQLAlchemy 

app = Flask(__name__)
# MAX 30 mb
app.config['MAX_CONTENT_LENGTH'] = 1024 * 1024 * 30 
app.config['UPLOAD_FOLDER'] = './uploads/'
app.config['SECRET_KEY'] = "supersecretkey_MUSTBECHANGED!"
app.config['SESSION_PERMANENT'] = True

# Add server-side sessionsß
# app.config['PERMANENT_SESSION_LIFETIME'] = timedelta(hours=5)
# app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///db.sqlite3'
# app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False 
# app.config['SESSION_TYPE'] = 'sqlalchemy'
# db = SQLAlchemy(app)
# app.config['SESSION_SQLALCHEMY'] = db
# sess = Session(app)
# db.create_all()

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
        os.remove('uploads/'+session["file"]+'.pkl')
        os.remove('uploads/'+session["file"])
        # remove pickles
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
                filename = filename + secrets.token_urlsafe(16)  
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
                            if temp[3] != "0" and temp[3] != "I":
                                temp[3] = temp[3] + temp[4]
                                del temp[-1]
                                rs_list.append(temp)
                if len(rs_list) > 0:
                    rs_df = pd.DataFrame(rs_list,columns=['rsid','chr','pos','geno'])
                    rs_df.replace({'chr': {"23": "X", "24": "Y", "25": "Y", "26" : "MT"}})
                    db = pd.read_pickle("./snp-db/snpDB_V1.pkl")
                    rs_df = rs_df.join(db.set_index('rsid'), on='rsid')
                    rs_df = rs_df.dropna()
                    rs_df = rs_df.reset_index(drop=True)
                    rs_df.to_pickle("./uploads/"+filename+".pkl")
                    return redirect(url_for("select"))
                else:
                    return render_template("input.html", error = "Error: File could not be parsed. Please check that this is a raw genotype file.")
            else:
                return render_template("input.html", error = "Error: Incorrect file type. Please upload a .txt file.")

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
    rs_df = pd.read_pickle("./uploads/"+filename+".pkl")
    rs_df_chr = rs_df.loc[rs_df['chr'] == chrom].values.tolist()
    cols = ['rsid', 'chromosome', 'position', 'genotype', 'substitution']
    df = pd.DataFrame(rs_df_chr, columns=cols)
    datatable = df.to_json(orient="table")
    jsontable = json.loads(datatable)
    return jsonify(jsontable)

@app.route('/nextsnp')
def next_snp():
    curr_snp = session["rsid"]
    filename = session["file"]
    rs_df = pd.read_pickle("./uploads/"+filename+".pkl")
    index = rs_df[rs_df['rsid']==curr_snp].index.values.astype(int)[0]
    if index == len(rs_df)-1:
        index = 0
    session["rsid"] = rs_df.loc[rs_df.index[index+1], 'rsid']
    session["genotype"] = rs_df.loc[rs_df.index[index+1], 'geno']
    return redirect(url_for("output"))

@app.route('/prevsnp')
def prev_snp():
    curr_snp = session["rsid"]
    filename = session["file"]
    rs_df = pd.read_pickle("./uploads/"+filename+".pkl")
    index = rs_df[rs_df['rsid']==curr_snp].index.values.astype(int)[0]
    if index == 0:
        index = len(rs_df)-1
    session["rsid"] = rs_df.loc[rs_df.index[index-1], 'rsid']
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
    print("here")
    # print(session["rsid"])
    if request.method == "POST":
        rsid, genotype, gene_name, chromosome, pdbs, first_pdb, first_aa, freq_kg, freq_hm, clinical, sorted_nuclist, sorted_aalist = session[
                "output"]
        if 'ddgcalc' in request.form:
            if first_aa != "N/A":
                ddgresults, chain = ddgcalcs(first_pdb, first_aa, gene_name, True)
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
