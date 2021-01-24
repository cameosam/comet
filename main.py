from flask import Flask, redirect, url_for, render_template, request, session, jsonify, flash
from werkzeug.utils import secure_filename
from Bio import Entrez
import itertools
import os
import json
import pandas as pd
# from flask_session import Session
# from datetime import timedelta

from helperfun import *
from ddgcalc import ddgcalcs

app = Flask(__name__)
app.secret_key = "secret"
app.config['UPLOAD_FOLDER'] = './uploads/'
# app.config['SESSION_PERMANENT'] = True
# app.config['SESSION_TYPE'] = 'filesystem'
# app.config['PERMANENT_SESSION_LIFETIME'] = timedelta(hours=1)
# app.config['SESSION_FILE_THRESHOLD'] = 100 
# app.config['SECRET_KEY'] = app.secret_key
# sess = Session()
# sess.init_app(app)

@app.route("/")
@app.route("/home")
def index():
    return render_template("index.html")

@app.route("/<string:page_name>")
def html_page(page_name):
    return render_template(page_name)

@app.route("/endsession")
def endsession():
    if session["file"]:
        os.remove('uploads/'+session["file"])
    for key in ['chrom', 'file', 'genotype', 'output', 'rsid']:
        session.pop(key)
    flash('Thanks for using COMET. All session data has been deleted')
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

@app.route("/getchoice", methods=["POST"])
def testid():
    if request.method == "POST":
        session["rsid"] = request.args.get('rsid', '')
        session["genotype"] = request.args.get('genotype', '')
        return redirect(url_for("output"))

@app.route("/output", methods=["POST", "GET"])
def output():
    if request.method == "POST":
        if 'pdbselect' in request.form:
            pdbselect = request.form['pdbselect']
            rsid, genotype, gene_name, chromosome, pdbs, first_pdb, first_aa, ddgresults, freq_kg, freq_hm, clinical, subs, chain = session[
                "output"]
            if first_aa != "N/A":
                ddgresults, chain = ddgcalcs(pdbselect, first_aa, gene_name)
            return render_template("output.html", snp=rsid, genotype = genotype, gene=gene_name, chr=chromosome, pdb=pdbs, pdbselect=pdbselect, aa1=first_aa, ddgresults=ddgresults, freq1000g=freq_kg, freqhapmap=freq_hm, clin=clinical, subs=subs, chain=chain, zip=zip, len=len)

    if "rsid" in session:
        rsid = session["rsid"]
        genotype = session["genotype"]
        gene_name, chromosome, freq_kg, freq_hm, clinical, subs, first_aa = getsnpinfo(rsid) 
        pdbs = findpdb(gene_name)
        # pdbs = pdbfilter(pdbs_all, gene_name)

        # calculate ddg
        if pdbs != "N/A":
            first_pdb = pdbs[0] if len(pdbs) > 1 else pdbs
            if first_aa != "N/A":
                ddgresults, chain = ddgcalcs(first_pdb, first_aa, gene_name)
            else:
                chain = "*"
                ddgresults = [["N/A" for i in range(2)] for j in range(4)]
        else:
            first_pdb = "N/A"

        session["output"] = [rsid, genotype, gene_name, chromosome, pdbs, first_pdb, first_aa, ddgresults, freq_kg, freq_hm, clinical, subs, chain]
 
        print(session["output"])
        return render_template("output.html", snp=rsid, genotype= genotype,gene=gene_name,
                               chr=chromosome, pdb=pdbs, pdbselect=first_pdb,
                               aa1=first_aa, ddgresults=ddgresults, freq1000g=freq_kg,
                               freqhapmap=freq_hm, clin=clinical,
                               subs=subs, chain=chain,
                               zip=zip, len=len)
    else:
        return redirect(url_for("select"))


if __name__ == "__main__":
    app.run(debug=True)
