"""Map NDE SRA/BioProject records to linked GEO series via NCBI E-utilities."""
import json, sys, time, urllib.parse, urllib.request
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
def call(tool, **p):
    time.sleep(0.4)
    with urllib.request.urlopen(E + tool + "?" + urllib.parse.urlencode({**p, "retmode": "json"}), timeout=60) as r:
        return json.load(r)
recs = json.load(open(sys.argv[1])); out = {}
for h in recs:
    acc = h.get("identifier") or h.get("_id")
    if acc.startswith("GSE"):
        out[acc] = [acc]; continue
    if acc.startswith("ncbi_sra_"):
        srp = acc.split("_")[-1].upper()
        ids = call("esearch.fcgi", db="sra", term=srp, retmax=1)["esearchresult"]["idlist"]
        bp = []
        if ids:
            ls = call("elink.fcgi", dbfrom="sra", db="bioproject", id=ids[0])["linksets"][0].get("linksetdbs", [])
            bp = [i for l in ls for i in l["links"]]
    else:
        bp = call("esearch.fcgi", db="bioproject", term=f"{acc}[Project Accession]")["esearchresult"]["idlist"]
    gse = []
    for b in bp[:1]:
        ls = call("elink.fcgi", dbfrom="bioproject", db="gds", id=b)["linksets"][0].get("linksetdbs", [])
        uids = [i for l in ls for i in l["links"]]
        if uids:
            summ = call("esummary.fcgi", db="gds", id=",".join(uids[:20]))["result"]
            gse += sorted({summ[u]["accession"] for u in summ.get("uids", []) if summ[u].get("accession", "").startswith("GSE")})
    out[acc] = gse
    print(acc, "->", gse or "no linked GEO series")
json.dump(out, open(sys.argv[2], "w"))
