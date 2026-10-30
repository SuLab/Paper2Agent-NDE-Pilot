"""For each GEO series: sample sheet summary, supplementary file names, NCBI raw-count availability."""
import collections, gzip, io, json, re, sys, urllib.request, urllib.error
out_dir, ids = sys.argv[1], json.load(open(sys.argv[2]))
def get(url, head=False):
    req = urllib.request.Request(url, method="HEAD" if head else "GET")
    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            return r.status, (b"" if head else r.read()), r.headers.get("Content-Type", "")
    except urllib.error.HTTPError as e:
        return e.code, b"", ""
for gse in ids:
    pre = gse[:-3] + "nnn"
    base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{pre}/{gse}"
    st, raw, _ = get(f"{base}/soft/{gse}_family.soft.gz")
    if st != 200:
        print(f"## {gse}: SOFT HTTP {st}"); continue
    open(f"{out_dir}/{gse}_family.soft.gz", "wb").write(raw)
    lines = gzip.decompress(raw).decode("utf-8", "replace").splitlines()
    samples, cur = [], None
    for l in lines:
        if l.startswith("^SAMPLE"):
            cur = {"chars": []}; samples.append(cur)
        elif cur is not None and "=" in l:
            k, v = l.split("=", 1); k, v = k.strip(), v.strip()
            if k == "!Sample_title": cur["title"] = v
            elif k == "!Sample_characteristics_ch1": cur["chars"].append(v)
            elif k == "!Sample_library_source": cur["src"] = v
            elif k == "!Sample_library_strategy": cur["strat"] = v
    title = next((l.split("=", 1)[1].strip() for l in lines if l.startswith("!Series_title")), "")
    _, listing, _ = get(f"{base}/suppl/")
    files = re.findall(r'href="([^"?/][^"]*)"', listing.decode("utf-8", "replace"))
    st2, fl, _ = get(f"{base}/suppl/filelist.txt")
    tar_files = [x.split("\t")[1] for x in fl.decode("utf-8", "replace").splitlines()[1:] if "\t" in x] if st2 == 200 else []
    stc, _, ctype = get(f"https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc={gse}&format=file&file={gse}_raw_counts_GRCh38.p13_NCBI.tsv.gz", head=True)
    print(f"## {gse} | n={len(samples)} | {title[:110]}")
    print("   source:", dict(collections.Counter(s.get("src") for s in samples)), "| strategy:", dict(collections.Counter(s.get("strat") for s in samples)))
    print(f"   NCBI raw counts: HTTP {stc} {ctype[:30]} | suppl: {[f for f in files if f != 'filelist.txt'][:6]}")
    if tar_files:
        exts = collections.Counter(re.sub(r"^.*?(\.[A-Za-z0-9]+(\.gz|\.bz2|\.zip)?)$", r"\1", f).lower() for f in tar_files)
        print(f"   RAW.tar: {len(tar_files)} files, extensions {dict(exts.most_common(6))}, e.g. {tar_files[:3]}")
    keys = collections.defaultdict(collections.Counter)
    for s in samples:
        for c in s["chars"]:
            k, _, v = c.partition(":"); keys[k.strip()][v.strip()] += 1
    for k, cnt in keys.items():
        if 1 < len(cnt) <= 12: print(f"   {k}: {dict(cnt)}")
        elif len(cnt) > 12: print(f"   {k}: {len(cnt)} distinct values")
    print("   titles:", "; ".join(s.get("title", "") for s in samples[:8])[:300])
