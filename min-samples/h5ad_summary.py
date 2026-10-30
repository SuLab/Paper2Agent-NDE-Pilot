import sys, numpy as np, anndata as ad
for p in sys.argv[1:]:
    a = ad.read_h5ad(p, backed="r")
    x = a.X[:500]
    x = x.toarray() if hasattr(x, "toarray") else np.asarray(x)
    print("==", p.split("/")[-2], a.shape, "| integer counts:", bool(np.all(x % 1 == 0)), "| max:", float(x.max()))
    print("   var head:", list(a.var_names[:5]))
    for c in a.obs.columns:
        v = a.obs[c]
        if v.nunique() <= 40:
            print(f"   obs {c} ({v.nunique()}):", dict(v.value_counts().head(8)))
        else:
            print(f"   obs {c}: {v.nunique()} unique")
