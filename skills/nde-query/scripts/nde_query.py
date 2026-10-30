"""Query the NIAID Data Ecosystem (NDE) Discovery API and save the raw JSON.

Usage examples:
  python nde_query.py -q '@type:"Dataset" AND infectiousAgent.name:"Zika virus"' --size 0 \
      --facets includedInDataCatalog.name,measurementTechnique.name.raw
  python nde_query.py -q '@type:"Dataset" AND ...' --fields identifier,name,sample.itemListElement.identifier \
      --size 1000 --out hits.json --table
  python nde_query.py -q '...' --all --fields identifier --out all.json   # scroll past the 1000-hit cap

Prints the total, any facet counts, and (with --table) one line per hit:
identifier, catalog, number of samples listed, name. Standard library only.
"""
import argparse
import json
import sys
import urllib.parse
import urllib.request

API = "https://api.data.niaid.nih.gov/v1/query"


def get(params: dict, base: str) -> dict:
    url = base + "?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=180) as r:
        return json.load(r)


def as_list(x):
    if x is None:
        return []
    return x if isinstance(x, list) else [x]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-q", "--query", required=True, help="Elasticsearch query_string, e.g. '@type:\"Dataset\" AND ...'")
    ap.add_argument("--size", type=int, default=10, help="hits to return (max 1000; 0 = counts/facets only)")
    ap.add_argument("--from", dest="offset", type=int, default=0, help="offset for paging")
    ap.add_argument("--fields", help="comma-separated fields to return (default: whole records)")
    ap.add_argument("--facets", help="comma-separated fields to aggregate on")
    ap.add_argument("--facet-size", type=int, default=20)
    ap.add_argument("--extra-filter", help="filter applied to hits and facets, e.g. 'species.name:(\"homo sapiens\")'")
    ap.add_argument("--all", action="store_true", help="fetch every hit with fetch_all/scroll (ignores --size/--from)")
    ap.add_argument("--out", help="write the raw JSON here")
    ap.add_argument("--table", action="store_true", help="print one line per hit")
    ap.add_argument("--base", default=API, help="API endpoint (e.g. the api-staging host)")
    a = ap.parse_args()

    params = {"q": a.query}
    if a.fields:
        params["fields"] = a.fields
    if a.extra_filter:
        params["extra_filter"] = a.extra_filter
    if a.all:
        params["fetch_all"] = "true"
        data = get(params, a.base)
        hits = list(data.get("hits", []))
        while data.get("_scroll_id") and len(hits) < data.get("total", 0):
            data = get({"scroll_id": data["_scroll_id"]}, a.base)
            if not data.get("hits"):
                break
            hits.extend(data["hits"])
        data = {"total": len(hits), "hits": hits}
    else:
        params.update({"size": a.size, "from": a.offset})
        if a.facets:
            params.update({"facets": a.facets, "facet_size": a.facet_size})
        data = get(params, a.base)

    if data.get("success") is False:
        sys.exit(f"API error: {json.dumps(data)[:500]}")
    if a.out:
        with open(a.out, "w") as fh:
            json.dump(data, fh)
    print("total:", data.get("total"))
    for name, facet in (data.get("facets") or {}).items():
        print(f"facet {name}:")
        for t in facet.get("terms", []):
            print(f"  {t['count']:>10}  {t['term']}")
    if a.table:
        for h in data.get("hits", []):
            cat = as_list(h.get("includedInDataCatalog"))
            n = len(as_list((h.get("sample") or {}).get("itemListElement")))
            print(f"{h.get('identifier', ''):<22} {(cat[0].get('name', '') if cat else '')[:14]:<14} "
                  f"n={n:<4} {str(h.get('name', ''))[:100]}")


if __name__ == "__main__":
    main()
