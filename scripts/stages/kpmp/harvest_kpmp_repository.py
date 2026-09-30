"""Harvest the KPMP Atlas Repository file index (public Elastic App Search engine used by
https://atlas.kpmp.org/repository/). The search key is the public read-only key embedded in the
repository front-end bundle. Queries are partitioned by experimental_strategy because App Search
returns at most 10,000 hits per query."""
import json, sys, time, urllib.request, csv, os
URL = "https://atlas.kpmp.org/spatial-viewer/search/api/as/v1/engines/atlas-repository/search.json"
KEY = "search-vwz67uj2sf8h83h4y8i8j6g3"
OUT = sys.argv[1] if len(sys.argv) > 1 else "."

def post(body):
    req = urllib.request.Request(URL, data=json.dumps(body).encode(), headers={
        "Authorization": "Bearer " + KEY, "Content-Type": "application/json"})
    for i in range(5):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)
        except Exception as e:
            time.sleep(3 * (i + 1)); err = e
    raise err

def flat(v):
    v = v.get("raw") if isinstance(v, dict) else v
    if isinstance(v, list): return "|".join(map(str, v))
    return "" if v is None else str(v)

strats = [x["value"] for x in post({"query": "", "page": {"size": 1},
          "facets": {"experimental_strategy": {"type": "value", "size": 250}}}
          )["facets"]["experimental_strategy"][0]["data"]]
rows, seen = [], set()
for s in strats:
    page = 1
    while True:
        d = post({"query": "", "page": {"size": 1000, "current": page},
                  "filters": {"all": [{"experimental_strategy": s}]},
                  "sort": [{"file_id": "asc"}]})
        for r in d["results"]:
            rec = {k: flat(v) for k, v in r.items() if k != "_meta"}
            if rec["id"] not in seen:
                seen.add(rec["id"]); rows.append(rec)
        tp = d["meta"]["page"]["total_pages"]
        print(s, page, tp, d["meta"]["page"]["total_results"], file=sys.stderr)
        if page >= tp: break
        page += 1
# documents with no experimental_strategy (e.g. OpenAccessClinicalData.csv)
d = post({"query": "", "page": {"size": 1000},
          "filters": {"none": [{"experimental_strategy": strats}]}})
for r in d["results"]:
    rec = {k: flat(v) for k, v in r.items() if k != "_meta"}
    if rec["id"] not in seen:
        seen.add(rec["id"]); rows.append(rec)
print("no-strategy docs", d["meta"]["page"]["total_results"], file=sys.stderr)
tot = post({"query": "", "page": {"size": 1}})["meta"]["page"]["total_results"]
cols = sorted({k for r in rows for k in r})
first = ["redcap_id", "access", "experimental_strategy", "workflow_type", "data_type",
         "data_format", "file_name", "file_size", "package_id", "file_id"]
cols = first + [c for c in cols if c not in first]
for r in rows:
    r["download_url"] = ("https://atlas.kpmp.org/api/v1/file/download/%s/%s" % (r["package_id"], r["file_name"])
                         if r.get("access") == "open" else "")
cols.append("download_url")
with open(os.path.join(OUT, "kpmp_repository_all_files.tsv"), "w", newline="") as f:
    w = csv.DictWriter(f, cols, delimiter="\t", extrasaction="ignore"); w.writeheader(); w.writerows(rows)
print("harvested", len(rows), "unique docs; engine total_results (capped at 10000):", tot, file=sys.stderr)
