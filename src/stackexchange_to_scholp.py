"""stackexchange_to_scholp.py — convert a Stack Exchange data-dump Posts.xml into ScHoLP-format
simplex streams (data/ScHoLP-Data/<name>/<name>-{nverts,simplices,times}.txt.gz) for src/closure.py.

Two streams per site, mirroring Benson's tags-* and threads-*:
  tags-<site>     : one question = its set of tags            (nodes = tags)
  threads-<site>  : one question = {asker} ∪ {answerers}      (nodes = users; time = question time)

usage: python3 src/stackexchange_to_scholp.py <dump>/matheducators/Posts.xml matheducators [--out data/ScHoLP-Data]

Copied from the CNA 2026 project (Benson/src/stackexchange_to_scholp.py) on 2026-10-01.
Change: CreationDate is read as UTC, so the times files do not depend on the machine's time zone
(the simplices and their order are unchanged).
"""
import gzip, os, sys, time
import xml.etree.ElementTree as ET
from collections import defaultdict, Counter
from datetime import datetime, timezone

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "ScHoLP-Data")


def ptime(s):
    return int(datetime.strptime(s[:19], "%Y-%m-%dT%H:%M:%S").replace(tzinfo=timezone.utc).timestamp())


def split_tags(s):
    # old dumps: "<a><b>", new dumps (2024+): "|a|b|"
    if s.startswith("<"):
        return [t for t in s.strip("<>").split("><") if t]
    return [t for t in s.split("|") if t]


def write_scholp(name, simplices):
    simplices.sort(key=lambda x: x[0])
    ids = {}
    d = os.path.join(OUT, name); os.makedirs(d, exist_ok=True)
    with gzip.open(f"{d}/{name}-nverts.txt.gz", "wt") as fn, \
         gzip.open(f"{d}/{name}-simplices.txt.gz", "wt") as fs, \
         gzip.open(f"{d}/{name}-times.txt.gz", "wt") as ft:
        for t, s in simplices:
            fn.write(f"{len(s)}\n"); ft.write(f"{t}\n")
            for v in s:
                fs.write(f"{ids.setdefault(v, len(ids) + 1)}\n")
    with gzip.open(f"{d}/{name}-nodes.txt.gz", "wt", encoding="utf-8") as fo:
        for v, i in sorted(ids.items(), key=lambda x: x[1]):
            fo.write(f"{i}\t{v}\n")
    sz = [len(s) for _, s in simplices]; c = Counter(min(k, 6) for k in sz)
    print(f"[{name}] simplices={len(simplices)} nodes={len(ids)} mean={sum(sz)/len(sz):.2f} max={max(sz)} "
          f"size(2..5,6+)={[c[k] for k in range(2,7)]} "
          f"span={datetime.fromtimestamp(simplices[0][0], timezone.utc).date()}..{datetime.fromtimestamp(simplices[-1][0], timezone.utc).date()}", flush=True)


def main(path, site):
    t0 = time.time()
    q_time, q_tags, q_owner, answerers = {}, {}, {}, defaultdict(set)
    n = 0
    for _, el in ET.iterparse(path, events=("end",)):
        if el.tag != "row":
            continue
        n += 1
        pt = el.get("PostTypeId")
        if pt == "1":
            qid = el.get("Id")
            q_time[qid] = ptime(el.get("CreationDate"))
            q_tags[qid] = split_tags(el.get("Tags", ""))
            if el.get("OwnerUserId"):
                q_owner[qid] = el.get("OwnerUserId")
        elif pt == "2" and el.get("OwnerUserId"):
            answerers[el.get("ParentId")].add(el.get("OwnerUserId"))
        el.clear()
    print(f"posts={n} questions={len(q_time)} ({time.time()-t0:.0f}s)")
    tags = [(q_time[q], sorted(set(tg))) for q, tg in q_tags.items() if len(set(tg)) >= 2]
    write_scholp(f"tags-{site}", tags)
    th = []
    for q, t in q_time.items():
        us = set(answerers.get(q, ())) | ({q_owner[q]} if q in q_owner else set())
        if len(us) >= 2:
            th.append((t, sorted(us)))
    write_scholp(f"threads-{site}", th)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("posts_xml"); ap.add_argument("site")
    ap.add_argument("--out", default=OUT, help="output root (default: <repo>/data/ScHoLP-Data)")
    a = ap.parse_args()
    OUT = a.out
    main(a.posts_xml, a.site)
