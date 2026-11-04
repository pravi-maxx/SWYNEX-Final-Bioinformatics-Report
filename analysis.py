"""SWYNEX Task 4 - Final Bioinformatics Report pipeline.
Reads sequences.fasta (headers: >id|group) if present; otherwise simulates
coding sequences (SYNTHETIC demo) so the pipeline runs end to end.
Writes charts + results.json used by REPORT.md."""
import os, json, itertools
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

B = "TCAG"
CODONS = ["".join(c) for c in itertools.product(B, repeat=3)]
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
TABLE = dict(zip(CODONS, AA))
SENSE = [c for c in CODONS if TABLE[c] != "*"]
rng = np.random.default_rng(7)

def read_fasta(path):
    recs, name, seq = [], None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name: recs.append((name, "".join(seq)))
            name, seq = line[1:], []
        elif line: seq.append(line.upper())
    if name: recs.append((name, "".join(seq)))
    return recs

def simulate():
    recs = []
    for grp, w in (("HighGC", 2.5), ("LowGC", 0.4)):
        weights = np.array([w if c[2] in "GC" else 1.0 for c in SENSE]); weights /= weights.sum()
        for i in range(40):
            n = int(rng.integers(100, 300))
            body = rng.choice(SENSE, n, p=weights)
            recs.append((f"{grp}_{i}|{grp}", "ATG" + "".join(body) + "TAA"))
    return recs

synthetic = not os.path.exists("sequences.fasta")
recs = simulate() if synthetic else read_fasta("sequences.fasta")

rows = []
for name, s in recs:
    grp = name.split("|")[-1] if "|" in name else "All"
    s = s[: len(s) // 3 * 3]
    cods = [s[i:i+3] for i in range(0, len(s), 3)]
    gc = (s.count("G") + s.count("C")) / len(s)
    gc3 = np.mean([c[2] in "GC" for c in cods])
    prot = "".join(TABLE.get(c, "X") for c in cods).rstrip("*")
    rows.append(dict(id=name.split("|")[0], group=grp, length_nt=len(s), gc=gc, gc3=gc3, prot=prot))
df = pd.DataFrame(rows)
df[["id","group","length_nt","gc","gc3"]].to_csv("sequence_stats.csv", index=False)

groups = sorted(df.group.unique())
# Chart 1: GC by group
plt.figure(figsize=(6,4))
plt.boxplot([df[df.group==g].gc for g in groups], tick_labels=groups)
plt.ylabel("GC content"); plt.title("GC content by group")
plt.tight_layout(); plt.savefig("fig1_gc_by_group.png", dpi=150); plt.close()
# Chart 2: GC vs GC3
plt.figure(figsize=(6,4))
for g in groups:
    d = df[df.group==g]; plt.scatter(d.gc3, d.gc, alpha=.7, label=g)
r = df.gc3.corr(df.gc)
plt.xlabel("GC3 (GC at 3rd codon position)"); plt.ylabel("Overall GC"); plt.title(f"GC vs GC3 (r={r:.2f})"); plt.legend()
plt.tight_layout(); plt.savefig("fig2_gc_vs_gc3.png", dpi=150); plt.close()
# Chart 3: amino-acid composition
aas = sorted(set(AA) - {"*"})
comp = {g: np.array([sum(p.count(a) for p in df[df.group==g].prot) for a in aas], float) for g in groups}
plt.figure(figsize=(8,4)); w = 0.8/len(groups)
for k, g in enumerate(groups):
    plt.bar(np.arange(len(aas)) + k*w, comp[g]/comp[g].sum()*100, w, label=g)
plt.xticks(np.arange(len(aas)) + w/2, aas); plt.ylabel("% of residues"); plt.title("Amino-acid composition"); plt.legend()
plt.tight_layout(); plt.savefig("fig3_aa_composition.png", dpi=150); plt.close()
# Chart 4: length vs GC
plt.figure(figsize=(6,4))
for g in groups:
    d = df[df.group==g]; plt.scatter(d.length_nt, d.gc, alpha=.7, label=g)
rl = df.length_nt.corr(df.gc)
plt.xlabel("Sequence length (nt)"); plt.ylabel("GC content"); plt.title(f"Length vs GC (r={rl:.2f})"); plt.legend()
plt.tight_layout(); plt.savefig("fig4_length_vs_gc.png", dpi=150); plt.close()

res = dict(synthetic=synthetic, n=len(df), groups=groups, r_gc_gc3=round(r,3), r_len_gc=round(rl,3))
for g in groups:
    d = df[df.group==g]
    res[g] = dict(n=len(d), gc_mean=round(d.gc.mean(),3), gc_sd=round(d.gc.std(),3), gc3_mean=round(d.gc3.mean(),3), len_mean=round(d.length_nt.mean(),1))
if len(groups) == 2:
    a, b = (df[df.group==g].gc.values for g in groups)
    obs = a.mean() - b.mean(); pool = np.concatenate([a, b]); cnt = 0
    for _ in range(10000):
        rng.shuffle(pool); cnt += abs(pool[:len(a)].mean() - pool[len(a):].mean()) >= abs(obs)
    res["diff_gc"] = round(obs,3); res["perm_p"] = (cnt+1)/10001
json.dump(res, open("results.json","w"), indent=1)
print(json.dumps(res, indent=1))
