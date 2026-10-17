"""Extras para Datamonkey: mapa columna->residuo, version deduplicada y arbol."""
import numpy as np, pandas as pd, itertools
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

def rd(p):
    d={};n=None
    for l in open(p):
        l=l.strip()
        if l.startswith(">"): n=l[1:]; d[n]=""
        elif l: d[n]+=l
    return d
cod=rd("results/selection/cl33274978_codones.fasta")
pal=rd("results/selection/cl33274978_proteina.fasta")
ks=list(cod); L=len(pal[ks[0]])
QUERY="GCA_001777825.1"      # representante de AgnostosDB: la numeracion del paper (E141, D121...)
assert QUERY in pal, "falta el genoma de referencia"

# --- mapa columna del alineamiento -> numeracion de residuo de la referencia ---
ref=pal[QUERY]; pos=0; rows=[]
for i,ch in enumerate(ref,1):
    if ch!="-": pos+=1
    rows.append(dict(columna_aln=i, codon_nt=f"{(i-1)*3+1}-{i*3}",
                     residuo_referencia=pos if ch!="-" else None, aa=ch))
m=pd.DataFrame(rows)
m.to_csv("results/selection/mapa_columnas.tsv",sep="\t",index=False)
INT=[119,121,141,147,157,161,166,168,110,126]
sel=m[m.residuo_referencia.isin(INT)]
print("residuos de interés y su columna en el alineamiento:")
for _,r in sel.iterrows():
    print(f"   {r.aa}{int(r.residuo_referencia):<4} -> columna {int(r.columna_aln)}  (nt {r.codon_nt})")

# --- version deduplicada al 99% ---
A=np.array([list(cod[k]) for k in ks]); keep=[]; drop={}
for i,k in enumerate(ks):
    dup=None
    for j in keep:
        m_=(A[i]!="-")&(A[j]!="-")
        if m_.sum() and (A[i][m_]==A[j][m_]).mean()>0.99: dup=ks[j]; break
    if dup is None: keep.append(i)
    else: drop[k]=dup
with open("results/selection/cl33274978_codones_dedup99.fasta","w") as o:
    for j in keep: o.write(f">{ks[j]}\n{cod[ks[j]]}\n")
print(f"\ndeduplicado al 99%: {len(keep)}/{len(ks)} secuencias")
for a,b in drop.items(): print(f"   {a} ≈ {b}")

# --- arbol NJ sobre el alineamiento de codones (ambas versiones) ---
for tag,idx in (("",list(range(len(ks)))),("_dedup99",keep)):
    names=[ks[i] for i in idx]
    D=[[0.0]*(r+1) for r in range(len(idx))]
    for a in range(len(idx)):
        for b in range(a):
            m_=(A[idx[a]]!="-")&(A[idx[b]]!="-")
            D[a][b]=float(1-(A[idx[a]][m_]==A[idx[b]][m_]).mean()) if m_.sum() else 1.0
    t=DistanceTreeConstructor().nj(DistanceMatrix(names,D))
    t.root_at_midpoint()
    for c in t.get_nonterminals(): c.name=None
    Phylo.write(t,f"results/selection/cl33274978{tag}.nwk","newick")
    print(f"-> results/selection/cl33274978{tag}.nwk")
print("-> results/selection/mapa_columnas.tsv")
print("-> results/selection/cl33274978_codones_dedup99.fasta")
