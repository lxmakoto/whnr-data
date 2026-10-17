"""Alineamiento de codones a escala de superfilo: retrotraduce el a3m de 742
homologos usando las CDS recuperadas de ENA/NCBI."""
import json, numpy as np, collections
B="TCAG"; AAS="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON={}; k=0
for a in B:
    for b in B:
        for c in B: CODON[a+b+c]=AAS[k]; k+=1
def tr(s): return "".join(CODON.get(s[i:i+3].upper(),"X") for i in range(0,len(s)-2,3))

cache=json.load(open("results/selection/cds_cache.json"))
pid_of=json.load(open("results/selection/pid_map.json"))["pid_of"]

rows=[];name=None
for line in open("results/msa/cl33274978.a3m"):
    line=line.rstrip("\n")
    if line.startswith(">"): rows.append([line[1:].split("\t")[0],[]])
    elif line: rows[-1][1].append(line)
a3m=[(n,"".join(s)) for n,s in rows]
print(f"filas del a3m: {len(a3m)}")

ok={}; why=collections.Counter()
for n,aln in a3m:
    if not n.startswith("UniRef100_"): continue
    acc=n.replace("UniRef100_","")
    pid=pid_of.get(acc)
    if not pid: why["sin protein_id"]+=1; continue
    nt=cache["cds"].get(pid)
    if not nt: why["sin CDS"]+=1; continue
    prot="".join(c for c in aln if c!="-").upper()
    p=tr(nt).rstrip("*")
    # El a3m trae la fila TRUNCADA a la region alineada: es subcadena de la CDS.
    if prot in p:
        off=p.index(prot)
    elif p==prot:
        off=0
    else:
        why["traducción discrepa"]+=1; continue
    if "*" in p[off:off+len(prot)]: why["parada interna"]+=1; continue
    ok[acc]=(aln,nt,off)
print(f"CDS verificadas: {len(ok)}")
for w,c in why.most_common(): print(f"   descartadas por {w}: {c}")

# retrotraduccion: columnas de coincidencia (MAYUSCULA/-); minuscula consume nt
cod={}
for acc,(aln,nt,off) in ok.items():
    i=off*3; out=[]
    for ch in aln:
        if ch==".": continue
        if ch=="-": out.append("---")
        elif ch.islower(): i+=3
        else:
            out.append(nt[i:i+3] if len(nt[i:i+3])==3 else "---"); i+=3
    cod[acc]="".join(out)
L=set(len(v) for v in cod.values())
assert len(L)==1, sorted(L)[:5]
Lc=L.pop(); print(f"alineamiento de codones: {len(cod)} x {Lc} nt ({Lc//3} codones)")

# comprobacion de coherencia codon/aminoacido
bad=0
for acc,(aln,_,_) in ok.items():
    pa="".join(c for c in aln if c=="-" or (c.isupper() and c!="."))
    for j in range(0,len(cod[acc]),3):
        c=cod[acc][j:j+3]; a=pa[j//3]
        if c=="---":
            if a!="-": bad+=1
        elif CODON.get(c.upper(),"X")!=a: bad+=1
print(f"codones incoherentes con su columna: {bad}")

with open("results/selection/familia_codones.fasta","w") as o:
    for acc,v in cod.items(): o.write(f">{acc}\n{v}\n")
print("-> results/selection/familia_codones.fasta")

# divergencia
ks=list(cod); A=np.array([list(cod[k]) for k in ks])
import random; random.seed(0)
idx=random.sample(range(len(ks)),min(120,len(ks)))
ids=[]
for a in range(len(idx)):
    for b in range(a):
        m=(A[idx[a]]!="-")&(A[idx[b]]!="-")
        if m.sum()>60: ids.append((A[idx[a]][m]==A[idx[b]][m]).mean())
ids=np.array(ids)
print(f"identidad nucleotídica por pares (muestra): mediana {np.median(ids):.3f} "
      f"p5 {np.percentile(ids,5):.3f} p95 {np.percentile(ids,95):.3f}")
print(f"pares >99% (casi clonales): {(ids>0.99).mean()*100:.1f}%")
