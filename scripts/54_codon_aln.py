"""Alineamiento de CODONES de los miembros del clúster cl33274978, para Datamonkey.

Datamonkey mide dN/dS: necesita nucleótido en marco, no proteína. Extraemos la
CDS de cada MAG por coordenadas, verificamos que traduce al producto almacenado,
alineamos las proteínas con un HMM de la familia y retrotraducimos a codones.
"""
import numpy as np, pandas as pd, pyhmmer, collections
from pyhmmer.easel import TextSequence, Alphabet

CL=33274978
nb=pd.read_parquet("data/interim/neighborhoods.parquet")
S=nb[(nb.cl_name==CL)&(nb.offset==0)].copy()
print(f"miembros del clúster: {len(S)}")

def load_fna(acc):
    d={}; name=None; buf=[]
    for line in open(f"data/genomes/{acc}_genomic.fna"):
        if line.startswith(">"):
            if name: d[name]="".join(buf)
            name=line[1:].split()[0]; buf=[]
        else: buf.append(line.strip())
    if name: d[name]="".join(buf)
    return d

def rc(s): return s.translate(str.maketrans("ACGTNacgtn","TGCANtgcan"))[::-1]
CODON={}
B="TCAG"
AAS="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
k=0
for a in B:
    for b in B:
        for c in B:
            CODON[a+b+c]=AAS[k]; k+=1
def tr(s):
    return "".join(CODON.get(s[i:i+3].upper(),"X") for i in range(0,len(s)-2,3))

cache={}; nt={}; aa={}
bad=[]
for _,r in S.iterrows():
    if r.acc not in cache:
        try: cache[r.acc]=load_fna(r.acc)
        except FileNotFoundError: bad.append((r.acc,"sin .fna")); continue
    ctg=cache[r.acc].get(r.ctg)
    if ctg is None: bad.append((r.acc,"contig ausente")); continue
    sub=ctg[int(r.beg)-1:int(r.end)]
    if r.strand=="-": sub=rc(sub)
    p=tr(sub).rstrip("*")
    ref=r.prot
    if p!=ref:
        # pyrodigal puede incluir el codón de parada o diferir en la Met inicial
        if p.lstrip("MVLI")==ref.lstrip("MVLI") and len(p)==len(ref): pass
        else:
            bad.append((r.acc,f"traducción no coincide ({len(p)} vs {len(ref)} aa)")); continue
    if "*" in p: bad.append((r.acc,"codón de parada interno")); continue
    nt[r.acc]=sub; aa[r.acc]=ref
print(f"CDS verificadas: {len(nt)}/{len(S)}")
for a,why in bad: print(f"   descartada {a}: {why}")

# --- alineamiento de proteína con HMM de la propia familia ---
abc=Alphabet.amino()
seqs=[TextSequence(name=k.encode(),sequence=v).digitize(abc) for k,v in aa.items()]
# El HMM se construye sobre el a3m de la familia (ya alineado, 191 columnas
# tras quitar las inserciones en minuscula). TextMSA exige longitud igual.
rows=[];nm=None
for line in open("results/msa/cl33274978.a3m"):
    line=line.rstrip("\n")
    if line.startswith(">"): rows.append([line[1:].split("\t")[0],[]])
    elif line: rows[-1][1].append(line)
fam=[(n,"".join(c for c in "".join(s) if not c.islower())) for n,s in rows]
Lq=len(fam[0][1]); fam=[(n,s) for n,s in fam if len(s)==Lq]
print(f"a3m de la familia para el HMM: {len(fam)} secuencias x {Lq}")
b=pyhmmer.plan7.Builder(abc); bg=pyhmmer.plan7.Background(abc)
msa=pyhmmer.easel.TextMSA(name=b"cl",
      sequences=[TextSequence(name=f"f{i}".encode(),sequence=s) for i,(n,s) in enumerate(fam)]
      ).digitize(abc)
hmm,_,_=b.build_msa(msa,bg)
aln=pyhmmer.hmmalign(hmm,seqs)
aln=aln.textize() if hasattr(aln,"textize") else aln
def dec(x): return x.decode() if isinstance(x,(bytes,bytearray)) else str(x)
raw={dec(s.name):str(s.sequence) for s in aln.sequences}

# hmmalign mezcla estados de coincidencia (MAYUSCULA / '-') con inserciones
# (minuscula / '.'), y pyhmmer no rellena las inserciones entre secuencias.
# Nos quedamos con los estados de coincidencia -> todas quedan a la longitud
# del HMM. Las inserciones consumen nucleotidos aunque no generen columna.
pal={}; cod={}
for k,ap in raw.items():
    s=nt[k]; i=0; pa=[]; cc=[]
    for ch in ap:
        if ch==".":                       # insercion sin residuo: nada
            continue
        if ch=="-":                       # delecion en estado de coincidencia
            pa.append("-"); cc.append("---")
        elif ch.islower():                # insercion con residuo: consume nt
            i+=3
        else:                             # estado de coincidencia
            pa.append(ch)
            cc.append(s[i:i+3] if len(s[i:i+3])==3 else "---"); i+=3
    pal[k]="".join(pa); cod[k]="".join(cc)
Lp=set(len(v) for v in pal.values()); Lc=set(len(v) for v in cod.values())
assert len(Lp)==1 and len(Lc)==1, (Lp,Lc)
print(f"alineamiento de proteína: {len(pal)} x {Lp.pop()}")
print(f"alineamiento de codones:  {len(cod)} x {Lc.pop()} nt")

with open("results/selection/cl33274978_codones.fasta","w") as o:
    for k,v in cod.items(): o.write(f">{k}\n{v}\n")
with open("results/selection/cl33274978_proteina.fasta","w") as o:
    for k,v in pal.items(): o.write(f">{k}\n{v}\n")

# comprobacion: cada codon traduce al aminoacido de su columna
errs=0
for k in cod:
    for j in range(0,len(cod[k]),3):
        c=cod[k][j:j+3]; a=pal[k][j//3]
        if c=="---":
            if a!="-": errs+=1
        elif CODON.get(c.upper(),"X")!=a.upper(): errs+=1
print(f"codones que no traducen a su columna: {errs}")

# --- divergencia: ¿hay senal suficiente para dN/dS? ---
ks=list(cod); A=np.array([list(cod[k]) for k in ks])
m=(A!="-").all(0)
sub=A[:,m]
ident=[]
for i in range(len(ks)):
    for j in range(i+1,len(ks)):
        ident.append((sub[i]==sub[j]).mean())
ident=np.array(ident)
print(f"\nidentidad nucleotídica por pares: mediana {np.median(ident):.3f}  "
      f"min {ident.min():.3f}  max {ident.max():.3f}")
print(f"pares con identidad > 0,99 (casi clonales): {(ident>0.99).mean()*100:.0f}%")
print("\n-> results/selection/cl33274978_codones.fasta")
print("-> results/selection/cl33274978_proteina.fasta")
