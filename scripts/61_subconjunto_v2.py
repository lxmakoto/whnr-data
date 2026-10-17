"""Busca subconjuntos COMPACTOS (todos los pares por encima del umbral) con
senal sinonima utilizable. El enlace simple encadenaba y saturaba el grupo."""
import numpy as np, json, random
def rd(p):
    d={};n=None
    for l in open(p):
        l=l.strip()
        if l.startswith(">"): n=l[1:]; d[n]=""
        elif l: d[n]+=l
    return d
cod=rd("results/selection/familia_codones.fasta")
ks=list(cod); A=np.array([list(c) for c in cod.values()]); n,L=A.shape
pos3=np.arange(2,L,3)
codes={c:i+1 for i,c in enumerate("ACGT")}
N=np.vectorize(lambda c: codes.get(c,0))(A).astype(np.int8); V=N>0
N3,V3=N[:,pos3],V[:,pos3]
ID3=np.zeros((n,n),dtype=np.float32)
for i in range(n):
    b=V3[i]&V3; s=(N3[i]==N3)&b
    c=b.sum(1); c[c==0]=1; ID3[i]=s.sum(1)/c

# referencia: el set de 8 genotipos que ya usamos
ids8=json.load(open("results/selection/subfamilia_ids.json")) if False else None
ocho=[l[1:].strip() for l in open("results/selection/cl33274978_codones_dedup99.fasta") if l.startswith(">")]
print(f"referencia (8 genotipos del clúster): identidad 3ª posición = ", end="")
# esos nombres son GCA_*, no estan en la familia; medir aparte
cod8=rd("results/selection/cl33274978_codones_dedup99.fasta")
A8=np.array([list(c) for c in cod8.values()]); p38=np.arange(2,A8.shape[1],3)
v=[]
for a in range(len(A8)):
    for b in range(a):
        m=(A8[a][p38]!="-")&(A8[b][p38]!="-")
        if m.sum(): v.append((A8[a][p38][m]==A8[b][p38][m]).mean())
print(f"{np.median(v):.3f}  <- esto es lo que buscamos igualar")

def greedy(seed,floor):
    """crece un grupo exigiendo que TODOS los pares superen floor en 3a posicion"""
    sel=[seed]
    cand=np.where(ID3[seed]>=floor)[0]
    cand=[c for c in cand if c!=seed]
    cand.sort(key=lambda c:-ID3[seed][c])
    for c in cand:
        if all(ID3[c][s]>=floor for s in sel): sel.append(c)
    return sel

print(f"\n{'umbral':>8}{'mayor grupo':>13}{'id 3ª mediana':>16}{'id global':>12}")
random.seed(0); seeds=random.sample(range(n),min(200,n))
best=None
for floor in (0.40,0.45,0.50,0.55):
    bg=[]
    for s in seeds:
        g=greedy(s,floor)
        if len(g)>len(bg): bg=g
    if len(bg)<5: print(f"{floor:>8.2f}{len(bg):>13}{'—':>16}{'—':>12}"); continue
    sub=ID3[np.ix_(bg,bg)]; iu=np.triu_indices(len(bg),1)
    med=float(np.median(sub[iu]))
    # identidad global
    IDg=np.zeros((len(bg),len(bg)))
    for a,x in enumerate(bg):
        for b,y in enumerate(bg):
            m=V[x]&V[y]; IDg[a,b]=((N[x]==N[y])&m).sum()/max(1,m.sum())
    print(f"{floor:>8.2f}{len(bg):>13}{med:>16.3f}{np.median(IDg[iu]):>12.3f}")
    if med>=0.45 and (best is None or len(bg)>len(best[1])): best=(floor,bg,med)

if best:
    floor,bg,med=best
    with open("results/selection/subfamilia_codones.fasta","w") as o:
        for i in bg: o.write(f">{ks[i]}\n{cod[ks[i]]}\n")
    print(f"\nELEGIDO: {len(bg)} secuencias, identidad 3ª pos. {med:.3f}")
    print("-> results/selection/subfamilia_codones.fasta")
else:
    print("\nNingún grupo compacto alcanza 0,45 en 3ª posición.")
