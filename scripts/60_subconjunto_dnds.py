"""Busca el subconjunto mas grande de las 716 CDS que conserve senal sinonima,
es decir, identidad en 3a posicion del codon por encima del umbral de saturacion."""
import numpy as np, collections, json

def rd(p):
    d={};n=None
    for l in open(p):
        l=l.strip()
        if l.startswith(">"): n=l[1:]; d[n]=""
        elif l: d[n]+=l
    return d
cod=rd("results/selection/familia_codones.fasta")
ks=list(cod); A=np.array([list(c) for c in cod.values()])
n,L=A.shape
pos3=np.arange(2,L,3)

# identidad aminoacidica por pares (sobre el alineamiento de codones traducido)
AA=np.array([[c for c in "".join(r)] for r in A])
print(f"secuencias: {n}, columnas: {L}")

# matriz de identidad nucleotidica completa (vectorizada)
num=np.where(A=="-",0,np.char.encode(A).view(np.uint8).reshape(A.shape) if False else 0)
codes={c:i+1 for i,c in enumerate("ACGT")}
N=np.vectorize(lambda c: codes.get(c,0))(A).astype(np.int8)
V=N>0
print("calculando matriz de identidad...",flush=True)
ID=np.zeros((n,n),dtype=np.float32); ID3=np.zeros((n,n),dtype=np.float32)
N3=N[:,pos3]; V3=V[:,pos3]
for i in range(n):
    both=V[i]&V; same=(N[i]==N)&both
    c=both.sum(1); c[c==0]=1; ID[i]=same.sum(1)/c
    b3=V3[i]&V3; s3=(N3[i]==N3)&b3
    c3=b3.sum(1); c3[c3==0]=1; ID3[i]=s3.sum(1)/c3

# agrupamiento aglomerativo por enlace simple a distintos umbrales de identidad
def grupos(thr):
    """componentes conexas con identidad >= thr"""
    seen=np.zeros(n,bool); out=[]
    adj=ID>=thr
    for i in range(n):
        if seen[i]: continue
        stack=[i]; comp=[]
        seen[i]=True
        while stack:
            x=stack.pop(); comp.append(x)
            for y in np.where(adj[x]&~seen)[0]:
                seen[y]=True; stack.append(y)
        out.append(comp)
    return sorted(out,key=len,reverse=True)

print(f"\n{'umbral':>8}{'nº grupos':>11}{'mayor':>8}{'id 3ª pos. del mayor':>24}{'veredicto':>14}")
best=None
for thr in (0.50,0.55,0.60,0.65,0.70,0.75):
    gs=grupos(thr); big=gs[0]
    if len(big)<6: 
        print(f"{thr:>8.2f}{len(gs):>11}{len(big):>8}{'—':>24}{'muy pequeño':>14}")
        continue
    sub=ID3[np.ix_(big,big)]
    med=float(np.median(sub[np.triu_indices(len(big),1)]))
    ver="USABLE" if med>=0.45 else ("límite" if med>=0.40 else "saturado")
    print(f"{thr:>8.2f}{len(gs):>11}{len(big):>8}{med:>24.3f}{ver:>14}")
    if med>=0.45 and (best is None or len(big)>len(best[1])): best=(thr,big,med)

if best:
    thr,big,med=best
    names=[ks[i] for i in big]
    with open("results/selection/subfamilia_codones.fasta","w") as o:
        for i in big: o.write(f">{ks[i]}\n{cod[ks[i]]}\n")
    print(f"\nSUBCONJUNTO ELEGIDO: {len(big)} secuencias (umbral {thr:.2f}), "
          f"identidad 3ª pos. mediana {med:.3f}")
    sub=ID[np.ix_(big,big)]; iu=np.triu_indices(len(big),1)
    print(f"   identidad nucleotídica global: mediana {np.median(sub[iu]):.3f}")
    print(f"   pares casi clonales (>0,99): {100*(sub[iu]>0.99).mean():.1f}%")
    json.dump(names,open("results/selection/subfamilia_ids.json","w"))
    print("-> results/selection/subfamilia_codones.fasta")
else:
    print("\nNingún subconjunto supera 0,45 de identidad en 3ª posición.")
