"""Control emparejado: conservacion de carboxilatos i/i+2 en una familia pol-beta
ACTIVA frente a la nuestra. Mismo servidor MSA, mismo codigo, mismo ponderado."""
import numpy as np, collections

def load_a3m(p):
    rows=[];
    for line in open(p):
        line=line.rstrip("\n")
        if line.startswith(">"): rows.append([])
        elif line: rows[-1].append(line)
    seqs=["".join(c for c in "".join(r) if not c.islower()) for r in rows]
    L=len(seqs[0]); seqs=[s for s in seqs if len(s)==L]
    return np.array([list(s) for s in seqs])

def weights(M):
    """1/(nº de secuencias a >=80% de identidad). Identico para ambas familias."""
    A=np.array([[ord(c) for c in r] for r in M],dtype=np.int16); gap=ord("-")
    v=(A!=gap); n=len(A); w=np.empty(n)
    for i in range(n):
        both=v[i]&v; same=(A[i]==A)&both
        cnt=both.sum(1); cnt[cnt==0]=1
        w[i]=1.0/ (( same.sum(1)/cnt >=0.80).sum())
    return w

def consw(M,w,i):
    col=M[:,i]; m=col!="-"
    if not m.any(): return "-",0.0,0.0,0.0
    ww=w[m]; tot=ww.sum()
    d=collections.defaultdict(float)
    for a,x in zip(col[m],ww): d[a]+=x
    aa=max(d,key=d.get)
    fDE=(d.get("D",0)+d.get("E",0))/tot
    return aa, d[aa]/tot, fDE, m.sum()/len(col)

CASES=[("pol-β ACTIVA (KNTasa P05057)","results/control_ntase/KNTase.a3m",50,52,"Asp50","Glu52"),
       ("nuestra familia (cl33274978)","results/msa/cl33274978.a3m",119,121,"Glu119","Asp121")]
out={}
for lab,path,i1,i2,n1,n2 in CASES:
    M=load_a3m(path); w=weights(M)
    neff=w.sum()
    print(f"\n=== {lab} ===")
    print(f"secuencias {len(M)}  ·  columnas {M.shape[1]}  ·  Neff {neff:.0f}")
    print(f"{'posición':>10}{'consulta':>10}{'aa cons.':>10}{'frac':>8}{'frac D/E':>10}{'cobert.':>9}")
    for p in range(i1-3,i2+4):
        if not (1<=p<=M.shape[1]): continue
        aa,f,fDE,cov=consw(M,w,p-1)
        mark="  <-- i" if p==i1 else ("  <-- i+2" if p==i2 else "")
        print(f"{p:>10}{M[0,p-1]:>10}{aa:>10}{f:>8.2f}{fDE:>10.2f}{cov:>9.2f}{mark}")
    a1=consw(M,w,i1-1); a2=consw(M,w,i2-1)
    out[lab]=(n1,a1[2],n2,a2[2])

print("\n\n################  CONTROL EMPAREJADO  ################")
print(f"{'familia':<32}{'carboxilato i':>18}{'carboxilato i+2':>20}{'¿par completo?':>18}")
for lab,(n1,f1,n2,f2) in out.items():
    ok="SÍ" if min(f1,f2)>0.8 else "NO"
    print(f"{lab:<32}{f'{n1} {f1:.2f}':>18}{f'{n2} {f2:.2f}':>20}{ok:>18}")
