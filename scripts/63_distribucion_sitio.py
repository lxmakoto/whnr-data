"""Control cuantitativo: distancia minima entre carboxilatos en el dominio
pol-beta de nucleotidiltransferasas ACTIVAS frente a WhnR. El par catalitico se
localiza en cada estructura mapeando el HMM de Pfam NTP_transf_2 (estados 33/35)."""
import json, os, itertools, urllib.request, urllib.parse
import numpy as np, pyhmmer
from pyhmmer.easel import TextSequence, Alphabet
from Bio.PDB import MMCIFParser

AA3={"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
"HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
"THR":"T","TRP":"W","TYR":"Y","VAL":"V","MSE":"M"}
abc=Alphabet.amino()
with pyhmmer.plan7.HMMFile("data/pfam/Pfam-A.hmm") as f:
    for h in f:
        n=h.name.decode() if isinstance(h.name,bytes) else str(h.name)
        if n=="NTP_transf_2": HMM=h; break
CAT=(33,35)        # estados de coincidencia catalíticos (verificado con KNTasa/6NML)
P=MMCIFParser(QUIET=True)
os.makedirs("data/pdb",exist_ok=True)

def fetch(pid):
    p=f"data/pdb/{pid}.cif"
    if not os.path.exists(p):
        try: urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pid}.cif",p)
        except Exception: return None
    return p

def carb(r):
    ids=("OD1","OD2") if r.get_resname()=="ASP" else ("OE1","OE2")
    at=[a for a in r if a.get_id() in ids]
    return np.mean([a.coord for a in at],axis=0) if at else None

def analiza(pid):
    p=fetch(pid)
    if not p: return None
    try: st=P.get_structure(pid,p)[0]
    except Exception: return None
    best=None
    for ch in st:
        res=[r for r in ch if r.id[0] in (" ","H_MSE") and r.get_resname() in AA3]
        if len(res)<80: continue
        seq="".join(AA3[r.get_resname()] for r in res)
        s=TextSequence(name=b"x",sequence=seq).digitize(abc)
        hits=list(pyhmmer.hmmsearch([HMM],[s],E=1.0))
        if not hits or not len(hits[0]): continue
        d=hits[0][0].domains[0]
        tp=d.alignment.target_from-1; hp=d.alignment.hmm_from-1; mapa={}
        for tc,hc in zip(d.alignment.target_sequence,d.alignment.hmm_sequence):
            if tc!="-": tp+=1
            if hc!=".": hp+=1
            if tc!="-" and hc!=".": mapa[hp]=tp
        if CAT[0] not in mapa or CAT[1] not in mapa: continue
        i1,i2=mapa[CAT[0]],mapa[CAT[1]]
        # ventana de 116 residuos centrada en el par, igual tamaño que nuestro dominio
        lo=max(1,i1-58); hi=min(len(res),i1+58)
        win=res[lo-1:hi]
        pts=[(r.id[1],r.get_resname(),carb(r)) for r in win
             if r.get_resname() in ("ASP","GLU") and carb(r) is not None]
        if len(pts)<4: continue
        ds=[float(np.linalg.norm(a[2]-b[2])) for a,b in itertools.combinations(pts,2)]
        cat_d=None
        if res[i1-1].get_resname() in ("ASP","GLU") and res[i2-1].get_resname() in ("ASP","GLU"):
            ca,cb=carb(res[i1-1]),carb(res[i2-1])
            if ca is not None and cb is not None: cat_d=float(np.linalg.norm(ca-cb))
        cand=dict(pdb=pid,cadena=ch.id,n_carb=len(pts),dmin=min(ds),
                  cat=f"{res[i1-1].get_resname()}{res[i1-1].id[1]}/{res[i2-1].get_resname()}{res[i2-1].id[1]}",
                  cat_d=cat_d,todas=ds)
        if best is None or cand["dmin"]<best["dmin"]: best=cand
    return best

IDS=['1F5A','1FA0','1JAJ','1JMS','1JQR','1KDH','1KEJ','1NO5','1PX5','1Q78','1R89',
     '1SZ1','1TFW','1UET','1YLQ','2DR5','2DVI','2HHP','2Q66','6NML','2ZH2']
out=[]
for pid in IDS:
    r=analiza(pid)
    if r:
        out.append(r)
        cd=f"{r['cat_d']:.2f}" if r['cat_d'] else "n/d"
        print(f"  {pid}  cad {r['cadena']}  carboxilatos={r['n_carb']:>3}  "
              f"par catalítico {r['cat']:<14} d={cd:>6}  d_mín={r['dmin']:.2f} Å",flush=True)
    else:
        print(f"  {pid}  sin dominio localizable",flush=True)
json.dump([{k:v for k,v in r.items() if k!="todas"} for r in out],
          open("results/selection/ntase_activas.json","w"),indent=1)
np.save("/tmp/ntase_dists.npy",np.array([r["todas"] for r in out],dtype=object),allow_pickle=True)
json.dump([r["pdb"] for r in out],open("/tmp/ntase_ids.json","w"))
dm=[r["dmin"] for r in out]
print(f"\nactivas: n={len(out)}  d_mín mediana {np.median(dm):.2f} Å  rango {min(dm):.2f}–{max(dm):.2f}")
