"""Recupera las CDS nucleotidicas de los 742 homologos para hacer dN/dS a escala
de superfilo, en lugar de los 8 genotipos de un solo orden.

Cadena: UniRef100_ACC -> UniParc UPI -> referencia cruzada EMBLWGS/EMBL
        -> protein_id -> NCBI efetch (rettype=fasta_cds_na) -> CDS.
Nota: no se envia ningun dato personal a NCBI; solo el parametro 'tool'.
"""
import json, os, time, urllib.request, urllib.parse, collections, sys
UA={"User-Agent":"cpr-dark-matter research"}
CACHE="results/selection/cds_cache.json"
cache=json.load(open(CACHE)) if os.path.exists(CACHE) else {"upi":{},"pid":{},"cds":{}}
def save(): json.dump(cache,open(CACHE,"w"))

def get(u,txt=False,tries=3):
    for a in range(tries):
        try:
            r=urllib.request.urlopen(urllib.request.Request(u,headers=UA),timeout=90)
            return r.read().decode() if txt else json.load(r)
        except Exception:
            if a==tries-1: return None
            time.sleep(2*(a+1))

accs=[l[1:].split("\t")[0].replace("UniRef100_","")
      for l in open("results/msa/cl33274978.a3m") if l.startswith(">UniRef100_")]
print(f"homólogos: {len(accs)}",flush=True)

# --- 1) accesion -> UPI ---
todo=[a for a in accs if a not in cache["upi"]]
print(f"UPI por resolver: {len(todo)}",flush=True)
for i,a in enumerate(todo,1):
    d=get(f"https://rest.uniprot.org/uniparc/search?query={a}&fields=upi&format=json")
    cache["upi"][a]=(d or {}).get("results",[{}])[0].get("uniParcId") if (d or {}).get("results") else None
    if i%50==0: print(f"   UPI {i}/{len(todo)}",flush=True); save()
save()
upis={a:u for a,u in cache["upi"].items() if u and a in set(accs)}
print(f"UPI resueltos: {len(upis)}/{len(accs)}",flush=True)

# --- 2) UPI -> protein_id de EMBL ---
todo=[u for u in set(upis.values()) if u not in cache["pid"]]
print(f"protein_id por resolver: {len(todo)}",flush=True)
for i,u in enumerate(todo,1):
    d=get(f"https://rest.uniprot.org/uniparc/{u}.json")
    pid=None
    for x in (d or {}).get("uniParcCrossReferences",[]):
        if x.get("database","").startswith("EMBL") and x.get("active"):
            pid=x.get("id"); break
    if pid is None:
        for x in (d or {}).get("uniParcCrossReferences",[]):
            if x.get("database","").startswith("EMBL"): pid=x.get("id"); break
    cache["pid"][u]=pid
    if i%50==0: print(f"   pid {i}/{len(todo)}",flush=True); save()
save()
pid_of={a:cache["pid"].get(u) for a,u in upis.items()}
have=[a for a,p in pid_of.items() if p]
print(f"con protein_id: {len(have)}/{len(accs)}",flush=True)

# --- 3) protein_id -> CDS (en lotes) ---
need=[pid_of[a] for a in have if pid_of[a] not in cache["cds"]]
need=list(dict.fromkeys(need))
print(f"CDS por descargar: {len(need)}",flush=True)
B=120
for i in range(0,len(need),B):
    batch=need[i:i+B]
    url=("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id="
         +",".join(batch)+"&rettype=fasta_cds_na&retmode=text&tool=cpr-dark-matter")
    txt=get(url,txt=True)
    if txt:
        name=None;buf=[]
        for line in txt.split("\n"):
            if line.startswith(">"):
                if name: cache["cds"][name]="".join(buf)
                m=[t for t in line.split("_cds_") if t]
                name=line.split("_cds_")[1].split("_")[0].split(".")[0] if "_cds_" in line else None
                buf=[]
            elif line.strip() and name: buf.append(line.strip())
        if name: cache["cds"][name]="".join(buf)
    print(f"   lote {i//B+1}/{(len(need)+B-1)//B}  acumuladas {len(cache['cds'])}",flush=True)
    save(); time.sleep(0.4)
print(f"\nCDS en caché: {len(cache['cds'])}")
json.dump({"pid_of":pid_of},open("results/selection/pid_map.json","w"))
print("-> results/selection/cds_cache.json, pid_map.json")
