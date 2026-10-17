"""Pide al MISMO servidor MSA de ColabFold los homologos de nucleotidiltransferasas
pol-beta ACTIVAS, para comparar conservacion con el mismo codigo y la misma fuente."""
import time, io, tarfile, urllib.request, urllib.parse, os, sys, json

API="https://api.colabfold.com"
UA={"User-Agent":"cpr-dark-matter research (contact via project repo)"}

def submit(seq,name):
    data=urllib.parse.urlencode({"q":f">{name}\n{seq}","mode":"env"}).encode()
    r=urllib.request.urlopen(urllib.request.Request(API+"/ticket/msa",data=data,headers=UA),timeout=120)
    return json.load(r)

def wait(tid,maxmin=25):
    t0=time.time()
    while time.time()-t0 < maxmin*60:
        r=urllib.request.urlopen(urllib.request.Request(f"{API}/ticket/{tid}",headers=UA),timeout=60)
        st=json.load(r).get("status")
        if st in ("COMPLETE","ERROR"): return st
        time.sleep(15)
    return "TIMEOUT"

def fetch(tid,out):
    r=urllib.request.urlopen(urllib.request.Request(f"{API}/result/download/{tid}",headers=UA),timeout=300)
    buf=io.BytesIO(r.read())
    with tarfile.open(fileobj=buf,mode="r:gz") as tar:
        for m in tar.getmembers():
            if m.name.endswith(".a3m"):
                open(out,"wb").write(tar.extractfile(m).read()); return True
    return False

TARGETS={
 # nucleotidiltransferasa de kanamicina: pol-beta ACTIVA, catalíticas Asp50/Glu52
 "KNTase": open("data/complexes/ctrlpos_KNTase_ATP.yaml").read().split("sequence: ")[1].split("\n")[0],
}
for name,seq in TARGETS.items():
    out=f"results/control_ntase/{name}.a3m"
    if os.path.exists(out) and os.path.getsize(out)>1000:
        print(f"{name}: ya existe"); continue
    print(f"{name}: {len(seq)} aa -> enviando al servidor MSA de ColabFold",flush=True)
    j=submit(seq,name); tid=j.get("id")
    print(f"   ticket {tid} ({j.get('status')})",flush=True)
    st=wait(tid)
    print(f"   estado final: {st}",flush=True)
    if st=="COMPLETE" and fetch(tid,out):
        n=sum(1 for l in open(out) if l.startswith(">"))
        print(f"   -> {out} ({n} secuencias)")
    else:
        print("   fallo"); sys.exit(1)
