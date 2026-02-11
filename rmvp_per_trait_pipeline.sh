#!/usr/bin/env bash
set -euo pipefail

# Per-trait genotype extraction + filtering + RMVP scaffold generation.
# Filters are applied AFTER trait-specific sample extraction.

VCF="${VCF:-lc7733all.PASS.ACGT.GTonly.merged.vcf.gz}"
PHE="${PHE:-lc3k.txt}"
OUT="${OUT:-lc3kRMVP_BATCH}"
THREADS="${THREADS:-80}"

# Per-trait hard filters requested by user
MAF_THR="${MAF_THR:-0.05}"          # remove MAF < 0.05
GENO_THR="${GENO_THR:-0.20}"        # remove missing rate > 0.20
HET_RATE_THR="${HET_RATE_THR:-0.10}"# remove heterozygosity rate > 0.10

MIN_N_BUILD_GENO="${MIN_N_BUILD_GENO:-60}"

need_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "[ERROR] missing command: $1" >&2; exit 1; }; }
for c in plink2 python3 awk sed wc date; do need_cmd "$c"; done

ts() { date '+%F %T'; }
step() { echo; echo "[$(ts)] [STEP] $*"; }
info() { echo "[$(ts)] [INFO] $*"; }
warn() { echo "[$(ts)] [WARN] $*" >&2; }

mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"
BASE="$OUT/_BASE"
TRAITS="$OUT/traits"
mkdir -p "$BASE" "$TRAITS"

info "VCF=$VCF"
info "PHE=$PHE"
info "OUT=$OUT"
info "THREADS=$THREADS"
info "Per-trait filters: MAF>=$MAF_THR missing<=$GENO_THR het_rate<=$HET_RATE_THR"

step "Convert VCF -> PLINK2 PGEN (once, no global variant filtering)"
if [[ ! -s "$BASE/base.pgen" || ! -s "$BASE/base.pvar" || ! -s "$BASE/base.psam" ]]; then
  plink2 \
    --vcf "$VCF" \
    --double-id \
    --allow-extra-chr \
    --max-alleles 2 \
    --snps-only just-acgt \
    --set-missing-var-ids @:#:\$r:\$a \
    --make-pgen \
    --out "$BASE/base" \
    --threads "$THREADS"
else
  info "Found existing base pgen; skip conversion."
fi

step "Parse phenotypes & create per-trait folders (keep/phe/run_rMVP.R)"
TRAITS_TSV="$OUT/traits.tsv"

python3 - "$PHE" "$TRAITS" "$TRAITS_TSV" <<'PY'
import sys, os, re, math
from statistics import NormalDist

phe_path, out_root, traits_tsv = sys.argv[1], sys.argv[2], sys.argv[3]
MISS = {"", "NA", "N/A", "NAN", "NULL", ".", "-9", "-99", "-999", "-9999"}

def is_missing(x):
    if x is None: return True
    s = str(x).strip()
    return s == "" or s.upper() in MISS

def safe_name(s):
    s = re.sub(r"[\/\\\:\*\?\"\<\>\|\s]+", "_", s.strip())
    s = re.sub(r"[^0-9A-Za-z_\-.]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s or "Trait"

def try_float(x):
    try: return float(x)
    except: return None

def detect_delim(first_line):
    if "\t" in first_line: return "\t"
    if "," in first_line and first_line.count(",") >= 2: return ","
    return None

def rank_INT(vals):
    n = len(vals)
    if n < 3: return [None]*n
    idx = sorted(range(n), key=lambda i: vals[i])
    ranks=[0.0]*n
    i=0
    while i<n:
      j=i
      while j+1<n and vals[idx[j+1]]==vals[idx[i]]: j+=1
      avg=(i+j+2)/2.0
      for k in range(i,j+1): ranks[idx[k]]=avg
      i=j+1
    nd = NormalDist()
    return [nd.inv_cdf(min(max((r-0.5)/n,1e-12),1-1e-12)) for r in ranks]

with open(phe_path, "r", encoding="utf-8", errors="replace") as f:
    head = f.readline()
    if not head:
        print("[ERROR] phenotype file empty", file=sys.stderr); sys.exit(2)
    delim = detect_delim(head)
    header = head.strip("\n").split(delim) if delim else re.split(r"\s+", head.strip())
    if len(header) < 2:
        print("[ERROR] phenotype header must have >=2 columns", file=sys.stderr); sys.exit(2)

    traits = header[1:]
    ids=[]
    cols=[[] for _ in traits]
    for line in f:
        line=line.strip("\n")
        if not line.strip():
            continue
        parts = line.split(delim) if delim else re.split(r"\s+", line.strip())
        sid = parts[0].strip() if parts else ""
        if not sid: continue
        ids.append(sid)
        vals = parts[1:] + [""]*(len(traits)-max(0,len(parts)-1))
        for j in range(len(traits)):
            cols[j].append(vals[j].strip() if j < len(vals) else "")

os.makedirs(out_root, exist_ok=True)
used={}
rows=[]
for j,tname in enumerate(traits):
    raw = cols[j]
    nm_idx=[i for i,x in enumerate(raw) if not is_missing(x)]
    N=len(nm_idx)
    tdir=safe_name(tname)
    if tdir in used:
        used[tdir]+=1
        tdir=f"{tdir}_{used[tdir]}"
    else:
        used[tdir]=0

    trait_dir=os.path.join(out_root, tdir)
    os.makedirs(trait_dir, exist_ok=True)

    with open(os.path.join(trait_dir,"keep.samples.txt"),"w",encoding="utf-8") as w:
        for i in nm_idx:
            sid=ids[i]
            w.write(f"{sid}\t{sid}\n")

    vals=[raw[i] for i in nm_idx]
    fv=[try_float(v) for v in vals]
    numeric_rate=sum(v is not None for v in fv)/(N if N else 1)

    phe_path_out=os.path.join(trait_dir,f"{tdir}.phe.tsv")
    n_cols=1
    ttype="missing_all"
    if N==0:
        with open(phe_path_out,"w",encoding="utf-8") as w:
            w.write("Taxa\t"+tname+"\n")
        ttype="missing_all"
    elif numeric_rate>=0.95:
        eff=[(ids[nm_idx[k]], fv[k]) for k in range(N) if fv[k] is not None]
        if len(eff)<3:
            ttype="too_few_numeric"
            with open(phe_path_out,"w",encoding="utf-8") as w:
                w.write("Taxa\t"+tname+"\n")
                for sid,v in eff:
                    w.write(f"{sid}\t{v:.6g}\n")
        else:
            eids=[x[0] for x in eff]
            evals=[x[1] for x in eff]
            z=rank_INT(evals)
            with open(phe_path_out,"w",encoding="utf-8") as w:
                w.write(f"Taxa\t{tname}\t{tname}__INT\n")
                for sid,v,zi in zip(eids,evals,z):
                    w.write(f"{sid}\t{v:.6g}\t{zi:.6g}\n")
            n_cols=2
            ttype="numeric"
            with open(os.path.join(trait_dir,"keep.samples.txt"),"w",encoding="utf-8") as w:
                for sid in eids:
                    w.write(f"{sid}\t{sid}\n")
    else:
        cats=sorted(set(vals))
        cmap={c:i+1 for i,c in enumerate(cats)}
        with open(os.path.join(trait_dir,"category_map.tsv"),"w",encoding="utf-8") as w:
            w.write("category\tcode\n")
            for c in cats: w.write(f"{c}\t{cmap[c]}\n")
        with open(phe_path_out,"w",encoding="utf-8") as w:
            w.write("Taxa\t"+tname+"\n")
            for i in nm_idx:
                w.write(f"{ids[i]}\t{cmap[raw[i]]}\n")
        ttype="categorical_str"

    with open(os.path.join(trait_dir,"run_rMVP.R"),"w",encoding="utf-8") as w:
        w.write(f'''library(rMVP)\nlibrary(bigmemory)\n\nbed_prefix <- "geno/{tdir}"\nphe_file   <- "{tdir}.phe.tsv"\n\nMVP.Data(\n  fileBed = bed_prefix,\n  filePhe = phe_file,\n  fileKin = FALSE,\n  filePC  = FALSE,\n  out     = "mvp.bed"\n)\n\ngenotypic_dat <- attach.big.matrix("mvp.bed.geno.desc")\nphenotype_dat <- read.table("mvp.bed.phe", head = TRUE, check.names = FALSE)\nmap_info      <- read.table("mvp.bed.geno.map", head = TRUE)\ndir.create("GWAS_out", showWarnings = FALSE)\nfor(i in 2:ncol(phenotype_dat)){{\n  trait_name <- colnames(phenotype_dat)[i]\n  rMVP::MVP(\n    phe  = phenotype_dat[, c(1, i)],\n    geno = genotypic_dat,\n    map  = map_info,\n    nPC.GLM = 10, nPC.MLM = 10, nPC.FarmCPU = 10,\n    maxLine = 100000,\n    ncpus = as.integer(Sys.getenv("NCPUS", "80")),\n    vc.method = "BRENT", maxLoop = 50, method.bin = "static", threshold = 0.05,\n    method = c("GLM", "MLM", "FarmCPU"),\n    memo = trait_name, outpath = "GWAS_out"\n  )\n  gc()\n}}\n''')

    with open(os.path.join(trait_dir,"trait_info.tsv"),"w",encoding="utf-8") as w:
        w.write("trait_orig\ttrait_dir\tN_raw_nonmissing\ttype\tn_cols\n")
        w.write(f"{tname}\t{tdir}\t{N}\t{ttype}\t{n_cols}\n")

    rows.append((tname,tdir,str(N),ttype,str(n_cols)))

with open(traits_tsv,"w",encoding="utf-8") as w:
    w.write("trait_orig\ttrait_dir\tN_raw_nonmissing\ttype\tn_cols\n")
    for r in rows:
        w.write("\t".join(r)+"\n")

print(f"[INFO] Wrote traits list: {traits_tsv} (n={len(rows)})", file=sys.stderr)
PY

step "Per-trait genotype extraction + per-trait filtering"
mapfile -t LINES < <(tail -n +2 "$TRAITS_TSV" || true)
NT=${#LINES[@]}
info "Total traits: $NT"
if [[ "$NT" -eq 0 ]]; then
  warn "No traits found in $TRAITS_TSV"
  exit 0
fi

for ((k=0; k<NT; k++)); do
  IFS=$'\t' read -r trait_orig trait_dir N_raw trait_type n_cols <<< "${LINES[$k]}"
  trait_path="$TRAITS/$trait_dir"
  keep="$trait_path/keep.samples.txt"
  info "[TRAIT $((k+1))/$NT] $trait_orig => $trait_dir (N=$N_raw, type=$trait_type)"

  if [[ ! -s "$keep" ]]; then
    warn "Missing keep file: $keep ; skip"
    continue
  fi
  if [[ "$N_raw" -lt "$MIN_N_BUILD_GENO" ]]; then
    warn "N=$N_raw < MIN_N_BUILD_GENO=$MIN_N_BUILD_GENO ; skip genotype build"
    continue
  fi

  mkdir -p "$trait_path/geno"
  subset="$trait_path/geno/${trait_dir}.subset"
  outpref="$trait_path/geno/${trait_dir}"
  basic="$trait_path/geno/${trait_dir}.pass_basic"
  gcount="$trait_path/geno/${trait_dir}.gcount"
  hetkeep="$trait_path/geno/${trait_dir}.pass_het.snplist"

  # Step A: subset samples first (no variant filtering)
  plink2 \
    --pfile "$BASE/base" \
    --allow-extra-chr \
    --keep "$keep" \
    --make-pgen \
    --out "$subset" \
    --threads "$THREADS"

  # Step B: per-trait MAF + missing
  plink2 \
    --pfile "$subset" \
    --allow-extra-chr \
    --maf "$MAF_THR" \
    --geno "$GENO_THR" \
    --write-snplist \
    --out "$basic" \
    --threads "$THREADS"

  if [[ ! -s "$basic.snplist" ]]; then
    warn "No variants survived MAF/geno for $trait_dir"
    : > "$hetkeep"
  else
    # Step C: heterozygosity filter (variant-level het rate)
    plink2 \
      --pfile "$subset" \
      --allow-extra-chr \
      --extract "$basic.snplist" \
      --geno-counts \
      --out "$gcount" \
      --threads "$THREADS"

    python3 - "$gcount.gcount" "$HET_RATE_THR" "$hetkeep" <<'PY'
import sys, re
infile, thr, out = sys.argv[1], float(sys.argv[2]), sys.argv[3]

def norm(s):
    return re.sub(r'[^A-Za-z0-9]+','',s).upper()

with open(infile, 'r', encoding='utf-8', errors='replace') as f:
    header=f.readline().strip()
    cols=re.split(r'[\t ]+', header)
    idx={norm(c):i for i,c in enumerate(cols)}

    def find(names):
        for n in names:
            nn=norm(n)
            if nn in idx: return idx[nn]
        for c,i in idx.items():
            for n in names:
                if norm(n) in c: return i
        return None

    id_i=find(["ID","SNP","VARIANT","RSID"])
    het_i=find(["HET_CT","HET","N_HET","HETCOUNT"])
    homr_i=find(["HOM_REF_CT","HOMREF_CT","N_HOM_REF","HOMREF"])
    homa_i=find(["HOM_ALT_CT","HOMALT_CT","N_HOM_ALT","HOMALT"])
    obs_i=find(["OBS_CT","NOBS","NONMISS_CT"])

    if id_i is None: id_i=0
    kept=0
    with open(out,'w',encoding='utf-8') as w:
        for line in f:
            p=re.split(r'[\t ]+', line.strip())
            if len(p)<=max(id_i, het_i if het_i is not None else 0):
                continue
            try:
                het=float(p[het_i])
                if obs_i is not None:
                    obs=float(p[obs_i])
                else:
                    homr=float(p[homr_i]) if homr_i is not None else 0.0
                    homa=float(p[homa_i]) if homa_i is not None else 0.0
                    obs=homr+homa+het
            except:
                continue
            if obs <= 0: continue
            if het/obs <= thr:
                w.write(p[id_i] + "\n")
                kept += 1
print(f"[INFO] pass_het variants={kept}", file=sys.stderr)
PY
  fi

  if [[ ! -s "$hetkeep" ]]; then
    warn "No variants survived heterozygosity filter for $trait_dir"
    continue
  fi

  # Step D: final per-trait bed/bim/fam
  plink2 \
    --pfile "$subset" \
    --allow-extra-chr \
    --extract "$hetkeep" \
    --make-bed \
    --out "$outpref" \
    --threads "$THREADS"

  ns=$(wc -l < "$outpref.fam" || echo 0)
  nv=$(wc -l < "$outpref.bim" || echo 0)
  {
    echo -e "trait_orig\t$trait_orig"
    echo -e "trait_dir\t$trait_dir"
    echo -e "N_raw_nonmissing\t$N_raw"
    echo -e "N_used\t$ns"
    echo -e "M_used\t$nv"
    echo -e "maf_threshold\t$MAF_THR"
    echo -e "geno_threshold\t$GENO_THR"
    echo -e "het_rate_threshold\t$HET_RATE_THR"
  } > "$trait_path/geno/geno_stats.tsv"

  info "Result: samples=$ns variants=$nv -> $outpref.(bed/bim/fam)"
done

step "DONE"
info "Per-trait files ready. Run GWAS manually in each trait folder: NCPUS=80 Rscript run_rMVP.R"
