# RMVP 批处理脚本过滤逻辑检查与建议

## 结论（先说重点）

你的直觉有一部分是对的：

- **按单性状样本提取后再过滤（MAF/缺失率）**，在统计意义上更合理。
- 但你当前脚本在单性状环节用的这一条 `plink2` 命令里：
  `--keep + --extract + --maf + --geno --make-bed`，其中 `--maf/--geno` 实际就是在 `--keep` 后的样本子集上计算的，**并不是先按全体样本算完再套到单性状**。
- 真正可能“过早过滤”的是前面的**全局预筛选**（`GLOBAL_SNPLIST` 的 MAF/geno 以及全局 het 过滤），它会在单性状阶段之前先把一部分位点永久排除。

## 当前流程里最需要注意的点

1. 你先在全体样本上构建 `GLOBAL_SNPLIST`（MAF>=0.01, missing<=0.20）。
2. 再基于该列表做全局 het 过滤，得到 `HET_KEEP`。
3. 单性状时 `--extract $EXTRACT_LIST`（通常是 `HET_KEEP` 或 `GLOBAL_SNPLIST`）后，再做本性状 `--maf/--geno`。

这意味着：

- 某个位点如果在全体样本 `MAF<0.01`，但在某个性状子集中 `MAF` 其实不低，**也会被前面的全局步骤提前丢掉**。
- 所以如果你的目标是“每个性状独立决定保留哪些位点”，当前的全局 MAF/geno 预筛选会带来偏保守。

## 推荐改法（优先级）

### 方案 A（最贴合你的诉求）

- 保留全局 het 过滤（可选）。
- **取消全局 MAF/geno 预筛选作为硬门槛**。
- 单性状时：
  1) 先 `--keep` 提取样本子集（可同时 `--extract` 到“仅het保留列表”）。
  2) 再在该子集上执行 `--maf/--geno` 过滤并输出 `bed`。

这样就实现了你说的“先提取该性状材料，再过滤”。

### 方案 B（兼顾速度）

- 保留你当前全局预筛（提速）。
- 但把全局阈值设得更宽（例如 `MIN_MAF_BASE=0.001`, `MAX_GENO_BASE=0.5`），主要作为降维，不作为严格过滤。
- 最终严格阈值仍由单性状 `--maf/--geno` 决定。

## 两阶段命令模板（建议替换单性状部分）

```bash
# 1) 先按性状样本提取（可选叠加het列表）
plink2 \
  --pfile "$BASE/base" \
  --allow-extra-chr \
  --keep "$keep" \
  --extract "$HET_KEEP" \
  --make-pgen \
  --out "$trait_path/geno/${trait_dir}.subset" \
  --threads "$THREADS"

# 2) 再在该子集上做本性状过滤
plink2 \
  --pfile "$trait_path/geno/${trait_dir}.subset" \
  --allow-extra-chr \
  --maf "$maf_thr" \
  --geno "$geno_thr" \
  --make-bed \
  --out "$trait_path/geno/$trait_dir" \
  --threads "$THREADS"
```

> 说明：你现在“一步命令”与这个“两步命令”在 MAF/geno 的统计样本上本质一致；两步法的优势是逻辑更清晰、结果更容易审计。

## 小样本阈值策略评价

你脚本中的这几条策略总体是合理的：

- `maf = max(MIN_MAF_BASE, MAC_MIN/(2N))`：对小样本防止过低 MAC，合理。
- 小样本放宽 `--geno`（N<200 用0.30，N<100 用0.40）：工程上常用。
- `MIN_N_BUILD_GENO=60`：避免极小样本浪费计算，合理。

可以再加一个保护：若过滤后位点数过少（如 <5000），记录 warning，避免后续 GWAS 不稳定。

## 实操建议（简洁版）

- 若你追求“统计最严格按性状独立过滤”：
  - 关掉全局 MAF/geno 硬筛，保留或放宽全局 het。
  - 单性状执行提取后过滤。
- 若你追求“速度优先”：
  - 保留全局预筛，但阈值放宽，不要太激进。

