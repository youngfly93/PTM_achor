# 修饰位点 × HLA 锚定位点耦合分析 Pipeline

> 目标：评估修饰是否 preferentially 出现在 HLA‑I 锚定位点（P2/PΩ），并比较 Tumor vs Normal 差异。

---

## 0. 依赖环境

```bash
conda create -n anchor_mod python=3.10 pandas pyteomics tqdm mhcflurry netMHCpan seaborn matplotlib scipy logomaker
conda activate anchor_mod
```

* `netMHCpan ≥ 4.1` 需在 `$PATH`（或改用 `mhcflurry-predict`）。

---

## 1. 组织元数据（一次性）

```python
# build_meta.py
import pandas as pd, pathlib

root = pathlib.Path("/path/to/project")
rows = []
for summ in root.glob("tumor_summary/*_summary.csv"):
    ds = summ.stem.replace("_summary", "")
    df = pd.read_csv(summ)
    for _, r in df.iterrows():
        rows.append({
            "dataset": ds,
            "sample":  r["Sample Name"],
            "type":    r["Type"].strip(),        # Tumor / Normal
            "spectra": str(root / f"{r['File Name']}.spectra")
        })
meta = pd.DataFrame(rows).query("type in ['Tumor','Normal']")
meta.to_csv("all_meta.tsv", sep="\t", index=False)
```

---

## 2. 提取 8–11 mer 肽段与修饰位点

```python
# extract_peptides.py
import pandas as pd, re

AA = re.compile(r"[A-Z]")

def parse_mod(s):
    """'2,Oxidation[M];5,Phospho[S]' → [(2,'Oxidation[M]'), …]"""
    if pd.isna(s) or not s:
        return []
    return [ (int(x.split(",")[0]), ",".join(x.split(",")[1:])) for x in s.split(";") ]

def load_one(path):
    df = pd.read_csv(path, sep="\t", usecols=["Sequence", "Modification"], low_memory=False)
    df["length"] = df.Sequence.str.len()
    df = df.query("length >= 8 and length <= 11")
    df["mod_list"] = df.Modification.apply(parse_mod)
    return df[["Sequence", "mod_list"]].drop_duplicates()
```

---

## 3. 预测 HLA 结合并标注锚位

### 3.1 获取样本 HLA‑I 型号

* 在 `all_meta.tsv` 添加列 `HLA_alleles`，格式如 `A*02:01,B*15:01,C*07:02`。
* 可用 OptiType 等工具从 WES/RNA‑seq 推断。

### 3.2 批量预测与锚位标注

```python
# predict_binding.py
from mhcflurry import Class1PresentationPredictor
pred = Class1PresentationPredictor.load()

def top_allele(seq, alleles):
    res = pred.predict(peptides=[seq], alleles=alleles)
    best = res.sort_values("presentation_score", ascending=False).iloc[0]
    return best["allele"], best["presentation_score"]

def annotate_anchor(seq, allele):
    """简化：P2 与 C‑terminal 为 primary anchors"""
    return {2, len(seq)}
```

---

## 4. 构建修饰‑锚位耦合记录

```python
# anchor_coupling.py
import pandas as pd
from predict_binding import top_allele, annotate_anchor

def tag_anchor_mod(df, alleles):
    recs = []
    for _, r in df.iterrows():
        allele, score = top_allele(r.Sequence, alleles)
        anchors = annotate_anchor(r.Sequence, allele)
        for pos, mod in r.mod_list:   # 可能多修饰
            recs.append({
                "seq":   r.Sequence,
                "allele": allele,
                "score":  score,
                "mod":    mod.split("[")[0],  # Oxidation / Phospho …
                "anchor_tag": "anchor" if pos in anchors else "non_anchor"
            })
    return pd.DataFrame(recs)
```

---

## 5. 统计检验与可视化

```python
# stats_plot.py
import pandas as pd, scipy.stats as ss, seaborn as sns, matplotlib.pyplot as plt

def enrichment(df):
    tbl = df.groupby(["mod", "anchor_tag"]).size().unstack(fill_value=0)
    rows = []
    for m in tbl.index:
        a, n = tbl.loc[m, "anchor"], tbl.loc[m, "non_anchor"]
        rest_a = tbl["anchor"].sum() - a
        rest_n = tbl["non_anchor"].sum() - n
        p = ss.fisher_exact([[a, n], [rest_a, rest_n]])[1]
        rows.append({"mod": m, "anchor": a, "non_anchor": n, "p": p})
    return pd.DataFrame(rows).sort_values("p")

def violin(df, mod_type):
    sns.violinplot(data=df.query("mod == @mod_type"),
                   x="anchor_tag", y="score", inner="box")
    plt.title(f"{mod_type}: anchor vs non‑anchor binding score")
    plt.savefig(f"{mod_type}_anchor_violin.png", dpi=300)
```

---

## 6. Tumor vs Normal 差异

```python
# compare_groups.py
import pandas as pd, scipy.stats as ss

def compare_groups(df_tumor, df_normal):
    def cnt(df):
        return (df.groupby(["mod", "anchor_tag"]).size()
                  .unstack(fill_value=0).reset_index())
    t = cnt(df_tumor).rename(columns={"anchor": "t_anchor", "non_anchor": "t_non"})
    n = cnt(df_normal).rename(columns={"anchor": "n_anchor", "non_anchor": "n_non"})
    m = pd.merge(t, n, on="mod", how="outer").fillna(0)
    pvals = []
    for _, r in m.iterrows():
        table = [[r.t_anchor, r.t_non], [r.n_anchor, r.n_non]]
        pvals.append(ss.fisher_exact(table)[1])
    m["p_group"] = pvals
    return m.sort_values("p_group")
```

---

## 7. 主运行脚本

```python
# run_pipeline.py
import pandas as pd, pathlib, tqdm
from extract_peptides import load_one
from anchor_coupling import tag_anchor_mod
from compare_groups import compare_groups

meta = pd.read_csv("all_meta.tsv", sep="\t")

rows_t, rows_n = [], []
for _, r in tqdm.tqdm(meta.iterrows(), total=len(meta)):
    alleles = r.HLA_alleles.split(",")
    df_pep = load_one(r.spectra)
    df_tag = tag_anchor_mod(df_pep, alleles)
    df_tag["sample"] = r.sample
    df_tag["group"]  = r.type
    (rows_t if r.type == "Tumor" else rows_n).append(df_tag)

tumor   = pd.concat(rows_t)
normal  = pd.concat(rows_n)

tumor.to_parquet("tumor_anchor_mod.parquet")
normal.to_parquet("normal_anchor_mod.parquet")

diff = compare_groups(tumor, normal)
diff.to_csv("anchor_mod_group_diff.csv", index=False)

# 绘图示例：磷酸化
from stats_plot import violin
violin(pd.concat([tumor, normal]), "Phospho")
```

---

## 8. 产出文件

| 文件                                                       | 说明                            |
| -------------------------------------------------------- | ----------------------------- |
| `all_meta.tsv`                                           | 样本 → `.spectra` 对照表 + HLA 型号  |
| `tumor_anchor_mod.parquet` / `normal_anchor_mod.parquet` | 每条肽的修饰‑锚位注释                   |
| `anchor_mod_group_diff.csv`                              | Tumor vs Normal Fisher P 值、比值 |
| `Phospho_anchor_violin.png` 等                            | 视觉化示例                         |

---

## 9. 可选进阶

1. **等位基因特异锚位**：使用 NetMHCpan positional weight matrix 动态确定 anchor。
2. **长度归一化**：将 8–11 mer 的 P2/PΩ 归为 *Start* / *End* 等类别。
3. **Sequence logo**：`logomaker` 展示修饰锚 vs 非修饰锚位点偏好。
4. **Motif 聚类**：NMDS + DBSCAN 复现非正典抗原文献中的 motif‑sharing 图。

---

## 10. 一句话指令（给 Claude Code）

> **读取 `all_meta.tsv`，预测 8–11 mer 修饰肽 HLA 结合，标记 P2/PΩ 锚位，统计 Tumor vs Normal 中各修饰在锚位的富集（Fisher exact），生成 `anchor_mod_group_diff.csv`，并为 `Phospho` 绘制 binding‑score violin‑box 图。**
