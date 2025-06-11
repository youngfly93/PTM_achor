# 修饰位点 × HLA 锚定位点耦合分析 Pipeline

这是一个用于分析蛋白质翻译后修饰（PTM）与HLA-I类分子锚定位点耦合关系的完整分析pipeline，特别关注肿瘤样本与正常样本之间的差异。

## 🎯 分析目标

评估修饰是否preferentially出现在HLA-I锚定位点（P2/PΩ），并比较Tumor vs Normal差异。

## 📋 功能特性

- ✅ 自动化肽段提取与修饰解析
- ✅ HLA-I结合预测与锚位标注
- ✅ 修饰-锚位耦合分析
- ✅ Fisher精确检验统计分析
- ✅ Tumor vs Normal差异比较
- ✅ 多重检验校正（Bonferroni, FDR）
- ✅ 富集分析与可视化
- ✅ 支持MHCflurry预测（可选）

## 🚀 快速开始

### 1. 环境设置

```bash
# 基础环境（使用内置工具）
bash setup_env.sh

# 或者手动安装依赖
conda create -n anchor_mod python=3.10 -y
conda activate anchor_mod
conda install pandas pyteomics tqdm seaborn matplotlib scipy -y
pip install mhcflurry logomaker  # 可选：用于更精确的HLA预测
```

### 2. 数据准备

```bash
# 构建元数据文件
python3 build_meta_simple.py
```

这将生成：
- `all_meta.tsv`：包含所有样本的元数据
- `hla_alleles_config.txt`：HLA等位基因配置模板

### 3. 运行分析

```bash
# 测试模式（前5个样本）
python3 run_pipeline.py --test

# 完整分析
python3 run_pipeline.py

# 使用MHCflurry预测（需要先安装）
python3 run_pipeline.py --mhcflurry

# 自定义参数
python3 run_pipeline.py --meta custom_meta.tsv --output custom_results --limit 100
```

## 📊 输出文件

| 文件名 | 描述 |
|--------|------|
| `tumor_anchor_mod.tsv` | 肿瘤样本的修饰-锚位注释 |
| `normal_anchor_mod.tsv` | 正常样本的修饰-锚位注释 |
| `all_anchor_mod.tsv` | 所有样本的修饰-锚位注释 |
| `anchor_mod_enrichment.tsv` | 修饰在锚位的富集分析结果 |
| `anchor_mod_group_diff.csv` | Tumor vs Normal差异分析 |
| `tumor_vs_normal_report.txt` | 详细比较报告 |
| `pipeline_summary.txt` | 运行总结 |

## 🔬 分析流程

### 第一步：数据提取
- 从`.spectra`文件提取8-11mer肽段
- 解析修饰信息（位置、类型）
- 过滤和去重

### 第二步：HLA结合预测
- 使用MHCflurry或简化算法预测结合
- 标注P2和C-terminal锚位
- 计算结合评分

### 第三步：耦合分析
- 分析修饰位置与锚位的重叠
- 生成修饰-锚位耦合记录
- 按样本类型分组

### 第四步：统计分析
- Fisher精确检验评估富集
- 多重检验校正
- 效应量计算（Odds Ratio）

### 第五步：组间比较
- Tumor vs Normal差异分析
- 生成可视化报告
- 显著性评估

## 📁 项目结构

```
cancer_datasets_links/
├── setup_env.sh              # 环境设置脚本
├── build_meta_simple.py      # 元数据构建工具
├── extract_peptides.py       # 肽段提取模块
├── predict_binding.py        # HLA结合预测模块
├── anchor_coupling.py        # 修饰-锚位耦合分析
├── stats_plot.py            # 统计检验和可视化
├── compare_groups.py         # 组间差异分析
├── run_pipeline.py          # 主运行脚本
├── all_meta.tsv             # 样本元数据
├── results/                 # 分析结果目录
└── [dataset]_human/         # 原始数据目录
    └── *.spectra           # pFind分析结果
```

## ⚙️ 配置选项

### HLA等位基因配置
更新`all_meta.tsv`中的`HLA_alleles`列，或使用HLA typing工具：

```bash
# 示例HLA配置
A*02:01,B*07:02,C*07:02  # 欧洲人群常见型号
A*01:01,B*08:01,C*07:01  # 其他常见组合
```

### 预测方法选择

1. **简化预测**（默认）：基于肽段长度的快速评分
2. **MHCflurry预测**：基于机器学习的精确HLA结合预测

## 📈 结果解读

### 富集分析
- **Anchor Percentage**：修饰在锚位的百分比
- **Enrichment Ratio**：富集倍数
- **P-value**：Fisher精确检验p值

### 组间比较
- **T_anchor_rate**：肿瘤组锚位修饰比例
- **N_anchor_rate**：正常组锚位修饰比例
- **Odds Ratio**：效应量
- **P_group**：组间差异p值

## 🔧 进阶使用

### 自定义锚位规则
在`predict_binding.py`中修改`get_hla_motif_anchors`函数：

```python
def get_hla_motif_anchors(allele, peptide_length):
    # 自定义等位基因特异性锚位
    if allele.startswith('A*02'):
        return {2, peptide_length, peptide_length-1}
    # 添加更多规则...
```

### 修饰类型过滤
在运行前过滤特定修饰类型：

```python
from extract_peptides import filter_by_modifications
filtered_peptides = filter_by_modifications(peptides, ['Phospho', 'Oxidation'])
```

## 🚨 注意事项

1. **HLA等位基因**：默认使用通用型号，建议使用实际HLA typing结果
2. **样本分类**：自动分类可能不准确，请手动验证`all_meta.tsv`中的`type`列
3. **内存使用**：大数据集可能需要较多内存，可使用`--limit`参数分批处理
4. **统计功效**：小样本量可能影响统计检验结果

## 📚 参考文献

- MHCflurry: https://github.com/openvax/mhcflurry
- NetMHCpan: https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.1
- HLA锚位分析相关文献请参考免疫肽组学领域研究

## 🤝 贡献

欢迎提出问题和改进建议！

## 📄 许可证

本项目用于学术研究目的。