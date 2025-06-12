# 批次分析功能说明

## 概述

本项目新增了按数据集（批次）进行分析的功能，支持：
- 按数据集分别计算修饰-锚位富集
- 使用随机效应模型进行meta分析
- 生成森林图展示效应量和异质性
- 评估批次间的异质性（I²和τ²）

## 新增文件

1. **batch_stats.py**: 批次级别的统计分析
   - `batch_enrich()`: 计算单批次的logOR和标准误
   - `batch_group_comparison()`: 批次内的组间比较
   - `summarize_batch_heterogeneity()`: 评估批次间异质性

2. **meta_merge.py**: Meta分析和可视化
   - `random_effect()`: DerSimonian-Laird随机效应模型
   - `fixed_effect()`: 固定效应模型
   - `forest_plot()`: 生成单个修饰的森林图
   - `create_summary_forest_plot()`: 生成汇总森林图

3. **run_pipeline_batch.py**: 支持批次分析的主程序
   - 新增`--mode`参数：`batch`（批次分析）或`traditional`（传统分析）
   - 自动按batch_id分组处理数据
   - 生成批次级别和meta级别的结果

## 使用方法

### 1. 更新元数据（添加batch_id列）

```bash
# 重新生成元数据，自动添加batch_id列
python3 build_meta.py
```

### 2. 运行批次分析

```bash
# 测试模式（5个样本）
python3 run_pipeline_batch.py --meta all_meta_test.tsv --test --mode batch

# 完整批次分析
python3 run_pipeline_batch.py --meta all_meta_updated.tsv --mode batch

# 使用MHCflurry预测
python3 run_pipeline_batch.py --meta all_meta_updated.tsv --mode batch --mhcflurry

# 传统模式（向后兼容）
python3 run_pipeline_batch.py --meta all_meta_updated.tsv --mode traditional
```

## 输出文件结构

```
results_batch/
├── batch_results/                      # 批次级别结果
│   ├── batch_PXD000394_tumor_enrichment.csv
│   ├── batch_PXD000394_normal_enrichment.csv
│   └── ...
├── forest_plots/                       # 森林图
│   ├── forest_Oxidation_tumor.png
│   ├── forest_Phospho_tumor.png
│   └── ...
├── per_batch_anchor_enrich.csv         # 所有批次的富集结果
├── meta_anchor_enrich_tumor.csv        # 肿瘤样本meta分析结果
├── meta_anchor_enrich_normal.csv       # 正常样本meta分析结果
├── meta_tumor_vs_normal.csv            # 组间比较meta分析
├── summary_forest_plot_tumor.png       # 汇总森林图
├── heterogeneity_report.txt            # 异质性评估报告
└── pipeline_summary.txt                # 运行总结
```

## 结果解读

### Meta分析结果（meta_anchor_enrich_tumor.csv）

| 列名 | 含义 |
|------|------|
| mod | 修饰类型 |
| n_batches | 包含该修饰的批次数 |
| OR | 合并的优势比 |
| OR_CI_lower/upper | 95%置信区间 |
| logOR | log优势比 |
| p_value | 显著性p值 |
| p_fdr | FDR校正后的p值 |
| I2 | 异质性百分比（0-100%） |
| tau2 | 批次间方差 |

### 异质性解读

- **I² < 25%**: 低异质性
- **I² 25-50%**: 中等异质性  
- **I² > 50%**: 高异质性
- **τ² > 0.1**: 实质性批次间差异

### 森林图解读

- 每个方块代表一个批次的效应量
- 方块大小反映权重（基于标准误）
- 菱形代表合并效应量
- 垂直虚线为零效应线（logOR=0）

## 高级功能

### 自定义批次定义

编辑元数据文件的`batch_id`列，可以按不同方式分组：
- 按数据集：`PXD000394`
- 按癌种：`BRCA_PXD000394`
- 按平台：`TMT_PXD000394`

### 调整异质性阈值

在`meta_merge.py`中修改：
```python
# 自定义高异质性阈值
HIGH_HETEROGENEITY_I2 = 75  # 默认50
```

### 批量生成森林图

```python
# 在run_pipeline_batch.py中修改
# 生成前10个显著修饰的森林图
for i, (_, row) in enumerate(meta_results.head(10).iterrows()):
    ...
```

## 注意事项

1. **最小批次数**：Meta分析至少需要2个批次
2. **样本量**：每个批次应有足够的样本量（建议>10）
3. **零值处理**：使用0.5校正避免零计数
4. **内存使用**：批次分析可能需要更多内存

## 故障排查

### 问题：森林图生成失败
- 检查是否安装matplotlib和seaborn
- 确保至少有2个批次包含该修饰

### 问题：Meta分析结果为空
- 检查批次数据是否正确生成
- 验证logOR和SE计算是否有效

### 问题：异质性过高
- 考虑按更细粒度分组（如癌种+平台）
- 检查是否有异常批次

## 引用

如使用批次分析功能，请引用：
- DerSimonian & Laird (1986) - 随机效应模型
- Higgins et al. (2003) - I²统计量