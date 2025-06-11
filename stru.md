# Spectra 数据集结构说明

## 项目概述

本项目是一个大规模的蛋白质组学数据分析仓库，包含多个数据集的质谱数据，主要聚焦于癌症（Cancer）和正常组织（Normal）样本的免疫肽组学分析。项目使用 pFind 软件进行肽段鉴定，并包含了对翻译后修饰（PTMs）的全面分析流程。

## 数据集元信息

### meta.txt.20250604 文件结构
- **文件格式**：制表符分隔的文本文件
- **列信息**：
  - `all_accession`：数据集编号（如 PXD020186）
  - `HLA(I/II)`：HLA 类型（I 型、II 型或混合）
  - `scenarios`：场景分类（Cancer、Normal、Autoimmune、Infection 等）
  - `Disease`：具体疾病名称（如 Lung Cancer、Melanoma 等）

### Cancer 数据集列表
根据 meta.txt.20250604，scenarios 为 Cancer 的数据集包括：
- PXD005084, PXD007596, PXD007860, PXD010372
- PXD013649, PXD014397, PXD015957, PXD021177
- PXD022020, PXD025716, PXD027766, PXD034820
- MSV000090437, PXD000394, PXD001898, PXD003790
- PXD004746, PXD004894, PXD005231, PXD006939
- PXD007635, PXD007935, PXD008984, PXD009602
- PXD009738, PXD009749, PXD009753, PXD009754
- PXD009755, PXD009925, PXD009935, PXD010450
- PXD010808, PXD011628, PXD012083, PXD012308
- PXD013057, PXD013831, PXD014017, PXD017824
- PXD018542, PXD020186, PXD020620, PXD022150
- PXD022949, PXD023044, PXD024871, PXD043989

## 目录结构

### 1. 数据集文件夹
**命名规则**：`[数据集编号]_[物种]/`（例如：`PXD020186_human/`）

**文件夹内容**：
- 包含 pFind 分析生成的 `.spectra` 文件
- 文件命名格式：`pFind-Filtered_res_openHLA_[样本信息].mgf_respFind.spectra`
- 每个文件夹可能包含多个技术重复（MS1、MS2、MS3）
- 不同的碎裂方法（CID、HCD）
- 不同的质量范围（400-650mz、400-1000mz）

### 2. tumor_summary 文件夹
**功能**：存储 Cancer 场景数据集的汇总信息

**文件格式**：`[数据集编号]_summary.csv`

**文件内容**：
- `File Name`：原始文件名
- `Sample Name`：样本名称
- `Total Unique Peptides`：总独特肽段数
- `Modified Peptides`：修饰肽段数量
- `Unmodified Peptides`：非修饰肽段数量
- `Modified Percentage`：修饰肽段百分比
- `Unmodified Percentage`：非修饰肽段百分比
- `Type`：样本类型（**Tumor** 或 **Normal**）- 这是最重要的列

## 文件格式详解

### .spectra 文件结构
**格式**：制表符分隔的文本文件

**主要列信息**：
- `Sequence`：肽段序列信息
- `Modification`：修饰信息
  - 空值：表示野生型肽段（无修饰）
  - 格式示例：`2,Oxidation[M]` - 表示序列第2位的M氨基酸发生氧化修饰
  - 多个修饰用分号(;)分隔：`2,Oxidation[M];5,Phospho[S]`
- 其他列：扫描信息、质量分数、蛋白质归属等

### 修饰信息解析规则
1. **位置标记**：数字表示修饰在序列中的位置（从1开始计数）
2. **修饰类型**：常见的修饰包括：
   - Oxidation[M]：甲硫氨酸氧化
   - Phospho[S/T/Y]：丝氨酸/苏氨酸/酪氨酸磷酸化
   - Acetyl[K]：赖氨酸乙酰化
   - Deamidated[N/Q]：天冬酰胺/谷氨酰胺脱酰胺化
3. **修饰目标**：方括号内为被修饰的氨基酸

## 样本类型识别

对于 Cancer 场景的数据集，通过 `tumor_summary/[数据集编号]_summary.csv` 文件中的 `Type` 列可以识别：
- **Tumor**：肿瘤样本
- **Normal**：正常样本（对照组）

这种分类对于后续的肿瘤vs正常样本对比分析至关重要。

## 数据处理流程

1. **原始质谱数据** → pFind 分析 → `.spectra` 文件
2. **样本级别汇总** → `tumor_summary/` 中的 CSV 文件
3. **修饰分析** → 识别和统计各类 PTMs
4. **差异分析** → 肿瘤 vs 正常样本的肽段和修饰差异

## 使用建议

1. **数据读取**：使用 pandas 读取 .spectra 文件时，建议使用 `sep='\t'` 参数
2. **修饰解析**：可以编写专门的函数解析 Modification 列，提取修饰位置和类型
3. **样本分组**：根据 tumor_summary 中的 Type 列进行样本分组分析
4. **批量处理**：可以根据 meta.txt.20250604 中的信息批量处理特定场景的数据集

## 注意事项

- 部分数据集可能同时包含肿瘤和正常样本
- 修饰信息的解析需要考虑多种修饰共存的情况
- 不同数据集的样本数量和质量可能存在差异，需要在分析时考虑