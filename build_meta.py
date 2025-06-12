#!/usr/bin/env python3
"""
构建元数据文件 all_meta.tsv
从现有数据结构提取样本信息，适配实际数据结构

注意：由于没有现成的tumor_summary文件夹，我们直接从spectra文件构建元数据
"""

import pandas as pd
import pathlib
from pathlib import Path
import re

def extract_sample_info_from_filename(filename):
    """从文件名提取样本信息和可能的肿瘤/正常分类"""
    # 移除前缀和后缀
    clean_name = filename.replace("pFind-Filtered_res_openHLA_", "").replace(".mgf_respFind.spectra", "")
    
    # 启发式判断肿瘤vs正常
    tumor_keywords = ["tumor", "cancer", "malignant", "HCC", "MM", "CPTAC"]
    normal_keywords = ["normal", "ctrl", "control", "Fib", "fibroblast"]
    
    sample_type = "Unknown"
    for keyword in tumor_keywords:
        if keyword.lower() in clean_name.lower():
            sample_type = "Tumor"
            break
    
    if sample_type == "Unknown":
        for keyword in normal_keywords:
            if keyword.lower() in clean_name.lower():
                sample_type = "Normal"
                break
    
    # 如果仍未确定，基于数据集特征判断
    if sample_type == "Unknown":
        # 对于细胞系，通常认为是tumor
        if any(cell_line in clean_name for cell_line in ["A375", "HCT116", "SupB15"]):
            sample_type = "Tumor"
        # 对于明确标记的对照组
        elif "ctrl" in clean_name.lower() or "control" in clean_name.lower():
            sample_type = "Normal"
        else:
            # 默认标记为Tumor（癌症研究中更常见）
            sample_type = "Tumor"
    
    return clean_name, sample_type

def build_metadata(root_path):
    """构建元数据文件"""
    root = Path(root_path)
    rows = []
    
    print("Scanning for human datasets...")
    
    # 扫描所有人类数据集
    for dataset_dir in root.glob("*_human"):
        dataset_name = dataset_dir.name.replace("_human", "")
        print(f"Processing dataset: {dataset_name}")
        
        # 扫描所有spectra文件
        for spectra_file in dataset_dir.glob("*.spectra"):
            sample_name, sample_type = extract_sample_info_from_filename(spectra_file.name)
            
            rows.append({
                "dataset": dataset_name,
                "sample": sample_name,
                "type": sample_type,
                "spectra": str(spectra_file),
                "batch_id": dataset_name,  # 使用dataset作为batch_id
                "HLA_alleles": "A*02:01,A*01:01,B*07:02"  # 默认HLA型号（仅A和B类），实际使用时需要更新
            })
    
    # 创建DataFrame
    meta_df = pd.DataFrame(rows)
    
    # 统计信息
    print(f"\nDataset summary:")
    print(f"Total samples: {len(meta_df)}")
    print(f"Datasets: {meta_df['dataset'].nunique()}")
    print(f"Sample type distribution:")
    print(meta_df['type'].value_counts())
    print(f"\nDataset breakdown:")
    print(meta_df.groupby(['dataset', 'type']).size().unstack(fill_value=0))
    
    return meta_df

def main():
    # 设置路径
    root_path = "/mnt/f/work/yang_ylab/cancer_datasets_links"
    
    # 构建元数据
    meta_df = build_metadata(root_path)
    
    # 保存文件
    output_file = Path(root_path) / "all_meta.tsv"
    meta_df.to_csv(output_file, sep="\t", index=False)
    print(f"\nMetadata saved to: {output_file}")
    
    # 创建一个示例HLA配置文件
    hla_config_file = Path(root_path) / "hla_alleles_config.txt"
    with open(hla_config_file, "w") as f:
        f.write("""# HLA Alleles Configuration
# Format: dataset_name:sample_pattern:HLA_alleles
# Example configurations (需要根据实际HLA typing结果更新)

# 默认配置 - 常见欧洲人群HLA型号（仅A和B类）
default:*:A*02:01,A*01:01,B*07:02

# 数据集特定配置示例（仅A和B类）
# PXD000394:HCC1143:A*01:01,A*03:01,B*08:01
# PXD000394:HCC1937:A*02:01,A*24:02,B*15:01

# 更多配置可以在这里添加...
""")
    print(f"HLA configuration template saved to: {hla_config_file}")
    
    print("\n" + "="*50)
    print("IMPORTANT NOTES:")
    print("1. HLA alleles are set to default values (A*02:01,A*01:01,B*07:02)")
    print("2. Only HLA-A and HLA-B classes are supported (HLA-C excluded)")
    print("3. Please update HLA_alleles column with actual HLA typing results")
    print("4. Sample type classification is heuristic - please verify manually")
    print("5. You can use OptiType, HLA-HD, or other tools for HLA typing from WES/RNA-seq data")
    print("="*50)

if __name__ == "__main__":
    main()