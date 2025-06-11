#!/usr/bin/env python3
"""
构建元数据文件 all_meta.tsv (简化版本，不依赖pandas)
"""

import os
import glob
from pathlib import Path

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
    rows = []
    dataset_counts = {}
    type_counts = {"Tumor": 0, "Normal": 0, "Unknown": 0}
    
    print("Scanning for human datasets...")
    
    # 扫描所有人类数据集
    for dataset_dir in glob.glob(os.path.join(root_path, "*_human")):
        if not os.path.isdir(dataset_dir):
            continue
            
        dataset_name = os.path.basename(dataset_dir).replace("_human", "")
        print(f"Processing dataset: {dataset_name}")
        
        dataset_counts[dataset_name] = {"Tumor": 0, "Normal": 0, "Unknown": 0}
        
        # 扫描所有spectra文件
        for spectra_file in glob.glob(os.path.join(dataset_dir, "*.spectra")):
            filename = os.path.basename(spectra_file)
            sample_name, sample_type = extract_sample_info_from_filename(filename)
            
            rows.append([
                dataset_name,
                sample_name,
                sample_type,
                spectra_file,
                "A*02:01,B*07:02,C*07:02"  # 默认HLA型号
            ])
            
            dataset_counts[dataset_name][sample_type] += 1
            type_counts[sample_type] += 1
    
    # 统计信息
    print(f"\nDataset summary:")
    print(f"Total samples: {len(rows)}")
    print(f"Datasets: {len(dataset_counts)}")
    print(f"Sample type distribution:")
    for type_name, count in type_counts.items():
        print(f"  {type_name}: {count}")
    
    print(f"\nDataset breakdown:")
    for dataset, counts in dataset_counts.items():
        print(f"  {dataset}: Tumor={counts['Tumor']}, Normal={counts['Normal']}, Unknown={counts['Unknown']}")
    
    return rows

def main():
    # 设置路径
    root_path = "/mnt/f/work/yang_ylab/cancer_datasets_links"
    
    # 构建元数据
    rows = build_metadata(root_path)
    
    # 保存文件
    output_file = os.path.join(root_path, "all_meta.tsv")
    with open(output_file, "w") as f:
        # 写入表头
        f.write("dataset\tsample\ttype\tspectra\tHLA_alleles\n")
        
        # 写入数据行
        for row in rows:
            f.write("\t".join(row) + "\n")
    
    print(f"\nMetadata saved to: {output_file}")
    
    # 创建一个示例HLA配置文件
    hla_config_file = os.path.join(root_path, "hla_alleles_config.txt")
    with open(hla_config_file, "w") as f:
        f.write("""# HLA Alleles Configuration
# Format: dataset_name:sample_pattern:HLA_alleles
# Example configurations (需要根据实际HLA typing结果更新)

# 默认配置 - 常见欧洲人群HLA型号
default:*:A*02:01,B*07:02,C*07:02

# 数据集特定配置示例
# PXD000394:HCC1143:A*01:01,B*08:01,C*07:01
# PXD000394:HCC1937:A*02:01,B*15:01,C*03:04

# 更多配置可以在这里添加...
""")
    print(f"HLA configuration template saved to: {hla_config_file}")
    
    print("\n" + "="*50)
    print("IMPORTANT NOTES:")
    print("1. HLA alleles are set to default values (A*02:01,B*07:02,C*07:02)")
    print("2. Please update HLA_alleles column with actual HLA typing results")
    print("3. Sample type classification is heuristic - please verify manually")
    print("4. You can use OptiType, HLA-HD, or other tools for HLA typing from WES/RNA-seq data")
    print("="*50)

if __name__ == "__main__":
    main()