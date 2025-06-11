#!/usr/bin/env python3
"""
更新HLA配置工具 - 根据HLA参考集更新元数据文件中的HLA等位基因
"""

import csv
import random
from hla_manager import HLAManager

def update_metadata_with_hla_reference(meta_file="all_meta.tsv", output_file="all_meta_updated.tsv"):
    """
    使用HLA参考集更新元数据文件中的HLA等位基因
    
    Args:
        meta_file: 输入元数据文件
        output_file: 输出元数据文件
    """
    # 初始化HLA管理器
    hla_manager = HLAManager()
    
    print("Updating HLA alleles in metadata file...")
    print(f"Available alleles: {len(hla_manager.get_all_alleles())}")
    
    # 准备不同人群的HLA配置
    population_configs = {
        'European': hla_manager.suggest_alleles_for_population('European'),
        'Asian': hla_manager.suggest_alleles_for_population('Asian'),
        'African': hla_manager.suggest_alleles_for_population('African')
    }
    
    print("\nPopulation-specific alleles:")
    for pop, alleles in population_configs.items():
        print(f"  {pop}: {alleles}")
    
    # 读取原始元数据
    try:
        with open(meta_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            metadata = list(reader)
    except FileNotFoundError:
        print(f"Error: File {meta_file} not found")
        return
    
    # 更新HLA等位基因
    updated_count = 0
    population_types = list(population_configs.keys())
    
    for record in metadata:
        # 随机分配人群类型（在实际应用中，这应该基于实际的人群信息）
        population = random.choice(population_types)
        suggested_alleles = population_configs[population]
        
        # 随机选择3-4个等位基因
        num_alleles = random.randint(3, min(4, len(suggested_alleles)))
        selected_alleles = random.sample(suggested_alleles, num_alleles)
        
        # 更新HLA等位基因字符串
        record['HLA_alleles'] = ','.join(selected_alleles)
        record['HLA_population'] = population  # 添加人群信息
        
        updated_count += 1
    
    # 写入更新后的元数据
    fieldnames = list(metadata[0].keys()) if metadata else []
    
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(metadata)
    
    print(f"\nUpdated {updated_count} records")
    print(f"Output saved to: {output_file}")
    
    # 统计更新后的HLA分布
    print("\nHLA allele distribution after update:")
    allele_counts = {}
    for record in metadata:
        alleles = record['HLA_alleles'].split(',')
        for allele in alleles:
            allele = allele.strip()
            allele_counts[allele] = allele_counts.get(allele, 0) + 1
    
    # 显示最常用的等位基因
    sorted_alleles = sorted(allele_counts.items(), key=lambda x: x[1], reverse=True)
    print("Top 10 most frequent alleles:")
    for allele, count in sorted_alleles[:10]:
        print(f"  {allele}: {count} samples")

def create_custom_hla_mapping(datasets, output_file="custom_hla_mapping.tsv"):
    """
    为特定数据集创建自定义HLA映射
    
    Args:
        datasets: 数据集列表
        output_file: 输出文件路径
    """
    hla_manager = HLAManager()
    
    # 为每个数据集分配不同的HLA配置
    dataset_configs = {}
    populations = ['European', 'Asian', 'African']
    
    for i, dataset in enumerate(datasets):
        population = populations[i % len(populations)]
        alleles = hla_manager.suggest_alleles_for_population(population)
        dataset_configs[dataset] = {
            'population': population,
            'alleles': alleles,
            'primary_alleles': alleles[:3]  # 主要等位基因
        }
    
    # 导出映射文件
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Dataset', 'Population', 'Primary_Alleles', 'All_Suggested_Alleles'])
        
        for dataset, config in dataset_configs.items():
            writer.writerow([
                dataset,
                config['population'],
                ','.join(config['primary_alleles']),
                ','.join(config['alleles'])
            ])
    
    print(f"Custom HLA mapping saved to: {output_file}")
    return dataset_configs

def main():
    """主函数"""
    print("HLA Configuration Update Tool")
    print("=" * 40)
    
    # 检查是否存在元数据文件
    import os
    if os.path.exists("all_meta.tsv"):
        print("Found all_meta.tsv, updating HLA alleles...")
        update_metadata_with_hla_reference()
        
        # 创建备份的简化版本（用于测试）
        print("\nCreating test version with first 10 samples...")
        
        with open("all_meta_updated.tsv", 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            test_records = list(reader)[:10]  # 只取前10个样本
        
        with open("all_meta_test.tsv", 'w', newline='', encoding='utf-8') as f:
            if test_records:
                writer = csv.DictWriter(f, fieldnames=test_records[0].keys(), delimiter='\t')
                writer.writeheader()
                writer.writerows(test_records)
        
        print("Test metadata saved to: all_meta_test.tsv")
        
    else:
        print("all_meta.tsv not found, creating example HLA mapping...")
        
        # 示例数据集
        example_datasets = [
            'MSV000090437', 'PXD000394', 'PXD001898', 'PXD003790', 'PXD004746'
        ]
        
        create_custom_hla_mapping(example_datasets)
    
    # 显示HLA统计信息
    hla_manager = HLAManager()
    hla_manager.export_statistics("hla_reference_stats.txt")

if __name__ == "__main__":
    main()