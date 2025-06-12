#!/usr/bin/env python3
"""
修饰位点 × HLA 锚定位点耦合分析 Pipeline - 批次分析版本
支持按数据集（batch）进行分析，然后进行meta分析
"""

import csv
import os
import sys
from pathlib import Path
import time
import pandas as pd
import numpy as np

# 导入自定义模块
from extract_peptides import load_one, get_modification_stats
from predict_binding import parse_allele_string, batch_predict_binding
from anchor_coupling import tag_anchor_modifications, export_coupling_records
from stats_plot import enrichment_analysis, export_enrichment_results, compare_anchor_vs_non_anchor, generate_text_violin_plot
from compare_groups import compare_groups, export_group_comparison_results, summarize_group_differences, create_comparison_report
from hla_manager import HLAManager
from batch_stats import batch_enrich, batch_group_comparison, export_batch_stats, summarize_batch_heterogeneity
from meta_merge import meta_analysis_by_mod, forest_plot, export_meta_results, create_summary_forest_plot


def load_metadata(meta_file):
    """
    加载元数据文件
    
    Args:
        meta_file: 元数据文件路径
    
    Returns:
        list: 元数据记录列表
    """
    metadata = []
    
    try:
        with open(meta_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                metadata.append(row)
    except FileNotFoundError:
        print(f"Error: Metadata file not found: {meta_file}")
        sys.exit(1)
    except Exception as e:
        print(f"Error loading metadata: {e}")
        sys.exit(1)
    
    return metadata


def process_sample(sample_record, use_mhcflurry=False, hla_manager=None):
    """
    处理单个样本（集成HLA管理器和batch_id）
    
    Args:
        sample_record: 样本元数据记录
        use_mhcflurry: 是否使用MHCflurry预测
        hla_manager: HLA管理器实例
    
    Returns:
        list: 该样本的耦合记录
    """
    if hla_manager is None:
        hla_manager = HLAManager()
    sample_name = sample_record['sample']
    spectra_path = sample_record['spectra']
    sample_type = sample_record['type']
    hla_alleles_str = sample_record['HLA_alleles']
    batch_id = sample_record.get('batch_id', sample_record['dataset'])  # 使用batch_id或dataset
    
    # 检查文件是否存在
    if not os.path.exists(spectra_path):
        print(f"        Warning: Spectra file not found, skipping")
        return []
    
    # 解析HLA等位基因（使用HLA管理器）
    alleles = parse_allele_string(hla_alleles_str, hla_manager)
    if not alleles:
        alleles = hla_manager.suggest_alleles_for_population('European')[:3]
    
    # 提取肽段
    peptides = load_one(spectra_path, min_length=8, max_length=11)
    if not peptides:
        print(f"        Warning: No peptides found")
        return []
    
    # 分析修饰-锚位耦合（使用HLA管理器）
    coupling_records = tag_anchor_modifications(peptides, alleles, use_mhcflurry, hla_manager)
    
    # 添加样本信息和batch_id
    for record in coupling_records:
        record['sample'] = sample_name
        record['dataset'] = sample_record['dataset']
        record['batch_id'] = batch_id
        record['group'] = sample_type
    
    return coupling_records


def run_pipeline_batch(meta_file="all_meta.tsv", 
                      use_mhcflurry=False,
                      sample_limit=None,
                      output_dir="results_batch",
                      mode='batch'):  # 'batch' or 'traditional'
    """
    运行批次分析pipeline
    
    Args:
        meta_file: 元数据文件路径
        use_mhcflurry: 是否使用MHCflurry预测
        sample_limit: 限制处理的样本数量（用于测试）
        output_dir: 输出目录
        mode: 运行模式 - 'batch' (按批次分析) 或 'traditional' (传统合并分析)
    """
    print("=" * 60)
    print("修饰位点 × HLA 锚定位点耦合分析 Pipeline - 批次分析版本")
    print("=" * 60)
    print()
    
    start_time = time.time()
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    batch_dir = os.path.join(output_dir, "batch_results")
    os.makedirs(batch_dir, exist_ok=True)
    
    # 1. 加载元数据
    print("Step 1: Loading metadata...")
    metadata = load_metadata(meta_file)
    
    if sample_limit:
        metadata = metadata[:sample_limit]
        print(f"  Limited to first {sample_limit} samples for testing")
    
    print(f"  Found {len(metadata)} samples")
    
    # 统计样本类型和批次
    type_counts = {}
    batch_counts = {}
    for record in metadata:
        sample_type = record['type']
        batch_id = record.get('batch_id', record['dataset'])
        type_counts[sample_type] = type_counts.get(sample_type, 0) + 1
        batch_counts[batch_id] = batch_counts.get(batch_id, 0) + 1
    
    print("  Sample type distribution:")
    for sample_type, count in type_counts.items():
        print(f"    {sample_type}: {count}")
    print(f"  Total batches: {len(batch_counts)}")
    print()
    
    # 2. 处理所有样本
    print("Step 2: Processing samples...")
    
    # 初始化HLA管理器
    hla_manager = HLAManager()
    print(f"  Loaded HLA reference set: {len(hla_manager.get_all_alleles())} alleles")
    print()
    
    # 按批次收集数据
    batch_tables = {}  # key = batch_id → list[DataFrame]
    all_coupling_records = []
    
    # 按数据集分组处理
    datasets = {}
    for record in metadata:
        dataset = record['dataset']
        if dataset not in datasets:
            datasets[dataset] = []
        datasets[dataset].append(record)
    
    sample_counter = 0
    for dataset_idx, (dataset_name, samples) in enumerate(datasets.items()):
        print(f"  Dataset {dataset_idx+1}/{len(datasets)}: {dataset_name} ({len(samples)} samples)")
        
        for sample_idx, sample_record in enumerate(samples):
            sample_counter += 1
            batch_id = sample_record.get('batch_id', dataset_name)
            print(f"    Sample {sample_idx+1}/{len(samples)} ({sample_counter}/{len(metadata)}): {sample_record['sample']} ({sample_record['type']})")
            
            coupling_records = process_sample(sample_record, use_mhcflurry, hla_manager)
            
            if coupling_records:
                all_coupling_records.extend(coupling_records)
                
                # 转换为DataFrame并按批次存储
                df_records = pd.DataFrame(coupling_records)
                batch_tables.setdefault(batch_id, []).append(df_records)
                
                print(f"      -> Generated {len(coupling_records)} coupling records")
            else:
                print(f"      -> No records generated")
        
        print(f"  Dataset {dataset_name} completed.")
        print()
    
    print(f"Total coupling records: {len(all_coupling_records)}")
    print()
    
    if mode == 'batch':
        # 3. 批次级别分析
        print("Step 3: Batch-level enrichment analysis...")
        
        per_batch_results = []
        per_batch_group_comparisons = []
        
        for batch_idx, (bid, tables) in enumerate(batch_tables.items()):
            print(f"  Batch {batch_idx+1}/{len(batch_tables)}: {bid}")
            
            # 合并该批次的所有数据
            df_batch = pd.concat(tables, ignore_index=True)
            
            # 分离肿瘤和正常样本
            df_tumor = df_batch[df_batch['group'] == 'Tumor']
            df_normal = df_batch[df_batch['group'] == 'Normal']
            
            print(f"    Tumor samples: {df_tumor['sample'].nunique()}, Normal samples: {df_normal['sample'].nunique()}")
            
            # 肿瘤内部富集分析
            if len(df_tumor) > 0:
                tumor_enrich = batch_enrich(df_tumor)
                tumor_enrich['batch_id'] = bid
                tumor_enrich['group'] = 'Tumor'
                per_batch_results.append(tumor_enrich)
                
                # 导出批次结果
                batch_file = os.path.join(batch_dir, f"batch_{bid}_tumor_enrichment.csv")
                tumor_enrich.to_csv(batch_file, index=False)
            
            # 正常样本内部富集分析（如果有足够样本）
            if len(df_normal) > 0:
                normal_enrich = batch_enrich(df_normal)
                normal_enrich['batch_id'] = bid
                normal_enrich['group'] = 'Normal'
                per_batch_results.append(normal_enrich)
                
                # 导出批次结果
                batch_file = os.path.join(batch_dir, f"batch_{bid}_normal_enrichment.csv")
                normal_enrich.to_csv(batch_file, index=False)
            
            # 组间比较（如果两组都有数据）
            if len(df_tumor) > 0 and len(df_normal) > 0:
                group_comp = batch_group_comparison(df_tumor, df_normal)
                group_comp['batch_id'] = bid
                per_batch_group_comparisons.append(group_comp)
        
        # 导出所有批次结果
        if per_batch_results:
            all_batch_df = pd.concat(per_batch_results, ignore_index=True)
            batch_results_file = os.path.join(output_dir, "per_batch_anchor_enrich.csv")
            export_batch_stats(per_batch_results, batch_results_file)
            
            # 4. Meta分析
            print("\nStep 4: Meta-analysis across batches...")
            
            # 分别对肿瘤和正常样本进行meta分析
            tumor_batch_data = all_batch_df[all_batch_df['group'] == 'Tumor']
            normal_batch_data = all_batch_df[all_batch_df['group'] == 'Normal']
            
            # 肿瘤meta分析
            if len(tumor_batch_data) > 0:
                print("  Performing meta-analysis for tumor samples...")
                tumor_meta_results = meta_analysis_by_mod(tumor_batch_data, method='random')
                tumor_meta_file = os.path.join(output_dir, "meta_anchor_enrich_tumor.csv")
                export_meta_results(tumor_meta_results, tumor_meta_file)
                
                # 生成森林图（前5个最显著的修饰）
                forest_dir = os.path.join(output_dir, "forest_plots")
                os.makedirs(forest_dir, exist_ok=True)
                
                for i, (_, row) in enumerate(tumor_meta_results.head(5).iterrows()):
                    if row['n_batches'] >= 2:  # 至少2个批次才能画森林图
                        forest_file = os.path.join(forest_dir, f"forest_{row['mod']}_tumor.png")
                        try:
                            forest_plot(tumor_meta_results, tumor_batch_data, row['mod'], forest_file)
                        except Exception as e:
                            print(f"    Warning: Could not create forest plot for {row['mod']}: {e}")
                
                # 生成汇总森林图
                summary_forest_file = os.path.join(output_dir, "summary_forest_plot_tumor.png")
                create_summary_forest_plot(tumor_meta_results, summary_forest_file)
            
            # 正常样本meta分析（如果有足够批次）
            if len(normal_batch_data) > 0 and normal_batch_data['batch_id'].nunique() >= 2:
                print("  Performing meta-analysis for normal samples...")
                normal_meta_results = meta_analysis_by_mod(normal_batch_data, method='random')
                normal_meta_file = os.path.join(output_dir, "meta_anchor_enrich_normal.csv")
                export_meta_results(normal_meta_results, normal_meta_file)
            
            # 组间比较的meta分析
            if per_batch_group_comparisons:
                print("  Performing meta-analysis for group comparisons...")
                group_comp_df = pd.concat(per_batch_group_comparisons, ignore_index=True)
                
                # 重命名列以适应meta分析函数
                group_comp_df['logOR'] = group_comp_df['group_logOR']
                group_comp_df['SE'] = group_comp_df['group_SE']
                
                group_meta_results = meta_analysis_by_mod(group_comp_df, method='random')
                group_meta_file = os.path.join(output_dir, "meta_tumor_vs_normal.csv")
                export_meta_results(group_meta_results, group_meta_file)
            
            # 5. 评估异质性
            print("\nStep 5: Assessing heterogeneity...")
            heterogeneity_summary = summarize_batch_heterogeneity(all_batch_df)
            
            # 导出异质性报告
            hetero_file = os.path.join(output_dir, "heterogeneity_report.txt")
            with open(hetero_file, 'w') as f:
                f.write("Batch Heterogeneity Assessment\n")
                f.write("=" * 60 + "\n\n")
                
                for mod, stats in sorted(heterogeneity_summary.items(), 
                                        key=lambda x: x[1]['I2'], reverse=True):
                    f.write(f"Modification: {mod}\n")
                    f.write(f"  Number of batches: {stats['n_batches']}\n")
                    f.write(f"  I² = {stats['I2']:.1f}%\n")
                    f.write(f"  τ² = {stats['tau2']:.4f}\n")
                    f.write(f"  Mean logOR = {stats['mean_logOR']:.3f}\n")
                    f.write(f"  Range logOR = [{stats['range_logOR'][0]:.3f}, {stats['range_logOR'][1]:.3f}]\n")
                    f.write("\n")
            
            print(f"  Heterogeneity report saved to: {hetero_file}")
    
    else:  # Traditional mode
        # 传统分析模式（所有数据合并）
        print("Step 3: Traditional analysis (all data combined)...")
        
        # 导出原始耦合记录
        all_df = pd.DataFrame(all_coupling_records)
        tumor_df = all_df[all_df['group'] == 'Tumor']
        normal_df = all_df[all_df['group'] == 'Normal']
        
        tumor_file = os.path.join(output_dir, "tumor_anchor_mod.tsv")
        normal_file = os.path.join(output_dir, "normal_anchor_mod.tsv")
        all_file = os.path.join(output_dir, "all_anchor_mod.tsv")
        
        if len(tumor_df) > 0:
            export_coupling_records(tumor_df.to_dict('records'), tumor_file)
        if len(normal_df) > 0:
            export_coupling_records(normal_df.to_dict('records'), normal_file)
        if len(all_df) > 0:
            export_coupling_records(all_df.to_dict('records'), all_file)
        
        # 富集分析
        if len(all_coupling_records) > 0:
            enrichment_results = enrichment_analysis(all_coupling_records)
            enrichment_file = os.path.join(output_dir, "anchor_mod_enrichment.tsv")
            export_enrichment_results(enrichment_results, enrichment_file)
        
        # Tumor vs Normal 比较
        if len(tumor_df) > 0 and len(normal_df) > 0:
            comparison_results = compare_groups(tumor_df.to_dict('records'), 
                                              normal_df.to_dict('records'))
            comparison_file = os.path.join(output_dir, "anchor_mod_group_diff.csv")
            export_group_comparison_results(comparison_results, comparison_file)
    
    # 6. 生成总结报告
    print("\nStep 6: Generating summary report...")
    
    summary_file = os.path.join(output_dir, "pipeline_summary.txt")
    with open(summary_file, 'w', encoding='utf-8') as f:
        f.write("修饰位点 × HLA 锚定位点耦合分析 Pipeline 总结报告\n")
        f.write("=" * 60 + "\n\n")
        
        f.write(f"分析模式: {'批次分析 (Batch Analysis)' if mode == 'batch' else '传统分析 (Traditional Analysis)'}\n")
        f.write(f"处理时间: {time.time() - start_time:.2f} 秒\n")
        f.write(f"处理样本数: {len(metadata)}\n")
        f.write(f"总耦合记录: {len(all_coupling_records)}\n")
        f.write(f"批次数量: {len(batch_counts)}\n\n")
        
        f.write("样本类型分布:\n")
        for sample_type, count in type_counts.items():
            f.write(f"  {sample_type}: {count}\n")
        f.write("\n")
        
        f.write("批次分布:\n")
        for batch_id, count in sorted(batch_counts.items()):
            f.write(f"  {batch_id}: {count} samples\n")
        f.write("\n")
        
        f.write("输出文件:\n")
        for root, dirs, files in os.walk(output_dir):
            level = root.replace(output_dir, '').count(os.sep)
            indent = ' ' * 2 * level
            folder_name = os.path.basename(root)
            if folder_name:
                f.write(f"{indent}{folder_name}/\n")
            sub_indent = ' ' * 2 * (level + 1)
            for file in sorted(files):
                if file.endswith(('.tsv', '.csv', '.txt', '.png')):
                    f.write(f"{sub_indent}{file}\n")
    
    print(f"Summary report saved to: {summary_file}")
    
    # 完成
    end_time = time.time()
    print("\n" + "=" * 60)
    print(f"Pipeline completed successfully in {end_time - start_time:.2f} seconds")
    print(f"Results saved to: {output_dir}/")
    print("=" * 60)


def main():
    """主函数"""
    import argparse
    
    parser = argparse.ArgumentParser(description='修饰位点 × HLA 锚定位点耦合分析 Pipeline - 批次分析版本')
    parser.add_argument('--meta', default='all_meta.tsv', help='元数据文件路径 (默认: all_meta.tsv)')
    parser.add_argument('--output', default='results_batch', help='输出目录 (默认: results_batch)')
    parser.add_argument('--mhcflurry', action='store_true', help='使用MHCflurry预测 (需要安装mhcflurry)')
    parser.add_argument('--limit', type=int, help='限制处理的样本数量 (用于测试)')
    parser.add_argument('--test', action='store_true', help='测试模式 (只处理前5个样本)')
    parser.add_argument('--mode', choices=['batch', 'traditional'], default='batch', 
                       help='分析模式: batch (按批次) 或 traditional (传统合并)')
    
    args = parser.parse_args()
    
    # 测试模式
    if args.test:
        args.limit = 5
        print("Running in test mode (first 5 samples only)")
    
    # 运行pipeline
    run_pipeline_batch(
        meta_file=args.meta,
        use_mhcflurry=args.mhcflurry,
        sample_limit=args.limit,
        output_dir=args.output,
        mode=args.mode
    )


if __name__ == "__main__":
    main()