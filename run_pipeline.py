#!/usr/bin/env python3
"""
修饰位点 × HLA 锚定位点耦合分析 Pipeline 主运行脚本

根据指导文档，执行以下步骤：
1. 读取 all_meta.tsv
2. 预测 8–11 mer 修饰肽 HLA 结合
3. 标记 P2/PΩ 锚位
4. 统计 Tumor vs Normal 中各修饰在锚位的富集（Fisher exact）
5. 生成 anchor_mod_group_diff.csv
6. 为 Phospho 绘制 binding‑score violin‑box 图
"""

import csv
import os
import sys
from pathlib import Path
import time

# 导入自定义模块
from extract_peptides import load_one, get_modification_stats
from predict_binding import parse_allele_string, batch_predict_binding
from anchor_coupling import tag_anchor_modifications, export_coupling_records
from stats_plot import enrichment_analysis, export_enrichment_results, compare_anchor_vs_non_anchor, generate_text_violin_plot
from compare_groups import compare_groups, export_group_comparison_results, summarize_group_differences, create_comparison_report
from hla_manager import HLAManager

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
    处理单个样本（集成HLA管理器）
    
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
    
    print(f"  Processing sample: {sample_name} ({sample_type})")
    print(f"    Spectra file: {Path(spectra_path).name}")
    
    # 检查文件是否存在
    if not os.path.exists(spectra_path):
        print(f"    Warning: Spectra file not found, skipping: {spectra_path}")
        return []
    
    # 解析HLA等位基因（使用HLA管理器）
    alleles = parse_allele_string(hla_alleles_str, hla_manager)
    if not alleles:
        print(f"    Warning: No valid HLA alleles found, using population defaults")
        alleles = hla_manager.suggest_alleles_for_population('European')[:3]
    
    # 提取肽段
    peptides = load_one(spectra_path, min_length=8, max_length=11)
    if not peptides:
        print(f"    Warning: No peptides found in {spectra_path}")
        return []
    
    print(f"    Found {len(peptides)} unique 8-11mer peptides")
    
    # 分析修饰-锚位耦合（使用HLA管理器）
    coupling_records = tag_anchor_modifications(peptides, alleles, use_mhcflurry, hla_manager)
    
    # 添加样本信息
    for record in coupling_records:
        record['sample'] = sample_name
        record['dataset'] = sample_record['dataset']
        record['group'] = sample_type
    
    # 显示修饰统计
    mod_stats = get_modification_stats(peptides)
    print(f"    Modified peptides: {mod_stats['modified_peptides']}")
    print(f"    Unmodified peptides: {mod_stats['unmodified_peptides']}")
    print(f"    Coupling records generated: {len(coupling_records)}")
    
    return coupling_records

def run_pipeline(meta_file="all_meta.tsv", 
                 use_mhcflurry=False,
                 sample_limit=None,
                 output_dir="results"):
    """
    运行完整的分析pipeline
    
    Args:
        meta_file: 元数据文件路径
        use_mhcflurry: 是否使用MHCflurry预测
        sample_limit: 限制处理的样本数量（用于测试）
        output_dir: 输出目录
    """
    print("=" * 60)
    print("修饰位点 × HLA 锚定位点耦合分析 Pipeline")
    print("=" * 60)
    print()
    
    start_time = time.time()
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. 加载元数据
    print("Step 1: Loading metadata...")
    metadata = load_metadata(meta_file)
    
    if sample_limit:
        metadata = metadata[:sample_limit]
        print(f"  Limited to first {sample_limit} samples for testing")
    
    print(f"  Found {len(metadata)} samples")
    
    # 统计样本类型
    type_counts = {}
    for record in metadata:
        sample_type = record['type']
        type_counts[sample_type] = type_counts.get(sample_type, 0) + 1
    
    print("  Sample type distribution:")
    for sample_type, count in type_counts.items():
        print(f"    {sample_type}: {count}")
    print()
    
    # 2. 处理所有样本
    print("Step 2: Processing samples...")
    
    # 初始化HLA管理器
    hla_manager = HLAManager()
    print(f"  Loaded HLA reference set: {len(hla_manager.get_all_alleles())} alleles")
    print(f"  Supported lengths: {hla_manager.get_supported_lengths()}")
    print()
    
    all_coupling_records = []
    tumor_records = []
    normal_records = []
    
    for i, sample_record in enumerate(metadata):
        print(f"  Processing sample {i+1}/{len(metadata)}:")
        
        coupling_records = process_sample(sample_record, use_mhcflurry, hla_manager)
        
        if coupling_records:
            all_coupling_records.extend(coupling_records)
            
            # 按组分类
            if sample_record['type'] == 'Tumor':
                tumor_records.extend(coupling_records)
            elif sample_record['type'] == 'Normal':
                normal_records.extend(coupling_records)
        
        print()
    
    print(f"Total coupling records: {len(all_coupling_records)}")
    print(f"Tumor records: {len(tumor_records)}")
    print(f"Normal records: {len(normal_records)}")
    print()
    
    # 3. 导出原始耦合记录
    print("Step 3: Exporting coupling records...")
    
    tumor_file = os.path.join(output_dir, "tumor_anchor_mod.tsv")
    normal_file = os.path.join(output_dir, "normal_anchor_mod.tsv")
    all_file = os.path.join(output_dir, "all_anchor_mod.tsv")
    
    if tumor_records:
        export_coupling_records(tumor_records, tumor_file)
    if normal_records:
        export_coupling_records(normal_records, normal_file)
    if all_coupling_records:
        export_coupling_records(all_coupling_records, all_file)
    
    # 4. 富集分析
    print("\nStep 4: Enrichment analysis...")
    
    if all_coupling_records:
        enrichment_results = enrichment_analysis(all_coupling_records)
        enrichment_file = os.path.join(output_dir, "anchor_mod_enrichment.tsv")
        export_enrichment_results(enrichment_results, enrichment_file)
        
        print("Top enriched modifications:")
        for result in enrichment_results[:5]:
            print(f"  {result['mod']}: {result['anchor_percentage']:.1f}% at anchor (p={result['p_value']:.6f})")
    
    # 5. Tumor vs Normal 比较
    print("\nStep 5: Tumor vs Normal comparison...")
    
    if tumor_records and normal_records:
        comparison_results = compare_groups(tumor_records, normal_records)
        comparison_file = os.path.join(output_dir, "anchor_mod_group_diff.csv")
        export_group_comparison_results(comparison_results, comparison_file)
        
        # 生成汇总报告
        summary = summarize_group_differences(comparison_results)
        report_file = os.path.join(output_dir, "tumor_vs_normal_report.txt")
        create_comparison_report(comparison_results, summary, report_file)
        
        print("Top differential modifications:")
        for result in comparison_results[:5]:
            print(f"  {result['mod']}: T={result['t_anchor_rate']:.1f}% vs N={result['n_anchor_rate']:.1f}% (p={result['p_group']:.6f})")
    else:
        print("  Insufficient tumor or normal samples for comparison")
    
    # 6. Phospho修饰的violin plot
    print("\nStep 6: Generating Phospho violin plot...")
    
    if all_coupling_records:
        phospho_comparison = compare_anchor_vs_non_anchor(all_coupling_records, 'Phospho')
        if phospho_comparison['violin_data']['anchor'] or phospho_comparison['violin_data']['non_anchor']:
            violin_file = os.path.join(output_dir, "Phospho_anchor_violin.txt")
            generate_text_violin_plot(phospho_comparison['violin_data'], 'Phospho', violin_file)
            print(f"  Phospho anchor vs non-anchor p-value: {phospho_comparison['p_value']:.6f}")
        else:
            print("  No Phospho modifications found for violin plot")
    
    # 7. 生成总结报告
    print("\nStep 7: Generating summary report...")
    
    summary_file = os.path.join(output_dir, "pipeline_summary.txt")
    with open(summary_file, 'w', encoding='utf-8') as f:
        f.write("修饰位点 × HLA 锚定位点耦合分析 Pipeline 总结报告\n")
        f.write("=" * 60 + "\n\n")
        
        f.write(f"处理时间: {time.time() - start_time:.2f} 秒\n")
        f.write(f"处理样本数: {len(metadata)}\n")
        f.write(f"总耦合记录: {len(all_coupling_records)}\n")
        f.write(f"肿瘤样本记录: {len(tumor_records)}\n")
        f.write(f"正常样本记录: {len(normal_records)}\n\n")
        
        f.write("样本类型分布:\n")
        for sample_type, count in type_counts.items():
            f.write(f"  {sample_type}: {count}\n")
        f.write("\n")
        
        f.write("输出文件:\n")
        for filename in os.listdir(output_dir):
            if filename.endswith(('.tsv', '.csv', '.txt')):
                f.write(f"  {filename}\n")
    
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
    
    parser = argparse.ArgumentParser(description='修饰位点 × HLA 锚定位点耦合分析 Pipeline')
    parser.add_argument('--meta', default='all_meta.tsv', help='元数据文件路径 (默认: all_meta.tsv)')
    parser.add_argument('--output', default='results', help='输出目录 (默认: results)')
    parser.add_argument('--mhcflurry', action='store_true', help='使用MHCflurry预测 (需要安装mhcflurry)')
    parser.add_argument('--limit', type=int, help='限制处理的样本数量 (用于测试)')
    parser.add_argument('--test', action='store_true', help='测试模式 (只处理前5个样本)')
    
    args = parser.parse_args()
    
    # 测试模式
    if args.test:
        args.limit = 5
        print("Running in test mode (first 5 samples only)")
    
    # 运行pipeline
    run_pipeline(
        meta_file=args.meta,
        use_mhcflurry=args.mhcflurry,
        sample_limit=args.limit,
        output_dir=args.output
    )

if __name__ == "__main__":
    main()