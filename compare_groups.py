#!/usr/bin/env python3
"""
Tumor vs Normal 差异分析模块
比较肿瘤样本和正常样本中修饰在锚位的分布差异
"""

import csv
from collections import defaultdict
from stats_plot import fisher_exact_test

def load_coupling_data(file_path):
    """
    从文件加载耦合数据
    
    Args:
        file_path: 耦合记录文件路径
    
    Returns:
        list: 耦合记录列表
    """
    coupling_records = []
    
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                # 转换数值字段
                row['score'] = float(row['score'])
                mod_pos = row['mod_position'].strip()
                row['mod_position'] = int(mod_pos) if mod_pos and mod_pos != 'None' else None
                row['sequence_length'] = int(row['sequence_length'])
                row['anchor_positions'] = [int(x) for x in row['anchor_positions'].split(',') if x.strip()]
                coupling_records.append(row)
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return []
    except Exception as e:
        print(f"Error loading coupling data: {e}")
        return []
    
    return coupling_records

def compare_groups(tumor_records, normal_records):
    """
    比较肿瘤组和正常组中修饰在锚位的分布差异
    
    Args:
        tumor_records: 肿瘤组耦合记录
        normal_records: 正常组耦合记录
    
    Returns:
        list: 组间比较结果，按p值排序
    """
    # 统计肿瘤组中每种修饰的分布
    tumor_counts = defaultdict(lambda: {"anchor": 0, "non_anchor": 0})
    for record in tumor_records:
        mod = record['mod']
        anchor_tag = record['anchor_tag']
        tumor_counts[mod][anchor_tag] += 1
    
    # 统计正常组中每种修饰的分布
    normal_counts = defaultdict(lambda: {"anchor": 0, "non_anchor": 0})
    for record in normal_records:
        mod = record['mod']
        anchor_tag = record['anchor_tag']
        normal_counts[mod][anchor_tag] += 1
    
    # 获取所有修饰类型
    all_mods = set(tumor_counts.keys()) | set(normal_counts.keys())
    
    results = []
    
    for mod in all_mods:
        t_anchor = tumor_counts[mod]["anchor"]
        t_non = tumor_counts[mod]["non_anchor"]
        n_anchor = normal_counts[mod]["anchor"]
        n_non = normal_counts[mod]["non_anchor"]
        
        # Fisher精确检验比较两组
        # [[tumor_anchor, tumor_non_anchor],
        #  [normal_anchor, normal_non_anchor]]
        p_value = fisher_exact_test(t_anchor, t_non, n_anchor, n_non)
        
        # 计算比例和富集比
        tumor_total = t_anchor + t_non
        normal_total = n_anchor + n_non
        
        tumor_anchor_rate = (t_anchor / tumor_total) * 100 if tumor_total > 0 else 0
        normal_anchor_rate = (n_anchor / normal_total) * 100 if normal_total > 0 else 0
        
        # 计算肿瘤vs正常的富集比
        if n_anchor > 0 and n_non > 0 and t_non > 0:
            enrichment_ratio = (t_anchor / t_non) / (n_anchor / n_non)
        elif t_anchor > 0 and n_anchor == 0:
            enrichment_ratio = float('inf')
        elif t_anchor == 0:
            enrichment_ratio = 0
        else:
            enrichment_ratio = 1
        
        results.append({
            'mod': mod,
            't_anchor': t_anchor,
            't_non': t_non,
            't_total': tumor_total,
            't_anchor_rate': tumor_anchor_rate,
            'n_anchor': n_anchor,
            'n_non': n_non,
            'n_total': normal_total,
            'n_anchor_rate': normal_anchor_rate,
            'enrichment_ratio': enrichment_ratio,
            'p_group': p_value
        })
    
    # 按p值排序
    return sorted(results, key=lambda x: x['p_group'])

def calculate_effect_size(t_anchor, t_non, n_anchor, n_non):
    """
    计算效应量（Odds Ratio）
    
    Args:
        t_anchor, t_non, n_anchor, n_non: 2x2列联表的值
    
    Returns:
        float: Odds Ratio
    """
    if t_non == 0 or n_anchor == 0:
        return float('inf') if t_anchor > 0 and n_non > 0 else 0
    
    if t_anchor == 0 or n_non == 0:
        return 0
    
    return (t_anchor * n_non) / (t_non * n_anchor)

def multiple_testing_correction(p_values, method='bonferroni'):
    """
    多重检验校正
    
    Args:
        p_values: p值列表
        method: 校正方法 ('bonferroni' 或 'fdr')
    
    Returns:
        list: 校正后的p值
    """
    n = len(p_values)
    
    if method == 'bonferroni':
        # Bonferroni校正
        corrected = [min(p * n, 1.0) for p in p_values]
    elif method == 'fdr':
        # Benjamini-Hochberg FDR校正（简化版）
        indexed_p = [(i, p) for i, p in enumerate(p_values)]
        indexed_p.sort(key=lambda x: x[1])  # 按p值排序
        
        corrected = [0] * n
        for rank, (original_index, p_val) in enumerate(indexed_p):
            corrected_p = p_val * n / (rank + 1)
            corrected[original_index] = min(corrected_p, 1.0)
        
        # 确保单调性
        for i in range(n - 2, -1, -1):
            if corrected[i] > corrected[i + 1]:
                corrected[i] = corrected[i + 1]
    else:
        corrected = p_values
    
    return corrected

def export_group_comparison_results(results, output_file):
    """
    导出组间比较结果
    
    Args:
        results: 比较结果列表
        output_file: 输出文件路径
    """
    fieldnames = [
        'mod', 't_anchor', 't_non', 't_total', 't_anchor_rate',
        'n_anchor', 'n_non', 'n_total', 'n_anchor_rate',
        'enrichment_ratio', 'odds_ratio', 'p_group', 'p_bonferroni', 'p_fdr'
    ]
    
    # 计算效应量和多重检验校正
    p_values = [r['p_group'] for r in results]
    p_bonferroni = multiple_testing_correction(p_values, 'bonferroni')
    p_fdr = multiple_testing_correction(p_values, 'fdr')
    
    for i, result in enumerate(results):
        result['odds_ratio'] = calculate_effect_size(
            result['t_anchor'], result['t_non'],
            result['n_anchor'], result['n_non']
        )
        result['p_bonferroni'] = p_bonferroni[i]
        result['p_fdr'] = p_fdr[i]
    
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(results)
    
    print(f"Group comparison results exported to: {output_file}")

def summarize_group_differences(results, significance_threshold=0.05):
    """
    汇总组间差异分析结果
    
    Args:
        results: 比较结果列表
        significance_threshold: 显著性阈值
    
    Returns:
        dict: 汇总信息
    """
    total_modifications = len(results)
    significant_raw = sum(1 for r in results if r['p_group'] < significance_threshold)
    significant_bonf = sum(1 for r in results if r.get('p_bonferroni', 1) < significance_threshold)
    significant_fdr = sum(1 for r in results if r.get('p_fdr', 1) < significance_threshold)
    
    # 找出最显著的差异
    most_significant = min(results, key=lambda x: x['p_group']) if results else None
    
    # 找出富集比最高的修饰
    enriched_mods = [r for r in results if r['enrichment_ratio'] > 1 and r['t_total'] > 0]
    depleted_mods = [r for r in results if r['enrichment_ratio'] < 1 and r['t_total'] > 0]
    
    return {
        'total_modifications': total_modifications,
        'significant_raw': significant_raw,
        'significant_bonferroni': significant_bonf,
        'significant_fdr': significant_fdr,
        'most_significant': most_significant,
        'enriched_in_tumor': len(enriched_mods),
        'depleted_in_tumor': len(depleted_mods)
    }

def create_comparison_report(results, summary, output_file):
    """
    创建比较分析报告
    
    Args:
        results: 比较结果
        summary: 汇总信息
        output_file: 输出文件路径
    """
    lines = []
    lines.append("Tumor vs Normal Comparison Analysis Report")
    lines.append("=" * 50)
    lines.append("")
    
    lines.append("Summary:")
    lines.append(f"  Total modifications analyzed: {summary['total_modifications']}")
    lines.append(f"  Significant differences (raw p < 0.05): {summary['significant_raw']}")
    lines.append(f"  Significant differences (Bonferroni corrected): {summary['significant_bonferroni']}")
    lines.append(f"  Significant differences (FDR corrected): {summary['significant_fdr']}")
    lines.append(f"  Enriched in tumor: {summary['enriched_in_tumor']}")
    lines.append(f"  Depleted in tumor: {summary['depleted_in_tumor']}")
    lines.append("")
    
    if summary['most_significant']:
        ms = summary['most_significant']
        lines.append("Most significant difference:")
        lines.append(f"  Modification: {ms['mod']}")
        lines.append(f"  P-value: {ms['p_group']:.6f}")
        lines.append(f"  Tumor anchor rate: {ms['t_anchor_rate']:.1f}%")
        lines.append(f"  Normal anchor rate: {ms['n_anchor_rate']:.1f}%")
        lines.append("")
    
    lines.append("Top 5 most significant modifications:")
    lines.append("-" * 40)
    for i, result in enumerate(results[:5]):
        lines.append(f"{i+1}. {result['mod']}")
        lines.append(f"   Tumor: {result['t_anchor']}/{result['t_total']} ({result['t_anchor_rate']:.1f}%) at anchor")
        lines.append(f"   Normal: {result['n_anchor']}/{result['n_total']} ({result['n_anchor_rate']:.1f}%) at anchor")
        lines.append(f"   P-value: {result['p_group']:.6f}")
        lines.append(f"   Enrichment ratio: {result['enrichment_ratio']:.2f}")
        lines.append("")
    
    report_text = "\n".join(lines)
    
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(report_text)
    
    print(f"Comparison report saved to: {output_file}")
    print(report_text)

def main():
    """测试函数"""
    # 创建模拟的肿瘤和正常样本数据用于测试
    print("Creating simulated tumor and normal samples for testing...")
    
    # 模拟肿瘤样本（更多Phospho修饰在锚位）
    tumor_records = [
        {'mod': 'Phospho', 'anchor_tag': 'anchor', 'score': 1.5},
        {'mod': 'Phospho', 'anchor_tag': 'anchor', 'score': 1.3},
        {'mod': 'Phospho', 'anchor_tag': 'non_anchor', 'score': 1.1},
        {'mod': 'Oxidation', 'anchor_tag': 'anchor', 'score': 0.9},
        {'mod': 'Oxidation', 'anchor_tag': 'non_anchor', 'score': 0.8},
        {'mod': 'Oxidation', 'anchor_tag': 'non_anchor', 'score': 0.7},
        {'mod': 'Acetyl', 'anchor_tag': 'anchor', 'score': 1.2},
        {'mod': 'Unmodified', 'anchor_tag': 'non_anchor', 'score': 0.5}
    ]
    
    # 模拟正常样本（Phospho修饰较少在锚位）
    normal_records = [
        {'mod': 'Phospho', 'anchor_tag': 'non_anchor', 'score': 1.0},
        {'mod': 'Phospho', 'anchor_tag': 'non_anchor', 'score': 0.9},
        {'mod': 'Oxidation', 'anchor_tag': 'anchor', 'score': 0.8},
        {'mod': 'Oxidation', 'anchor_tag': 'non_anchor', 'score': 0.7},
        {'mod': 'Acetyl', 'anchor_tag': 'non_anchor', 'score': 1.0},
        {'mod': 'Unmodified', 'anchor_tag': 'non_anchor', 'score': 0.4}
    ]
    
    print(f"Tumor samples: {len(tumor_records)}")
    print(f"Normal samples: {len(normal_records)}")
    
    # 进行组间比较
    comparison_results = compare_groups(tumor_records, normal_records)
    
    print("\nGroup Comparison Results:")
    print("-" * 60)
    for result in comparison_results:
        print(f"Modification: {result['mod']}")
        print(f"  Tumor: {result['t_anchor']}/{result['t_total']} ({result['t_anchor_rate']:.1f}%) at anchor")
        print(f"  Normal: {result['n_anchor']}/{result['n_total']} ({result['n_anchor_rate']:.1f}%) at anchor")
        print(f"  Enrichment ratio: {result['enrichment_ratio']:.2f}")
        print(f"  P-value: {result['p_group']:.6f}")
        print()
    
    # 导出结果
    export_group_comparison_results(comparison_results, "test_group_comparison.tsv")
    
    # 生成汇总报告
    summary = summarize_group_differences(comparison_results)
    create_comparison_report(comparison_results, summary, "test_comparison_report.txt")

if __name__ == "__main__":
    main()