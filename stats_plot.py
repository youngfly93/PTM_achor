#!/usr/bin/env python3
"""
统计检验和可视化模块
"""

import math
import csv
from collections import defaultdict

def fisher_exact_test(a, b, c, d):
    """
    计算Fisher精确检验的p值
    
    参数：
    a, b, c, d: 2x2列联表的四个值
    [[a, b],
     [c, d]]
    
    返回：
    p_value: 双尾检验的p值
    """
    
    def factorial(n):
        if n <= 1:
            return 1
        result = 1
        for i in range(2, n + 1):
            result *= i
        return result
    
    def combinations(n, k):
        if k > n or k < 0:
            return 0
        if k == 0 or k == n:
            return 1
        return factorial(n) // (factorial(k) * factorial(n - k))
    
    # 计算边际总数
    row1_total = a + b
    row2_total = c + d
    col1_total = a + c
    col2_total = b + d
    total = a + b + c + d
    
    if total == 0:
        return 1.0
    
    # 计算观测值的概率
    observed_prob = (combinations(row1_total, a) * combinations(row2_total, c)) / combinations(total, col1_total)
    
    # 计算所有可能情况的概率，找出极端情况
    p_value = 0.0
    
    # 确定a的可能范围
    min_a = max(0, col1_total - row2_total)
    max_a = min(row1_total, col1_total)
    
    for a_test in range(min_a, max_a + 1):
        b_test = row1_total - a_test
        c_test = col1_total - a_test
        d_test = row2_total - c_test
        
        if b_test >= 0 and c_test >= 0 and d_test >= 0:
            prob = (combinations(row1_total, a_test) * combinations(row2_total, c_test)) / combinations(total, col1_total)
            
            # 双尾检验：累加概率小于等于观测概率的所有情况
            if prob <= observed_prob + 1e-10:  # 加小量避免浮点误差
                p_value += prob
    
    return min(p_value, 1.0)

def enrichment_analysis(coupling_records):
    """
    修饰在锚位的富集分析（Fisher精确检验）
    
    Args:
        coupling_records: 耦合记录列表
    
    Returns:
        list: 富集分析结果，按p值排序
    """
    # 统计每种修饰在锚位和非锚位的分布
    mod_counts = defaultdict(lambda: {"anchor": 0, "non_anchor": 0})
    
    for record in coupling_records:
        mod = record['mod']
        anchor_tag = record['anchor_tag']
        mod_counts[mod][anchor_tag] += 1
    
    results = []
    
    for mod, counts in mod_counts.items():
        anchor = counts['anchor']
        non_anchor = counts['non_anchor']
        
        # 计算其他修饰的统计
        other_anchor = sum(c['anchor'] for m, c in mod_counts.items() if m != mod)
        other_non_anchor = sum(c['non_anchor'] for m, c in mod_counts.items() if m != mod)
        
        # Fisher精确检验
        # [[当前修饰_锚位, 当前修饰_非锚位],
        #  [其他修饰_锚位, 其他修饰_非锚位]]
        p_value = fisher_exact_test(anchor, non_anchor, other_anchor, other_non_anchor)
        
        # 计算富集倍数
        if non_anchor > 0 and other_non_anchor > 0 and other_anchor > 0:
            enrichment_fold = (anchor / non_anchor) / (other_anchor / other_non_anchor)
        elif anchor > 0 and non_anchor == 0:
            enrichment_fold = float('inf')
        elif anchor == 0:
            enrichment_fold = 0
        else:
            enrichment_fold = 1
        
        results.append({
            'mod': mod,
            'anchor': anchor,
            'non_anchor': non_anchor,
            'total': anchor + non_anchor,
            'anchor_percentage': (anchor / (anchor + non_anchor)) * 100 if (anchor + non_anchor) > 0 else 0,
            'enrichment_fold': enrichment_fold,
            'p_value': p_value
        })
    
    # 按p值排序
    return sorted(results, key=lambda x: x['p_value'])

def create_simple_violin_data(coupling_records, mod_type):
    """
    为指定修饰类型创建violin plot数据（简化版本）
    
    Args:
        coupling_records: 耦合记录列表
        mod_type: 修饰类型
    
    Returns:
        dict: 包含anchor和non_anchor组的结合评分数据
    """
    data = {'anchor': [], 'non_anchor': []}
    
    for record in coupling_records:
        if record['mod'] == mod_type:
            group = record['anchor_tag']
            score = record['score']
            data[group].append(score)
    
    return data

def simple_violin_stats(data):
    """
    计算violin plot的简单统计信息
    
    Args:
        data: 数据列表
    
    Returns:
        dict: 统计信息
    """
    if not data:
        return {'count': 0, 'mean': 0, 'median': 0, 'std': 0}
    
    data_sorted = sorted(data)
    n = len(data)
    mean_val = sum(data) / n
    
    # 中位数
    if n % 2 == 0:
        median_val = (data_sorted[n//2 - 1] + data_sorted[n//2]) / 2
    else:
        median_val = data_sorted[n//2]
    
    # 标准差
    variance = sum((x - mean_val) ** 2 for x in data) / n
    std_val = math.sqrt(variance)
    
    return {
        'count': n,
        'mean': mean_val,
        'median': median_val,
        'std': std_val,
        'min': min(data),
        'max': max(data)
    }

def compare_anchor_vs_non_anchor(coupling_records, mod_type):
    """
    比较特定修饰在锚位vs非锚位的结合评分
    
    Args:
        coupling_records: 耦合记录列表
        mod_type: 修饰类型
    
    Returns:
        dict: 比较结果
    """
    violin_data = create_simple_violin_data(coupling_records, mod_type)
    
    anchor_stats = simple_violin_stats(violin_data['anchor'])
    non_anchor_stats = simple_violin_stats(violin_data['non_anchor'])
    
    # 简单的t检验（假设方差相等）
    def welch_t_test(data1, data2):
        if len(data1) < 2 or len(data2) < 2:
            return 1.0  # 样本太小，返回不显著
        
        mean1 = sum(data1) / len(data1)
        mean2 = sum(data2) / len(data2)
        
        var1 = sum((x - mean1) ** 2 for x in data1) / (len(data1) - 1)
        var2 = sum((x - mean2) ** 2 for x in data2) / (len(data2) - 1)
        
        if var1 == 0 and var2 == 0:
            return 1.0 if mean1 == mean2 else 0.0
        
        se = math.sqrt(var1 / len(data1) + var2 / len(data2))
        if se == 0:
            return 1.0
        
        t_stat = abs(mean1 - mean2) / se
        
        # 简化的p值估算（正态近似）
        # 这里使用简化公式，实际应用中建议使用scipy.stats
        p_approx = 2 * (1 - min(0.5 + t_stat / 6, 0.95))  # 非常粗糙的近似
        
        return p_approx
    
    p_value = welch_t_test(violin_data['anchor'], violin_data['non_anchor'])
    
    return {
        'modification': mod_type,
        'anchor_stats': anchor_stats,
        'non_anchor_stats': non_anchor_stats,
        'p_value': p_value,
        'violin_data': violin_data
    }

def generate_text_violin_plot(violin_data, mod_type, output_file=None):
    """
    生成文本版本的violin plot
    
    Args:
        violin_data: violin plot数据
        mod_type: 修饰类型
        output_file: 输出文件路径（可选）
    """
    anchor_data = violin_data['anchor']
    non_anchor_data = violin_data['non_anchor']
    
    anchor_stats = simple_violin_stats(anchor_data)
    non_anchor_stats = simple_violin_stats(non_anchor_data)
    
    plot_text = []
    plot_text.append(f"Violin Plot for {mod_type} Modification")
    plot_text.append("=" * 50)
    plot_text.append("")
    
    plot_text.append("Anchor Position:")
    plot_text.append(f"  Count: {anchor_stats['count']}")
    plot_text.append(f"  Mean: {anchor_stats['mean']:.3f}")
    plot_text.append(f"  Median: {anchor_stats['median']:.3f}")
    plot_text.append(f"  Std: {anchor_stats['std']:.3f}")
    plot_text.append(f"  Range: [{anchor_stats['min']:.3f}, {anchor_stats['max']:.3f}]")
    plot_text.append("")
    
    plot_text.append("Non-Anchor Position:")
    plot_text.append(f"  Count: {non_anchor_stats['count']}")
    plot_text.append(f"  Mean: {non_anchor_stats['mean']:.3f}")
    plot_text.append(f"  Median: {non_anchor_stats['median']:.3f}")
    plot_text.append(f"  Std: {non_anchor_stats['std']:.3f}")
    plot_text.append(f"  Range: [{non_anchor_stats['min']:.3f}, {non_anchor_stats['max']:.3f}]")
    plot_text.append("")
    
    # 简单的ASCII图形展示
    if anchor_data and non_anchor_data:
        plot_text.append("Distribution (simplified):")
        plot_text.append("Anchor:     " + "*" * min(anchor_stats['count'], 50))
        plot_text.append("Non-anchor: " + "*" * min(non_anchor_stats['count'], 50))
    
    result_text = "\n".join(plot_text)
    
    if output_file:
        with open(output_file, 'w') as f:
            f.write(result_text)
        print(f"Text violin plot saved to: {output_file}")
    else:
        print(result_text)

def export_enrichment_results(enrichment_results, output_file):
    """
    导出富集分析结果
    
    Args:
        enrichment_results: 富集分析结果列表
        output_file: 输出文件路径
    """
    fieldnames = [
        'mod', 'anchor', 'non_anchor', 'total', 'anchor_percentage',
        'enrichment_fold', 'p_value'
    ]
    
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(enrichment_results)
    
    print(f"Enrichment results exported to: {output_file}")

def main():
    """测试函数"""
    # 读取测试数据
    test_file = "test_coupling_records.tsv"
    
    try:
        coupling_records = []
        with open(test_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                # 转换数值字段
                row['score'] = float(row['score'])
                mod_pos = row['mod_position'].strip()
                row['mod_position'] = int(mod_pos) if mod_pos and mod_pos != 'None' else None
                row['sequence_length'] = int(row['sequence_length'])
                row['anchor_positions'] = [int(x) for x in row['anchor_positions'].split(',') if x.strip()]
                coupling_records.append(row)
        
        print(f"Loaded {len(coupling_records)} coupling records for testing")
        
        # 富集分析
        enrichment_results = enrichment_analysis(coupling_records)
        print("\nEnrichment Analysis Results:")
        print("-" * 70)
        for result in enrichment_results:
            print(f"Modification: {result['mod']}")
            print(f"  Anchor: {result['anchor']}, Non-anchor: {result['non_anchor']}")
            print(f"  Anchor percentage: {result['anchor_percentage']:.1f}%")
            print(f"  Enrichment fold: {result['enrichment_fold']:.2f}")
            print(f"  P-value: {result['p_value']:.6f}")
            print()
        
        # 导出富集结果
        export_enrichment_results(enrichment_results, "test_enrichment_results.tsv")
        
        # Violin plot分析（如果有Oxidation修饰）
        oxidation_records = [r for r in coupling_records if r['mod'] == 'Oxidation']
        if oxidation_records:
            print("Generating violin plot for Oxidation modification...")
            comparison = compare_anchor_vs_non_anchor(coupling_records, 'Oxidation')
            print(f"Anchor vs Non-anchor comparison p-value: {comparison['p_value']:.6f}")
            
            generate_text_violin_plot(comparison['violin_data'], 'Oxidation', 'oxidation_violin_plot.txt')
    
    except FileNotFoundError:
        print(f"Test file {test_file} not found. Please run anchor_coupling.py first.")

if __name__ == "__main__":
    main()