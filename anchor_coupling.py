#!/usr/bin/env python3
"""
修饰-锚位耦合分析模块
分析修饰位点与HLA锚位的耦合关系
"""

from extract_peptides import load_one, extract_mod_type
from predict_binding import parse_allele_string, batch_predict_binding

def tag_anchor_modifications(peptides, alleles, use_mhcflurry=True):
    """
    为肽段标注修饰与锚位的耦合关系
    
    Args:
        peptides: 肽段列表 [{'Sequence': str, 'mod_list': list}, ...]
        alleles: HLA等位基因列表
        use_mhcflurry: 是否使用MHCflurry预测
    
    Returns:
        list: 修饰-锚位耦合记录
    """
    # 首先进行HLA结合预测
    binding_results = batch_predict_binding(peptides, alleles, use_mhcflurry)
    
    coupling_records = []
    
    for i, result in enumerate(binding_results):
        seq = result['sequence']
        allele = result['best_allele'] 
        score = result['binding_score']
        anchors = result['anchor_positions']
        modifications = result['modifications']
        
        # 如果肽段有修饰，分析每个修饰与锚位的关系
        if modifications:
            for pos, mod in modifications:
                mod_type = extract_mod_type(mod)
                
                # 判断修饰是否在锚位
                is_anchor = pos in anchors
                anchor_tag = "anchor" if is_anchor else "non_anchor"
                
                coupling_records.append({
                    "seq": seq,
                    "allele": allele,
                    "score": score,
                    "mod": mod_type,
                    "mod_position": pos,
                    "mod_full": mod,
                    "anchor_tag": anchor_tag,
                    "anchor_positions": sorted(list(anchors)),
                    "sequence_length": len(seq)
                })
        else:
            # 无修饰的肽段，记录为无修饰记录
            coupling_records.append({
                "seq": seq,
                "allele": allele,
                "score": score,
                "mod": "Unmodified",
                "mod_position": None,
                "mod_full": None,
                "anchor_tag": "non_anchor",  # 无修饰时标记为非锚位
                "anchor_positions": sorted(list(anchors)),
                "sequence_length": len(seq)
            })
    
    return coupling_records

def analyze_anchor_enrichment(coupling_records):
    """
    分析修饰在锚位的富集情况
    
    Args:
        coupling_records: 耦合记录列表
    
    Returns:
        dict: 富集分析结果
    """
    # 统计每种修饰在锚位和非锚位的分布
    mod_counts = {}
    
    for record in coupling_records:
        mod = record['mod']
        anchor_tag = record['anchor_tag']
        
        if mod not in mod_counts:
            mod_counts[mod] = {"anchor": 0, "non_anchor": 0}
        
        mod_counts[mod][anchor_tag] += 1
    
    # 计算富集比例
    enrichment_results = {}
    
    for mod, counts in mod_counts.items():
        anchor_count = counts['anchor']
        non_anchor_count = counts['non_anchor']
        total = anchor_count + non_anchor_count
        
        if total > 0:
            anchor_percentage = (anchor_count / total) * 100
            enrichment_ratio = anchor_count / max(non_anchor_count, 1)  # 避免除零
            
            enrichment_results[mod] = {
                'anchor_count': anchor_count,
                'non_anchor_count': non_anchor_count,
                'total_count': total,
                'anchor_percentage': anchor_percentage,
                'enrichment_ratio': enrichment_ratio
            }
    
    return enrichment_results

def filter_by_modification_type(coupling_records, target_mods):
    """
    根据修饰类型过滤耦合记录
    
    Args:
        coupling_records: 耦合记录列表
        target_mods: 目标修饰类型列表
    
    Returns:
        list: 过滤后的记录
    """
    return [record for record in coupling_records if record['mod'] in target_mods]

def get_position_specific_analysis(coupling_records):
    """
    获取位置特异性分析
    
    Args:
        coupling_records: 耦合记录列表
    
    Returns:
        dict: 位置特异性分析结果
    """
    position_analysis = {}
    
    for record in coupling_records:
        if record['mod'] == 'Unmodified':
            continue
            
        mod = record['mod']
        position = record['mod_position']
        seq_len = record['sequence_length']
        
        # 计算相对位置（N-terminal, Middle, C-terminal）
        if position <= 2:
            relative_pos = "N-terminal"
        elif position >= seq_len - 1:
            relative_pos = "C-terminal"
        else:
            relative_pos = "Middle"
        
        key = f"{mod}_{relative_pos}"
        
        if key not in position_analysis:
            position_analysis[key] = {"anchor": 0, "non_anchor": 0}
        
        position_analysis[key][record['anchor_tag']] += 1
    
    return position_analysis

def summarize_by_allele(coupling_records):
    """
    按HLA等位基因汇总分析结果
    
    Args:
        coupling_records: 耦合记录列表
    
    Returns:
        dict: 按等位基因分组的结果
    """
    allele_summary = {}
    
    for record in coupling_records:
        allele = record['allele']
        mod = record['mod']
        anchor_tag = record['anchor_tag']
        
        if allele not in allele_summary:
            allele_summary[allele] = {}
        
        if mod not in allele_summary[allele]:
            allele_summary[allele][mod] = {"anchor": 0, "non_anchor": 0}
        
        allele_summary[allele][mod][anchor_tag] += 1
    
    return allele_summary

def export_coupling_records(coupling_records, output_file):
    """
    导出耦合记录到TSV文件
    
    Args:
        coupling_records: 耦合记录列表
        output_file: 输出文件路径
    """
    import csv
    
    if not coupling_records:
        print("No coupling records to export.")
        return
    
    fieldnames = [
        'seq', 'allele', 'score', 'mod', 'mod_position', 'mod_full',
        'anchor_tag', 'anchor_positions', 'sequence_length', 'sample', 'dataset', 'group'
    ]
    
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        
        for record in coupling_records:
            # 将anchor_positions列表转换为字符串
            record_copy = record.copy()
            record_copy['anchor_positions'] = ','.join(map(str, record['anchor_positions']))
            writer.writerow(record_copy)
    
    print(f"Coupling records exported to: {output_file}")

def main():
    """测试函数"""
    # 测试数据
    test_peptides = [
        {'Sequence': 'KFKESFAEM', 'mod_list': [(9, 'Oxidation[M]')]},
        {'Sequence': 'RILEMNDKYVK', 'mod_list': [(5, 'Oxidation[M]')]},
        {'Sequence': 'ADMAHISGL', 'mod_list': [(3, 'Oxidation[M]')]},
        {'Sequence': 'KLNPYAKTMR', 'mod_list': [(9, 'Oxidation[M]')]},
        {'Sequence': 'KLNPYAKTMR', 'mod_list': [(5, 'Phospho[Y]'), (9, 'Oxidation[M]')]},
        {'Sequence': 'ADMAHISGL', 'mod_list': []},  # 无修饰
    ]
    
    alleles = ['A*02:01', 'B*07:02', 'C*07:02']
    
    print("Testing anchor-modification coupling analysis...")
    
    # 分析修饰-锚位耦合
    coupling_records = tag_anchor_modifications(test_peptides, alleles, use_mhcflurry=False)
    
    print(f"\nGenerated {len(coupling_records)} coupling records:")
    for i, record in enumerate(coupling_records[:5]):  # 显示前5个
        print(f"{i+1}. {record['seq']}")
        print(f"   Modification: {record['mod']} at position {record['mod_position']}")
        print(f"   Anchor tag: {record['anchor_tag']}")
        print(f"   Anchor positions: {record['anchor_positions']}")
        print(f"   Score: {record['score']:.3f}")
        print()
    
    # 富集分析
    enrichment = analyze_anchor_enrichment(coupling_records)
    print("Enrichment analysis:")
    for mod, stats in enrichment.items():
        print(f"{mod}: {stats['anchor_count']} anchor, {stats['non_anchor_count']} non-anchor")
        print(f"  Anchor percentage: {stats['anchor_percentage']:.1f}%")
        print(f"  Enrichment ratio: {stats['enrichment_ratio']:.2f}")
        print()
    
    # 位置特异性分析
    position_analysis = get_position_specific_analysis(coupling_records)
    print("Position-specific analysis:")
    for pos_key, counts in position_analysis.items():
        print(f"{pos_key}: anchor={counts['anchor']}, non_anchor={counts['non_anchor']}")
    
    # 导出记录
    output_file = "test_coupling_records.tsv"
    export_coupling_records(coupling_records, output_file)

if __name__ == "__main__":
    main()