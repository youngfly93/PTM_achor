#!/usr/bin/env python3
"""
肽段提取模块 - 从spectra文件提取8-11mer肽段与修饰位点信息
"""

import re
import csv
from pathlib import Path

# 氨基酸正则表达式
AA = re.compile(r"[A-Z]")

def parse_mod(s):
    """
    解析修饰信息
    输入格式: '2,Oxidation[M];5,Phospho[S]' 
    输出: [(2,'Oxidation[M]'), (5,'Phospho[S]')]
    """
    if not s or s.strip() == '' or s.strip() == 'nan':
        return []
    
    modifications = []
    # 分号分隔多个修饰
    for mod_str in s.split(';'):
        if ',' in mod_str:
            parts = mod_str.split(',', 1)  # 只分割第一个逗号
            try:
                pos = int(parts[0])
                mod_type = parts[1]
                modifications.append((pos, mod_type))
            except ValueError:
                # 忽略格式错误的修饰
                continue
    
    return modifications

def extract_mod_type(mod_string):
    """
    从修饰字符串中提取修饰类型
    例如: 'Oxidation[M]' -> 'Oxidation'
         'Phospho[S]' -> 'Phospho'
    """
    if '[' in mod_string:
        return mod_string.split('[')[0]
    return mod_string

def load_one(spectra_path, min_length=8, max_length=11):
    """
    从单个spectra文件加载肽段数据
    
    Args:
        spectra_path: spectra文件路径
        min_length: 最小肽段长度
        max_length: 最大肽段长度
    
    Returns:
        list: 包含去重后的肽段信息 [{'Sequence': str, 'mod_list': list}, ...]
    """
    peptides = []
    seen_peptides = set()  # 用于去重
    
    try:
        with open(spectra_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            
            for row in reader:
                sequence = row.get('Sequence', '').strip()
                modification = row.get('Modification', '').strip()
                
                # 检查肽段长度
                if len(sequence) < min_length or len(sequence) > max_length:
                    continue
                
                # 检查是否为有效氨基酸序列
                if not sequence or not all(c in 'ACDEFGHIKLMNPQRSTVWY' for c in sequence):
                    continue
                
                # 解析修饰信息
                mod_list = parse_mod(modification)
                
                # 创建唯一标识符用于去重
                mod_str = ';'.join(f"{pos},{mod}" for pos, mod in mod_list)
                peptide_id = f"{sequence}:{mod_str}"
                
                if peptide_id not in seen_peptides:
                    seen_peptides.add(peptide_id)
                    peptides.append({
                        'Sequence': sequence,
                        'mod_list': mod_list
                    })
    
    except Exception as e:
        print(f"Error reading {spectra_path}: {e}")
        return []
    
    return peptides

def load_multiple(spectra_paths, min_length=8, max_length=11):
    """
    从多个spectra文件加载肽段数据
    
    Args:
        spectra_paths: spectra文件路径列表
        min_length: 最小肽段长度
        max_length: 最大肽段长度
    
    Returns:
        list: 合并去重后的肽段信息
    """
    all_peptides = []
    seen_peptides = set()  # 全局去重
    
    for path in spectra_paths:
        print(f"Processing: {Path(path).name}")
        peptides = load_one(path, min_length, max_length)
        
        for peptide in peptides:
            # 全局去重
            mod_str = ';'.join(f"{pos},{mod}" for pos, mod in peptide['mod_list'])
            peptide_id = f"{peptide['Sequence']}:{mod_str}"
            
            if peptide_id not in seen_peptides:
                seen_peptides.add(peptide_id)
                all_peptides.append(peptide)
    
    return all_peptides

def filter_by_modifications(peptides, target_modifications=None):
    """
    根据修饰类型过滤肽段
    
    Args:
        peptides: 肽段列表
        target_modifications: 目标修饰类型列表，如['Phospho', 'Oxidation']
                             如果为None，返回所有肽段
    
    Returns:
        list: 过滤后的肽段列表
    """
    if target_modifications is None:
        return peptides
    
    filtered = []
    for peptide in peptides:
        for pos, mod in peptide['mod_list']:
            mod_type = extract_mod_type(mod)
            if mod_type in target_modifications:
                filtered.append(peptide)
                break  # 找到一个匹配的修饰就足够了
    
    return filtered

def get_modification_stats(peptides):
    """
    获取修饰统计信息
    
    Args:
        peptides: 肽段列表
    
    Returns:
        dict: 修饰统计信息
    """
    stats = {
        'total_peptides': len(peptides),
        'modified_peptides': 0,
        'unmodified_peptides': 0,
        'modification_counts': {},
        'length_distribution': {}
    }
    
    for peptide in peptides:
        seq_len = len(peptide['Sequence'])
        stats['length_distribution'][seq_len] = stats['length_distribution'].get(seq_len, 0) + 1
        
        if peptide['mod_list']:
            stats['modified_peptides'] += 1
            for pos, mod in peptide['mod_list']:
                mod_type = extract_mod_type(mod)
                stats['modification_counts'][mod_type] = stats['modification_counts'].get(mod_type, 0) + 1
        else:
            stats['unmodified_peptides'] += 1
    
    return stats

def main():
    """测试函数"""
    # 测试单个文件
    test_file = "/mnt/f/work/yang_ylab/cancer_datasets_links/PXD000394_human/pFind-Filtered_res_openHLA_20120321_EXQ1_MiBa_SA_HCC1143_1_HCDFT.mgf_respFind.spectra"
    
    print("Testing peptide extraction...")
    peptides = load_one(test_file)
    
    print(f"Found {len(peptides)} unique peptides (8-11 mer)")
    
    # 显示前5个肽段
    print("\nFirst 5 peptides:")
    for i, peptide in enumerate(peptides[:5]):
        print(f"{i+1}. {peptide['Sequence']} (length: {len(peptide['Sequence'])})")
        if peptide['mod_list']:
            for pos, mod in peptide['mod_list']:
                print(f"   Modification at position {pos}: {mod}")
        else:
            print("   No modifications")
    
    # 统计信息
    stats = get_modification_stats(peptides)
    print(f"\nStatistics:")
    print(f"Total peptides: {stats['total_peptides']}")
    print(f"Modified: {stats['modified_peptides']}")
    print(f"Unmodified: {stats['unmodified_peptides']}")
    print(f"Length distribution: {stats['length_distribution']}")
    print(f"Modification counts: {stats['modification_counts']}")

if __name__ == "__main__":
    main()