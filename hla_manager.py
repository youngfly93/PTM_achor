#!/usr/bin/env python3
"""
HLA管理模块 - 处理HLA参考位点文件和等位基因管理
"""

import csv
import re
from collections import defaultdict
from pathlib import Path

class HLAManager:
    """HLA等位基因和位点管理器"""
    
    def __init__(self, hla_ref_file="hla_ref_set.class_i.txt"):
        """
        初始化HLA管理器
        
        Args:
            hla_ref_file: HLA参考位点文件路径
        """
        self.hla_ref_file = hla_ref_file
        self.hla_alleles = set()
        self.supported_lengths = set()
        self.allele_length_pairs = set()
        self.allele_families = defaultdict(list)
        
        self._load_hla_reference()
    
    def _load_hla_reference(self):
        """加载HLA参考位点文件"""
        try:
            with open(self.hla_ref_file, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if line and ',' in line:
                        allele, length = line.split(',')
                        
                        # 标准化等位基因名称（移除HLA-前缀）
                        if allele.startswith('HLA-'):
                            allele = allele[4:]
                        
                        length = int(length)
                        
                        self.hla_alleles.add(allele)
                        self.supported_lengths.add(length)
                        self.allele_length_pairs.add((allele, length))
                        
                        # 按基因家族分组 (A, B, C)
                        family = allele.split('*')[0]
                        if allele not in self.allele_families[family]:
                            self.allele_families[family].append(allele)
        
        except FileNotFoundError:
            print(f"Warning: HLA reference file not found: {self.hla_ref_file}")
            self._load_default_alleles()
        except Exception as e:
            print(f"Error loading HLA reference file: {e}")
            self._load_default_alleles()
    
    def _load_default_alleles(self):
        """加载默认的HLA等位基因（仅A和B类）"""
        # 从hla_ref_set.class_i.txt文件中提取的所有A和B类等位基因
        default_alleles = [
            # HLA-A类等位基因
            'A*01:01', 'A*02:01', 'A*02:03', 'A*02:06', 'A*03:01',
            'A*11:01', 'A*23:01', 'A*24:02', 'A*26:01', 'A*30:01',
            'A*30:02', 'A*31:01', 'A*32:01', 'A*33:01', 'A*68:01', 'A*68:02',
            # HLA-B类等位基因
            'B*07:02', 'B*08:01', 'B*15:01', 'B*35:01', 'B*40:01',
            'B*44:02', 'B*44:03', 'B*51:01', 'B*53:01', 'B*57:01', 'B*58:01'
        ]
        
        for allele in default_alleles:
            self.hla_alleles.add(allele)
            family = allele.split('*')[0]
            if allele not in self.allele_families[family]:
                self.allele_families[family].append(allele)
        
        self.supported_lengths = {8, 9, 10, 11}
        
        # 生成所有组合
        for allele in default_alleles:
            for length in self.supported_lengths:
                self.allele_length_pairs.add((allele, length))
    
    def get_all_alleles(self):
        """获取所有支持的HLA等位基因"""
        return sorted(list(self.hla_alleles))
    
    def get_supported_lengths(self):
        """获取支持的肽段长度"""
        return sorted(list(self.supported_lengths))
    
    def get_alleles_by_family(self, family):
        """
        按基因家族获取等位基因
        
        Args:
            family: 基因家族 ('A', 'B', 'C')
        
        Returns:
            list: 该家族的等位基因列表
        """
        return sorted(self.allele_families.get(family, []))
    
    def is_supported_combination(self, allele, length):
        """
        检查等位基因和长度组合是否被支持
        
        Args:
            allele: HLA等位基因
            length: 肽段长度
        
        Returns:
            bool: 是否支持该组合
        """
        # 标准化等位基因名称
        if allele.startswith('HLA-'):
            allele = allele[4:]
        
        return (allele, length) in self.allele_length_pairs
    
    def validate_allele_string(self, allele_string):
        """
        验证和标准化HLA等位基因字符串（仅支持A和B类）
        
        Args:
            allele_string: 逗号分隔的HLA等位基因字符串
        
        Returns:
            tuple: (valid_alleles, invalid_alleles)
        """
        alleles = [allele.strip() for allele in allele_string.split(',')]
        valid_alleles = []
        invalid_alleles = []
        
        for allele in alleles:
            # 移除HLA-前缀（如果存在）
            if allele.startswith('HLA-'):
                allele = allele[4:]
            
            # 过滤C类等位基因
            if allele.startswith('C*'):
                invalid_alleles.append(allele)
                continue
            
            if allele in self.hla_alleles and (allele.startswith('A*') or allele.startswith('B*')):
                valid_alleles.append(allele)
            else:
                invalid_alleles.append(allele)
        
        return valid_alleles, invalid_alleles
    
    def get_representative_alleles(self, max_per_family=3):
        """
        获取每个基因家族的代表性等位基因
        
        Args:
            max_per_family: 每个家族最大等位基因数
        
        Returns:
            list: 代表性等位基因列表
        """
        representative = []
        
        for family in ['A', 'B']:  # 仅支持A和B类
            family_alleles = self.get_alleles_by_family(family)
            representative.extend(family_alleles[:max_per_family])
        
        return representative
    
    def filter_alleles_for_length(self, alleles, length):
        """
        过滤出支持指定长度的等位基因
        
        Args:
            alleles: 等位基因列表
            length: 目标肽段长度
        
        Returns:
            list: 支持该长度的等位基因
        """
        filtered = []
        for allele in alleles:
            if self.is_supported_combination(allele, length):
                filtered.append(allele)
        
        return filtered
    
    def suggest_alleles_for_population(self, population='European'):
        """
        根据人群推荐HLA等位基因
        
        Args:
            population: 目标人群
        
        Returns:
            list: 推荐的等位基因列表
        """
        populations = {
            'European': [
                'A*02:01', 'A*01:01', 'A*03:01', 'A*24:02',
                'B*07:02', 'B*08:01', 'B*44:02', 'B*35:01'
            ],
            'Asian': [
                'A*24:02', 'A*02:01', 'A*11:01', 'A*33:01',
                'B*58:01', 'B*15:01', 'B*44:03', 'B*51:01'
            ],
            'African': [
                'A*30:01', 'A*68:01', 'A*02:01', 'A*23:01',
                'B*53:01', 'B*15:01', 'B*58:01', 'B*07:02'
            ]
        }
        
        suggested = populations.get(population, populations['European'])
        
        # 过滤出实际支持的等位基因（仅A和B类）
        valid_alleles = []
        for allele in suggested:
            if allele in self.hla_alleles and (allele.startswith('A*') or allele.startswith('B*')):
                valid_alleles.append(allele)
        
        return valid_alleles
    
    def export_statistics(self, output_file):
        """
        导出HLA参考集统计信息
        
        Args:
            output_file: 输出文件路径
        """
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write("HLA Reference Set Statistics\n")
            f.write("=" * 40 + "\n\n")
            
            f.write(f"Total alleles: {len(self.hla_alleles)}\n")
            f.write(f"Supported lengths: {sorted(self.supported_lengths)}\n")
            f.write(f"Total allele-length combinations: {len(self.allele_length_pairs)}\n\n")
            
            f.write("Alleles by family:\n")
            for family in sorted(self.allele_families.keys()):
                alleles = self.allele_families[family]
                f.write(f"  {family}: {len(alleles)} alleles\n")
                for allele in sorted(alleles):
                    f.write(f"    {allele}\n")
                f.write("\n")
            
            f.write("Supported length distribution:\n")
            length_counts = defaultdict(int)
            for allele, length in self.allele_length_pairs:
                length_counts[length] += 1
            
            for length in sorted(length_counts.keys()):
                f.write(f"  {length}mer: {length_counts[length]} alleles\n")

def main():
    """测试函数"""
    # 创建HLA管理器
    hla_manager = HLAManager()
    
    print("HLA Reference Set Loaded:")
    print(f"Total alleles: {len(hla_manager.get_all_alleles())}")
    print(f"Supported lengths: {hla_manager.get_supported_lengths()}")
    print()
    
    # 按家族显示等位基因
    for family in ['A', 'B', 'C']:
        alleles = hla_manager.get_alleles_by_family(family)
        print(f"HLA-{family} family ({len(alleles)} alleles):")
        for allele in alleles:
            print(f"  {allele}")
        print()
    
    # 测试等位基因验证
    test_string = "HLA-A*02:01,B*07:02,C*07:02,A*99:99"
    valid, invalid = hla_manager.validate_allele_string(test_string)
    print(f"Validation test for: {test_string}")
    print(f"Valid: {valid}")
    print(f"Invalid: {invalid}")
    print()
    
    # 测试人群特异性推荐
    for population in ['European', 'Asian', 'African']:
        suggested = hla_manager.suggest_alleles_for_population(population)
        print(f"{population} population suggested alleles:")
        print(f"  {', '.join(suggested)}")
    print()
    
    # 测试长度过滤
    all_alleles = hla_manager.get_all_alleles()[:10]  # 前10个等位基因
    for length in [8, 9, 10, 11]:
        filtered = hla_manager.filter_alleles_for_length(all_alleles, length)
        print(f"Alleles supporting {length}mer: {len(filtered)}/{len(all_alleles)}")
    
    # 导出统计信息
    hla_manager.export_statistics("hla_reference_stats.txt")
    print("\nStatistics exported to: hla_reference_stats.txt")

if __name__ == "__main__":
    main()