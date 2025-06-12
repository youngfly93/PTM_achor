#!/usr/bin/env python3
"""
测试HLA设置 - 验证修改后的HLA管理和预测功能
"""

from hla_manager import HLAManager
from predict_binding import parse_allele_string, validate_hla_allele

def test_hla_manager():
    """测试HLA管理器"""
    print("=" * 60)
    print("测试HLA管理器")
    print("=" * 60)
    
    hm = HLAManager()
    
    print(f"总等位基因数: {len(hm.get_all_alleles())}")
    print(f"A类等位基因数: {len(hm.get_alleles_by_family('A'))}")
    print(f"B类等位基因数: {len(hm.get_alleles_by_family('B'))}")
    print(f"C类等位基因数: {len(hm.get_alleles_by_family('C'))}")
    print(f"支持的肽段长度: {sorted(hm.get_supported_lengths())}")
    
    print("\nA类等位基因:")
    for allele in hm.get_alleles_by_family('A'):
        print(f"  {allele}")
    
    print("\nB类等位基因:")
    for allele in hm.get_alleles_by_family('B'):
        print(f"  {allele}")

def test_allele_validation():
    """测试等位基因验证"""
    print("\n" + "=" * 60)
    print("测试等位基因验证")
    print("=" * 60)
    
    hm = HLAManager()
    
    test_cases = [
        "A*02:01,B*07:02,C*07:02",  # 包含C类
        "A*02:01,B*07:02",          # 仅A和B类
        "A*02:01,A*01:01,B*07:02,B*08:01",  # 多个A和B类
        "C*07:02,C*07:01",          # 仅C类
        "A*99:99,B*99:99",          # 不存在的等位基因
        ""                          # 空字符串
    ]
    
    for i, allele_string in enumerate(test_cases, 1):
        print(f"\n测试 {i}: '{allele_string}'")
        valid, invalid = hm.validate_allele_string(allele_string)
        print(f"  有效: {valid}")
        print(f"  无效: {invalid}")

def test_population_suggestions():
    """测试人群特异性推荐"""
    print("\n" + "=" * 60)
    print("测试人群特异性推荐")
    print("=" * 60)
    
    hm = HLAManager()
    
    populations = ['European', 'Asian', 'African']
    
    for pop in populations:
        suggested = hm.suggest_alleles_for_population(pop)
        print(f"\n{pop}人群推荐等位基因:")
        for allele in suggested:
            print(f"  {allele}")

def test_parsing_function():
    """测试解析函数"""
    print("\n" + "=" * 60)
    print("测试解析函数")
    print("=" * 60)
    
    test_strings = [
        "A*02:01,B*07:02,C*07:02",
        "A*02:01,B*07:02",
        "C*07:02"
    ]
    
    for allele_string in test_strings:
        print(f"\n解析: '{allele_string}'")
        result = parse_allele_string(allele_string)
        print(f"  结果: {result}")

def test_individual_validation():
    """测试单个等位基因验证"""
    print("\n" + "=" * 60)
    print("测试单个等位基因验证")
    print("=" * 60)
    
    test_alleles = [
        "A*02:01",   # 有效A类
        "B*07:02",   # 有效B类
        "C*07:02",   # C类（应该无效）
        "A*99:99",   # 无效格式
        "invalid"    # 完全无效
    ]
    
    for allele in test_alleles:
        is_valid = validate_hla_allele(allele)
        print(f"  {allele}: {'有效' if is_valid else '无效'}")

def main():
    """运行所有测试"""
    print("HLA设置测试开始...\n")
    
    test_hla_manager()
    test_allele_validation()
    test_population_suggestions()
    test_parsing_function()
    test_individual_validation()
    
    print("\n" + "=" * 60)
    print("✅ 所有测试完成！")
    print("现在系统仅支持HLA-A和HLA-B类等位基因，")
    print("自动过滤HLA-C类等位基因，避免MHCflurry兼容性问题。")
    print("=" * 60)

if __name__ == "__main__":
    main()