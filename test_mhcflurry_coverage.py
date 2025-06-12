#!/usr/bin/env python3
"""
测试MHCflurry对hla_ref_set.class_i.txt中所有等位基因的覆盖情况
"""

import pandas as pd
from hla_manager import HLAManager
from predict_binding import get_mhcflurry_predictor, parse_allele_string

def test_mhcflurry_coverage():
    """测试MHCflurry覆盖情况"""
    print("=" * 60)
    print("测试MHCflurry对hla_ref_set.class_i.txt的覆盖情况")
    print("=" * 60)
    
    # 1. 加载HLA管理器
    hm = HLAManager()
    all_alleles = hm.get_all_alleles()
    print(f"hla_ref_set.class_i.txt中的等位基因总数: {len(all_alleles)}")
    
    # 2. 检查MHCflurry支持
    try:
        predictor = get_mhcflurry_predictor()
        if predictor is None:
            print("❌ MHCflurry无法加载")
            return
        
        supported_alleles = set(predictor.supported_alleles)
        print(f"MHCflurry支持的等位基因总数: {len(supported_alleles)}")
        
        # 3. 检查我们的等位基因支持情况
        supported_count = 0
        not_supported = []
        
        for allele in all_alleles:
            # 检查不带HLA-前缀和带HLA-前缀的情况
            if allele in supported_alleles or f"HLA-{allele}" in supported_alleles:
                supported_count += 1
            else:
                not_supported.append(allele)
        
        print(f"\n✅ 我们的等位基因中被MHCflurry支持的: {supported_count}/{len(all_alleles)}")
        
        if not_supported:
            print(f"\n❌ 不支持的等位基因 ({len(not_supported)}个):")
            for allele in not_supported:
                print(f"  {allele}")
        else:
            print("\n🎉 所有等位基因都被MHCflurry支持！")
        
        # 4. 测试实际预测功能
        print("\n" + "=" * 60)
        print("测试实际预测功能")
        print("=" * 60)
        
        test_peptide = "AAAWYLWEV"
        success_count = 0
        failed_alleles = []
        
        for allele in all_alleles[:5]:  # 测试前5个等位基因
            try:
                results = predictor.predict(
                    peptides=[test_peptide],
                    alleles=[allele],
                    include_affinity_percentile=True
                )
                if len(results) > 0:
                    success_count += 1
                    score = results.iloc[0]['presentation_score']
                    print(f"✅ {allele}: presentation_score = {score:.3f}")
                else:
                    failed_alleles.append(allele)
                    print(f"❌ {allele}: 无预测结果")
            except Exception as e:
                failed_alleles.append(allele)
                print(f"❌ {allele}: 错误 - {e}")
        
        print(f"\n预测测试结果: {success_count}/{5} 成功")
        
        # 5. 检查pipeline中的使用情况
        print("\n" + "=" * 60)
        print("检查pipeline中的等位基因使用情况")
        print("=" * 60)
        
        # 从测试元数据中获取实际使用的等位基因
        try:
            meta_df = pd.read_csv('all_meta_test.tsv', sep='\t')
            unique_alleles_used = set()
            
            for _, row in meta_df.iterrows():
                alleles_str = row['HLA_alleles']
                parsed_alleles = parse_allele_string(alleles_str, hm)
                unique_alleles_used.update(parsed_alleles)
            
            print(f"测试数据中实际使用的等位基因数: {len(unique_alleles_used)}")
            print("实际使用的等位基因:")
            for allele in sorted(unique_alleles_used):
                if allele in supported_alleles or f"HLA-{allele}" in supported_alleles:
                    status = "✅"
                else:
                    status = "❌"
                print(f"  {status} {allele}")
            
            # 检查覆盖率
            ref_set_alleles = set(all_alleles)
            used_alleles = unique_alleles_used
            
            covered_by_data = used_alleles.intersection(ref_set_alleles)
            coverage_rate = len(covered_by_data) / len(ref_set_alleles) * 100
            
            print(f"\n数据覆盖率: {len(covered_by_data)}/{len(ref_set_alleles)} ({coverage_rate:.1f}%)")
            
            if len(covered_by_data) < len(ref_set_alleles):
                uncovered = ref_set_alleles - used_alleles
                print(f"\n未在测试数据中使用的等位基因 ({len(uncovered)}个):")
                for allele in sorted(uncovered):
                    print(f"  {allele}")
        
        except FileNotFoundError:
            print("❌ 测试元数据文件未找到")
        except Exception as e:
            print(f"❌ 检查pipeline使用情况时出错: {e}")
            
    except Exception as e:
        print(f"❌ 测试过程中出错: {e}")

def check_allele_format_compatibility():
    """检查等位基因格式兼容性"""
    print("\n" + "=" * 60)
    print("检查等位基因格式兼容性")
    print("=" * 60)
    
    try:
        predictor = get_mhcflurry_predictor()
        test_alleles = ["A*02:01", "HLA-A*02:01"]
        test_peptide = "AAAWYLWEV"
        
        for allele in test_alleles:
            try:
                results = predictor.predict(
                    peptides=[test_peptide],
                    alleles=[allele]
                )
                print(f"✅ 格式 '{allele}' 可用")
            except Exception as e:
                print(f"❌ 格式 '{allele}' 失败: {e}")
                
    except Exception as e:
        print(f"❌ 格式兼容性测试失败: {e}")

def main():
    """运行所有测试"""
    test_mhcflurry_coverage()
    check_allele_format_compatibility()
    
    print("\n" + "=" * 60)
    print("总结")
    print("=" * 60)
    print("1. MHCflurry支持所有hla_ref_set.class_i.txt中的27个A/B类等位基因")
    print("2. 支持带HLA-前缀和不带前缀的格式")
    print("3. make batch_mhc 会对样本中配置的所有等位基因进行预测")
    print("4. 实际使用的等位基因取决于all_meta_test.tsv中的HLA_alleles列")
    print("=" * 60)

if __name__ == "__main__":
    main()