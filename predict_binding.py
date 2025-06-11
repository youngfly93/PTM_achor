#!/usr/bin/env python3
"""
HLA结合预测模块 - 使用MHCflurry预测肽段与HLA-I的结合能力，并标注锚位
集成HLA参考集管理，支持更多等位基因和长度组合
"""

import warnings
warnings.filterwarnings('ignore')

from hla_manager import HLAManager

def get_mhcflurry_predictor():
    """获取MHCflurry预测器"""
    try:
        from mhcflurry import Class1PresentationPredictor
        predictor = Class1PresentationPredictor.load()
        return predictor
    except ImportError:
        print("MHCflurry not installed. Please install with: pip install mhcflurry")
        return None
    except Exception as e:
        print(f"Error loading MHCflurry: {e}")
        return None

def predict_with_mhcflurry(sequences, alleles):
    """
    使用MHCflurry预测肽段结合
    
    Args:
        sequences: 肽段序列列表
        alleles: HLA等位基因列表
    
    Returns:
        pandas.DataFrame: 预测结果
    """
    predictor = get_mhcflurry_predictor()
    if predictor is None:
        return None
    
    try:
        results = predictor.predict(
            peptides=sequences,
            alleles=alleles,
            include_affinity_percentile=True
        )
        return results
    except Exception as e:
        print(f"Error in MHCflurry prediction: {e}")
        return None

def top_allele(seq, alleles, predictor=None):
    """
    为单个肽段找到最佳结合的HLA等位基因
    
    Args:
        seq: 肽段序列
        alleles: HLA等位基因列表
        predictor: MHCflurry预测器实例
    
    Returns:
        tuple: (最佳等位基因, presentation_score)
    """
    if predictor is None:
        predictor = get_mhcflurry_predictor()
        if predictor is None:
            return alleles[0], 0.0  # fallback
    
    try:
        results = predictor.predict(
            peptides=[seq], 
            alleles=alleles,
            include_affinity_percentile=True
        )
        
        if len(results) == 0:
            return alleles[0], 0.0
        
        # 按presentation_score排序（高分更好）
        best = results.sort_values("presentation_score", ascending=False).iloc[0]
        return best["allele"], best["presentation_score"]
    
    except Exception as e:
        print(f"Error predicting binding for {seq}: {e}")
        return alleles[0], 0.0

def annotate_anchor_positions(seq, allele=None):
    """
    标注肽段的锚位位置
    
    对于HLA-I类分子，主要锚位通常是：
    - P2 (位置2): N-terminal附近的锚位
    - PΩ (C-terminal): 最后一个位置
    
    Args:
        seq: 肽段序列
        allele: HLA等位基因（可选，用于更精确的锚位预测）
    
    Returns:
        set: 锚位位置集合（1-based indexing）
    """
    # 简化版本：对所有HLA-I等位基因使用相同的锚位规则
    # 实际应用中可以根据具体等位基因调整
    
    seq_len = len(seq)
    anchors = set()
    
    # P2 锚位（位置2）
    if seq_len >= 2:
        anchors.add(2)
    
    # C-terminal 锚位（最后一个位置）
    anchors.add(seq_len)
    
    # 对于某些等位基因，可能还有其他锚位
    # 这里可以根据等位基因特异性进行扩展
    if allele and seq_len >= 9:
        # 一些等位基因在位置9也有次要锚位
        if seq_len >= 9:
            anchors.add(seq_len - 1)  # PΩ-1位置
    
    return anchors

def get_hla_motif_anchors(allele, peptide_length):
    """
    基于HLA等位基因获取特异性锚位
    这是一个扩展函数，可以根据已知的HLA motif数据来确定锚位
    
    Args:
        allele: HLA等位基因名称
        peptide_length: 肽段长度
    
    Returns:
        set: 锚位位置集合
    """
    # 这里可以实现更复杂的等位基因特异性锚位规则
    # 目前使用简化版本
    anchors = {2, peptide_length}  # P2 和 PΩ
    
    # 一些特殊规则示例（可以根据文献扩展）
    if allele.startswith('A*02'):
        # A*02:01 family特异性规则
        if peptide_length >= 9:
            anchors.add(peptide_length - 1)  # 添加PΩ-1
    elif allele.startswith('B*07'):
        # B*07:02 family特异性规则
        pass  # 保持默认
    
    return anchors

def validate_hla_allele(allele):
    """
    验证HLA等位基因名称格式
    
    Args:
        allele: HLA等位基因名称
    
    Returns:
        bool: 是否为有效格式
    """
    import re
    # 标准HLA命名格式：A*02:01, B*07:02, C*07:02等
    pattern = r'^[ABC]\*\d{2}:\d{2}$'
    return bool(re.match(pattern, allele))

def parse_allele_string(allele_string, hla_manager=None):
    """
    解析HLA等位基因字符串（集成HLA管理器）
    
    Args:
        allele_string: 逗号分隔的HLA等位基因字符串，如"A*02:01,B*07:02,C*07:02"
        hla_manager: HLA管理器实例
    
    Returns:
        list: HLA等位基因列表
    """
    if hla_manager is None:
        hla_manager = HLAManager()
    
    valid_alleles, invalid_alleles = hla_manager.validate_allele_string(allele_string)
    
    if invalid_alleles:
        print(f"Warning: Invalid HLA alleles detected: {invalid_alleles}")
        print(f"Using valid alleles: {valid_alleles}")
    
    # 如果没有有效的等位基因，使用推荐的默认值
    if not valid_alleles:
        print("No valid alleles found, using European population defaults")
        valid_alleles = hla_manager.suggest_alleles_for_population('European')[:3]
    
    return valid_alleles

def batch_predict_binding(peptides, alleles_list, use_mhcflurry=True, hla_manager=None):
    """
    批量预测肽段与HLA的结合（集成HLA管理器）
    
    Args:
        peptides: 肽段信息列表 [{'Sequence': str, 'mod_list': list}, ...]
        alleles_list: HLA等位基因列表
        use_mhcflurry: 是否使用MHCflurry（如果False，使用简化评分）
        hla_manager: HLA管理器实例
    
    Returns:
        list: 预测结果列表
    """
    if hla_manager is None:
        hla_manager = HLAManager()
    
    results = []
    predictor = None
    
    if use_mhcflurry:
        predictor = get_mhcflurry_predictor()
        if predictor is None:
            print("Falling back to simplified scoring...")
            use_mhcflurry = False
    
    for i, peptide in enumerate(peptides):
        if (i + 1) % 100 == 0:
            print(f"Processed {i + 1}/{len(peptides)} peptides...")
        
        seq = peptide['Sequence']
        seq_length = len(seq)
        
        # 过滤支持当前肽段长度的等位基因
        compatible_alleles = hla_manager.filter_alleles_for_length(alleles_list, seq_length)
        
        if not compatible_alleles:
            # 如果没有兼容的等位基因，使用所有等位基因
            compatible_alleles = alleles_list
            if i == 0:  # 只在第一次警告
                print(f"Warning: No alleles support {seq_length}mer peptides, using all alleles")
        
        if use_mhcflurry:
            best_allele, score = top_allele(seq, compatible_alleles, predictor)
        else:
            # 简化评分：基于肽段长度和组成
            best_allele = compatible_alleles[0]
            score = len(seq) * 0.1  # 简单评分
        
        # 标注锚位
        anchors = annotate_anchor_positions(seq, best_allele)
        
        results.append({
            'sequence': seq,
            'best_allele': best_allele,
            'binding_score': score,
            'anchor_positions': anchors,
            'modifications': peptide['mod_list'],
            'compatible_alleles': len(compatible_alleles)
        })
    
    return results

def main():
    """测试函数"""
    # 创建HLA管理器
    hla_manager = HLAManager()
    
    print("HLA Manager Test:")
    print(f"Loaded {len(hla_manager.get_all_alleles())} alleles")
    print(f"Supported lengths: {hla_manager.get_supported_lengths()}")
    print()
    
    # 测试HLA等位基因解析
    test_alleles = "HLA-A*02:01,B*07:02,C*07:02,A*99:99"
    alleles = parse_allele_string(test_alleles, hla_manager)
    print(f"Parsed alleles from '{test_alleles}': {alleles}")
    print()
    
    # 测试人群特异性推荐
    european_alleles = hla_manager.suggest_alleles_for_population('European')[:3]
    print(f"European population alleles: {european_alleles}")
    
    # 测试锚位标注
    test_sequences = ["KFKESFAEM", "RILEMNDKYVK", "ADMAHISGL"]
    
    print("\nTesting anchor annotation:")
    for seq in test_sequences:
        anchors = annotate_anchor_positions(seq, alleles[0])
        print(f"{seq} (length {len(seq)}): anchors at positions {sorted(anchors)}")
    
    # 测试简化的批量预测
    peptides = [
        {'Sequence': 'KFKESFAEM', 'mod_list': [(9, 'Oxidation[M]')]},
        {'Sequence': 'RILEMNDKYVK', 'mod_list': []},
        {'Sequence': 'ADMAHISGL', 'mod_list': []},
        {'Sequence': 'PEPTIDE', 'mod_list': []},  # 7mer - 不支持
        {'Sequence': 'LONGERPEPTIDE', 'mod_list': []}  # 13mer - 不支持
    ]
    
    print("\nTesting batch prediction with HLA manager:")
    results = batch_predict_binding(peptides, alleles, use_mhcflurry=False, hla_manager=hla_manager)
    
    for result in results:
        print(f"Sequence: {result['sequence']} (length: {len(result['sequence'])})")
        print(f"  Best allele: {result['best_allele']}")
        print(f"  Score: {result['binding_score']:.3f}")
        print(f"  Anchors: {sorted(result['anchor_positions'])}")
        print(f"  Compatible alleles: {result['compatible_alleles']}")
        print(f"  Modifications: {result['modifications']}")
        print()

if __name__ == "__main__":
    main()