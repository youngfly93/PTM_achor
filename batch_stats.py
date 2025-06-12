#!/usr/bin/env python3
"""
批次统计模块 - 按数据集计算修饰-锚位富集的效应量
用于后续meta分析
"""

import pandas as pd
import numpy as np
import scipy.stats as ss
from collections import defaultdict


def batch_enrich(df):
    """
    输入单批DataFrame，输出logOR、SE等meta-analysis所需量
    
    Args:
        df: 单个批次的DataFrame，包含mod和anchor_tag列
    
    Returns:
        pd.DataFrame: 包含每种修饰的效应量统计
            - mod: 修饰类型
            - a: 修饰在锚位的数量
            - n: 修饰在非锚位的数量
            - rest_a: 其他修饰在锚位的数量
            - rest_n: 其他修饰在非锚位的数量
            - logOR: log odds ratio
            - SE: 标准误
            - anchor_pct: 修饰在锚位的百分比
    """
    # 构建2x2列联表
    tbl = df.groupby(['mod', 'anchor_tag']).size().unstack(fill_value=0)
    
    out = []
    
    for mod in tbl.index:
        # 当前修饰的锚位/非锚位计数
        a = tbl.loc[mod, 'anchor'] if 'anchor' in tbl.columns else 0
        n = tbl.loc[mod, 'non_anchor'] if 'non_anchor' in tbl.columns else 0
        
        # 其他所有修饰的锚位/非锚位计数
        rest_a = tbl['anchor'].sum() - a if 'anchor' in tbl.columns else 0
        rest_n = tbl['non_anchor'].sum() - n if 'non_anchor' in tbl.columns else 0
        
        # 计算odds ratio和标准误
        # 添加0.5校正避免零值
        a_corr = a + 0.5
        n_corr = n + 0.5
        rest_a_corr = rest_a + 0.5
        rest_n_corr = rest_n + 0.5
        
        or_ = (a_corr / n_corr) / (rest_a_corr / rest_n_corr)
        se = np.sqrt(1/a_corr + 1/n_corr + 1/rest_a_corr + 1/rest_n_corr)
        
        # 计算Fisher精确检验p值
        p_value = ss.fisher_exact([[a, n], [rest_a, rest_n]])[1]
        
        # 计算锚位百分比
        anchor_pct = (a / (a + n)) * 100 if (a + n) > 0 else 0
        
        out.append({
            'mod': mod,
            'a': a,
            'n': n,
            'rest_a': rest_a,
            'rest_n': rest_n,
            'logOR': np.log(or_),
            'SE': se,
            'anchor_pct': anchor_pct,
            'p_value': p_value,
            'total': a + n
        })
    
    return pd.DataFrame(out)


def batch_group_comparison(tumor_df, normal_df):
    """
    比较单个批次内的肿瘤vs正常差异
    
    Args:
        tumor_df: 肿瘤样本的DataFrame
        normal_df: 正常样本的DataFrame
    
    Returns:
        pd.DataFrame: 包含组间比较结果
    """
    # 分别计算肿瘤和正常的富集
    tumor_stats = batch_enrich(tumor_df)
    normal_stats = batch_enrich(normal_df)
    
    # 合并结果
    merged = pd.merge(
        tumor_stats, 
        normal_stats, 
        on='mod', 
        how='outer',
        suffixes=('_tumor', '_normal')
    ).fillna(0)
    
    # 计算组间差异
    results = []
    for _, row in merged.iterrows():
        # 组间Fisher精确检验
        table = [
            [row['a_tumor'], row['n_tumor']],
            [row['a_normal'], row['n_normal']]
        ]
        
        if sum([sum(r) for r in table]) > 0:
            p_group = ss.fisher_exact(table)[1]
        else:
            p_group = 1.0
        
        # 计算组间效应量（log odds ratio）
        a_t = row['a_tumor'] + 0.5
        n_t = row['n_tumor'] + 0.5
        a_n = row['a_normal'] + 0.5
        n_n = row['n_normal'] + 0.5
        
        group_or = (a_t / n_t) / (a_n / n_n)
        group_logOR = np.log(group_or)
        group_se = np.sqrt(1/a_t + 1/n_t + 1/a_n + 1/n_n)
        
        results.append({
            'mod': row['mod'],
            'tumor_anchor_pct': row['anchor_pct_tumor'],
            'normal_anchor_pct': row['anchor_pct_normal'],
            'group_logOR': group_logOR,
            'group_SE': group_se,
            'p_group': p_group,
            'total_tumor': row['total_tumor'],
            'total_normal': row['total_normal']
        })
    
    return pd.DataFrame(results)


def export_batch_stats(batch_results, output_file):
    """
    导出批次统计结果
    
    Args:
        batch_results: 批次统计结果列表
        output_file: 输出文件路径
    """
    # 合并所有批次结果
    all_results = pd.concat(batch_results, ignore_index=True)
    
    # 按批次和修饰类型排序
    all_results = all_results.sort_values(['batch_id', 'mod'])
    
    # 导出到CSV
    all_results.to_csv(output_file, index=False)
    
    print(f"Batch statistics exported to: {output_file}")
    print(f"Total batches: {all_results['batch_id'].nunique()}")
    print(f"Total modifications: {all_results['mod'].nunique()}")
    
    return all_results


def summarize_batch_heterogeneity(batch_results):
    """
    评估批次间的异质性
    
    Args:
        batch_results: 包含所有批次结果的DataFrame
    
    Returns:
        dict: 异质性统计摘要
    """
    summary = {}
    
    # 按修饰类型分组
    for mod, group in batch_results.groupby('mod'):
        if len(group) < 2:
            continue
            
        # 计算批次间变异
        logOR_values = group['logOR'].values
        se_values = group['SE'].values
        
        # 计算Q统计量
        weights = 1 / (se_values ** 2)
        weighted_mean = np.sum(weights * logOR_values) / np.sum(weights)
        Q = np.sum(weights * (logOR_values - weighted_mean) ** 2)
        
        # 计算I²
        df = len(group) - 1
        I2 = max(0, (Q - df) / Q) * 100 if Q > 0 else 0
        
        # 计算τ²（批次间方差）
        tau2 = max(0, (Q - df) / (np.sum(weights) - np.sum(weights**2) / np.sum(weights)))
        
        summary[mod] = {
            'n_batches': len(group),
            'Q': Q,
            'I2': I2,
            'tau2': tau2,
            'mean_logOR': weighted_mean,
            'range_logOR': (logOR_values.min(), logOR_values.max())
        }
    
    return summary