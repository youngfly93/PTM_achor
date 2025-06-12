#!/usr/bin/env python3
"""
Meta分析模块 - 合并多个批次的效应量
实现固定效应和随机效应模型
"""

import pandas as pd
import numpy as np
import scipy.stats as ss
import matplotlib.pyplot as plt
import seaborn as sns


def fixed_effect(df):
    """
    固定效应模型meta分析
    
    Args:
        df: DataFrame包含logOR和SE列
    
    Returns:
        tuple: (合并logOR, 标准误, Z值, p值)
    """
    # 计算权重（方差的倒数）
    weights = 1 / (df['SE'] ** 2)
    
    # 加权平均
    pooled_logOR = np.sum(weights * df['logOR']) / np.sum(weights)
    
    # 合并标准误
    pooled_SE = np.sqrt(1 / np.sum(weights))
    
    # Z值和p值
    z = pooled_logOR / pooled_SE
    p_value = 2 * (1 - ss.norm.cdf(abs(z)))
    
    return pooled_logOR, pooled_SE, z, p_value


def random_effect(df):
    """
    DerSimonian-Laird随机效应模型，返回全局logOR, 95%CI, I²
    
    Args:
        df: DataFrame包含logOR和SE列
    
    Returns:
        tuple: (合并logOR, 标准误, I²百分比, τ², Q统计量, p值)
    """
    if len(df) == 0:
        return np.nan, np.nan, 0, 0, 0, 1.0
    
    if len(df) == 1:
        # 只有一个研究时，直接返回该研究的结果
        return df.iloc[0]['logOR'], df.iloc[0]['SE'], 0, 0, 0, 1.0
    
    # 先计算固定效应估计
    w = 1 / (df['SE'] ** 2)
    fixed = np.sum(w * df['logOR']) / w.sum()
    
    # 计算异质性Q统计量
    Q = np.sum(w * (df['logOR'] - fixed) ** 2)
    df_deg = len(df) - 1
    
    # 计算τ²（研究间方差）
    tau2 = max(0, (Q - df_deg) / (w.sum() - (w**2).sum() / w.sum()))
    
    # 计算随机效应权重
    w_star = 1 / (df['SE'] ** 2 + tau2)
    
    # 随机效应合并估计
    random = np.sum(w_star * df['logOR']) / w_star.sum()
    se_re = np.sqrt(1 / w_star.sum())
    
    # 计算I²（异质性百分比）
    I2 = max(0, (Q - df_deg) / Q) * 100 if Q > 0 else 0
    
    # 计算p值
    z = random / se_re
    p_value = 2 * (1 - ss.norm.cdf(abs(z)))
    
    return random, se_re, I2, tau2, Q, p_value


def meta_analysis_by_mod(batch_data, method='random'):
    """
    按修饰类型进行meta分析
    
    Args:
        batch_data: 包含所有批次数据的DataFrame
        method: 'fixed' 或 'random'
    
    Returns:
        pd.DataFrame: meta分析结果
    """
    meta_results = []
    
    # 按修饰类型分组
    for mod, group in batch_data.groupby('mod'):
        # 过滤掉缺失值
        valid_group = group.dropna(subset=['logOR', 'SE'])
        
        if len(valid_group) == 0:
            continue
        
        if method == 'random':
            logOR, SE, I2, tau2, Q, p_value = random_effect(valid_group)
            
            meta_results.append({
                'mod': mod,
                'n_batches': len(valid_group),
                'logOR': logOR,
                'SE': SE,
                'CI_lower': logOR - 1.96 * SE,
                'CI_upper': logOR + 1.96 * SE,
                'OR': np.exp(logOR),
                'OR_CI_lower': np.exp(logOR - 1.96 * SE),
                'OR_CI_upper': np.exp(logOR + 1.96 * SE),
                'p_value': p_value,
                'I2': I2,
                'tau2': tau2,
                'Q': Q,
                'method': 'random'
            })
        else:  # fixed effect
            logOR, SE, z, p_value = fixed_effect(valid_group)
            
            meta_results.append({
                'mod': mod,
                'n_batches': len(valid_group),
                'logOR': logOR,
                'SE': SE,
                'CI_lower': logOR - 1.96 * SE,
                'CI_upper': logOR + 1.96 * SE,
                'OR': np.exp(logOR),
                'OR_CI_lower': np.exp(logOR - 1.96 * SE),
                'OR_CI_upper': np.exp(logOR + 1.96 * SE),
                'p_value': p_value,
                'method': 'fixed'
            })
    
    results_df = pd.DataFrame(meta_results)
    
    # 多重检验校正
    if len(results_df) > 0:
        # Bonferroni校正
        results_df['p_bonferroni'] = results_df['p_value'] * len(results_df)
        results_df['p_bonferroni'] = results_df['p_bonferroni'].clip(upper=1.0)
        
        # FDR校正（Benjamini-Hochberg）
        from statsmodels.stats.multitest import multipletests
        _, results_df['p_fdr'], _, _ = multipletests(
            results_df['p_value'], method='fdr_bh'
        )
    
    # 按p值排序
    results_df = results_df.sort_values('p_value')
    
    return results_df


def forest_plot(meta_results, batch_data, mod_name, output_file=None):
    """
    生成森林图
    
    Args:
        meta_results: meta分析结果（包含该修饰的合并效应）
        batch_data: 原始批次数据（用于显示各批次效应）
        mod_name: 修饰类型名称
        output_file: 输出文件路径（可选）
    """
    # 筛选特定修饰的数据
    mod_batch_data = batch_data[batch_data['mod'] == mod_name].copy()
    mod_meta = meta_results[meta_results['mod'] == mod_name].iloc[0]
    
    # 准备绘图数据
    n_studies = len(mod_batch_data)
    
    # 创建图形
    fig, ax = plt.subplots(1, 1, figsize=(10, 6 + n_studies * 0.3))
    
    # Y轴位置
    y_positions = list(range(n_studies + 2))  # +2 for meta result and spacing
    
    # 绘制各批次的效应量和置信区间
    for i, (_, row) in enumerate(mod_batch_data.iterrows()):
        y_pos = y_positions[i]
        
        # 计算置信区间
        ci_lower = row['logOR'] - 1.96 * row['SE']
        ci_upper = row['logOR'] + 1.96 * row['SE']
        
        # 绘制置信区间线
        ax.plot([ci_lower, ci_upper], [y_pos, y_pos], 'k-', linewidth=1)
        
        # 绘制点估计（大小反映权重）
        weight = 1 / (row['SE'] ** 2)
        marker_size = np.sqrt(weight) * 50 / np.sqrt(np.mean(1 / (mod_batch_data['SE'] ** 2)))
        ax.scatter(row['logOR'], y_pos, s=marker_size, c='blue', marker='s')
        
        # 添加批次标签
        ax.text(-3, y_pos, f"{row['batch_id']}", ha='right', va='center', fontsize=9)
    
    # 添加合并效应（菱形）
    y_meta = y_positions[-1]
    diamond_width = (mod_meta['CI_upper'] - mod_meta['CI_lower']) / 2
    diamond = plt.Polygon([
        (mod_meta['logOR'] - diamond_width, y_meta),
        (mod_meta['logOR'], y_meta + 0.2),
        (mod_meta['logOR'] + diamond_width, y_meta),
        (mod_meta['logOR'], y_meta - 0.2)
    ], fc='red', ec='red')
    ax.add_patch(diamond)
    
    # 添加标签
    ax.text(-3, y_meta, 'Overall', ha='right', va='center', fontweight='bold')
    
    # 添加零线
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    
    # 设置轴标签和标题
    ax.set_xlabel('Log Odds Ratio', fontsize=12)
    ax.set_title(f'Forest Plot: {mod_name} at Anchor Positions', fontsize=14)
    
    # 设置Y轴
    ax.set_ylim(-1, n_studies + 2)
    ax.set_yticks([])
    
    # 添加异质性统计信息
    heterogeneity_text = f"I² = {mod_meta['I2']:.1f}%, τ² = {mod_meta['tau2']:.3f}"
    if 'Q' in mod_meta:
        heterogeneity_text += f", Q = {mod_meta['Q']:.2f}"
    ax.text(0.02, 0.02, heterogeneity_text, transform=ax.transAxes, 
            fontsize=10, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # 调整布局
    plt.tight_layout()
    
    # 保存或显示
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Forest plot saved to: {output_file}")
    else:
        plt.show()
    
    plt.close()


def export_meta_results(meta_results, output_file):
    """
    导出meta分析结果
    
    Args:
        meta_results: meta分析结果DataFrame
        output_file: 输出文件路径
    """
    # 选择要导出的列
    export_columns = [
        'mod', 'n_batches', 'OR', 'OR_CI_lower', 'OR_CI_upper',
        'logOR', 'SE', 'p_value', 'p_bonferroni', 'p_fdr',
        'I2', 'tau2', 'method'
    ]
    
    # 确保所有列都存在
    available_columns = [col for col in export_columns if col in meta_results.columns]
    
    # 导出
    meta_results[available_columns].to_csv(output_file, index=False)
    
    print(f"Meta-analysis results exported to: {output_file}")
    print(f"Total modifications analyzed: {len(meta_results)}")
    
    # 打印摘要
    significant = meta_results[meta_results['p_fdr'] < 0.05]
    if len(significant) > 0:
        print(f"\nSignificant modifications (FDR < 0.05):")
        for _, row in significant.iterrows():
            print(f"  - {row['mod']}: OR = {row['OR']:.2f} "
                  f"({row['OR_CI_lower']:.2f}-{row['OR_CI_upper']:.2f}), "
                  f"p_fdr = {row['p_fdr']:.4f}, I² = {row.get('I2', 0):.1f}%")
    else:
        print("\nNo modifications reached significance after FDR correction.")


def create_summary_forest_plot(meta_results, output_file=None, top_n=20):
    """
    创建汇总森林图，显示所有修饰的meta分析结果
    
    Args:
        meta_results: meta分析结果DataFrame
        output_file: 输出文件路径
        top_n: 显示前N个最显著的修饰
    """
    # 选择要显示的修饰（按p值排序）
    plot_data = meta_results.nsmallest(top_n, 'p_value').sort_values('logOR')
    
    # 创建图形
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    
    # Y轴位置
    y_pos = np.arange(len(plot_data))
    
    # 绘制置信区间和点估计
    for i, (_, row) in enumerate(plot_data.iterrows()):
        # 置信区间
        ax.plot([row['CI_lower'], row['CI_upper']], [y_pos[i], y_pos[i]], 
                'k-', linewidth=1)
        
        # 点估计（颜色表示显著性）
        color = 'red' if row['p_fdr'] < 0.05 else 'blue'
        ax.scatter(row['logOR'], y_pos[i], s=100, c=color, marker='o')
        
        # 添加I²信息
        i2_text = f"I²={row.get('I2', 0):.0f}%"
        ax.text(row['CI_upper'] + 0.1, y_pos[i], i2_text, 
                va='center', fontsize=8, color='gray')
    
    # 添加零线
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    
    # 设置Y轴标签
    ax.set_yticks(y_pos)
    ax.set_yticklabels(plot_data['mod'])
    
    # 设置轴标签和标题
    ax.set_xlabel('Log Odds Ratio (95% CI)', fontsize=12)
    ax.set_title('Meta-Analysis: Modification Enrichment at Anchor Positions', fontsize=14)
    
    # 添加图例
    red_patch = plt.Line2D([0], [0], marker='o', color='w', 
                          markerfacecolor='r', markersize=10, label='FDR < 0.05')
    blue_patch = plt.Line2D([0], [0], marker='o', color='w', 
                           markerfacecolor='b', markersize=10, label='FDR ≥ 0.05')
    ax.legend(handles=[red_patch, blue_patch], loc='lower right')
    
    # 调整布局
    plt.tight_layout()
    
    # 保存或显示
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Summary forest plot saved to: {output_file}")
    else:
        plt.show()
    
    plt.close()