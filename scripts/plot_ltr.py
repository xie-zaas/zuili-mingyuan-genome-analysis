#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""绘制 LTR/Gypsy vs LTR/Copia 的插入时间分布图"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# 数据 (来自 LTR 分析)
gypsy_data = [19268, 8952, 7739, 8791, 7921, 5447, 2563, 0]
copia_data = [3218, 4360, 5802, 6346, 6009, 3727, 1116, 0]
bin_labels = ['0–5%', '5–10%', '10–15%', '15–20%', '20–25%', '25–30%', '30–50%', '>50%']
bin_centers = [2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 40, 60]  # 用于插入时间换算
mu = 1.3e-9  # substitutions/site/year

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# 图1: PercDiv 分布
ax = axes[0]
x = np.arange(len(bin_labels))
width = 0.35
ax.bar(x - width/2, gypsy_data, width, label='LTR/Gypsy', color='#E74C3C')
ax.bar(x + width/2, copia_data, width, label='LTR/Copia', color='#3498DB')
ax.set_xlabel('PercDiv (% divergence)', fontsize=12)
ax.set_ylabel('Number of LTR-RT elements', fontsize=12)
ax.set_title('LTR/Gypsy vs LTR/Copia divergence distribution', fontsize=13, weight='bold')
ax.set_xticks(x)
ax.set_xticklabels(bin_labels, rotation=30)
ax.legend(fontsize=11)
ax.grid(axis='y', alpha=0.3)

# 图2: 推算的插入时间
ax = axes[1]
# 推算每个bin的中心插入年龄 (Mya)
times = [c / (2 * mu * 100) / 1e6 for c in bin_centers]  # Mya
ax.bar(x - width/2, gypsy_data, width, label='LTR/Gypsy', color='#E74C3C')
ax.bar(x + width/2, copia_data, width, label='LTR/Copia', color='#3498DB')
ax.set_xlabel('Estimated insertion age (Mya)', fontsize=12)
ax.set_ylabel('Number of LTR-RT elements', fontsize=12)
ax.set_title('Inferred insertion time of LTR-RT in Zuili plum', fontsize=13, weight='bold')
ax.set_xticks(x)
ax.set_xticklabels([f'{t:.1f}' for t in times], rotation=30)
ax.legend(fontsize=11)
ax.grid(axis='y', alpha=0.3)
ax.axvline(x=2.5-0.5, color='gray', linestyle='--', alpha=0.5)
ax.text(2.7, ax.get_ylim()[1]*0.7, 'Recent burst\n(<5 Mya)', color='gray', fontsize=10)

plt.tight_layout()
img_dir = r"D:\xie's files\Xie's Data\李子数据\李子基因组、转录组测序信息\槜李-名媛\论文图表_v10"
out1 = os.path.join(img_dir, 'Fig18_LTR_divergence.png')
plt.savefig(out1, dpi=200, bbox_inches='tight')
plt.close()
print(f'Saved: {out1}')

# 图2: LTR 类别柱状图（Gypsy/Copia/Cassandra 等）
fig, ax = plt.subplots(figsize=(10, 6))
classes = ['Gypsy', 'Copia', 'Cassandra', 'ERV1', 'ERVK', 'Caulimovirus', 'Pao']
counts = [60681, 30578, 1786, 2134, 1094, 971, 601]
colors = ['#E74C3C', '#3498DB', '#F39C12', '#9B59B6', '#1ABC9C', '#34495E', '#E67E22']
bars = ax.bar(classes, counts, color=colors)
ax.set_ylabel('Number of elements', fontsize=12)
ax.set_title('LTR-RT subfamily composition in Zuili plum genome', fontsize=13, weight='bold')
ax.set_yscale('log')
for b, v in zip(bars, counts):
    ax.text(b.get_x() + b.get_width()/2, v, f'{v:,}', ha='center', va='bottom', fontsize=10)
ax.grid(axis='y', alpha=0.3)
plt.xticks(rotation=15)
plt.tight_layout()
out2 = os.path.join(img_dir, 'Fig19_LTR_subfamily.png')
plt.savefig(out2, dpi=200, bbox_inches='tight')
plt.close()
print(f'Saved: {out2}')

# 图3: 关键基因家族汇总柱状图
fig, ax = plt.subplots(figsize=(10, 6))
families = ['GH28\n(Polygalacturonase)', 'CE8\n(PME)', 'PL1\n(Pectate lyase)', 'GH\n(Total)', 'PL\n(Total)',
            'LOX', 'AAT-like\n(BAHD)', 'TPS', 'ADH', 'SUT', 'INV', 'SUS', 'SPS']
counts = [70, 61, 17, 408, 26, 12, 8, 17, 24, 37, 8, 9, 11]
colors_fam = ['#E74C3C', '#F39C12', '#9B59B6', '#3498DB', '#1ABC9C',
              '#16A085', '#27AE60', '#2980B9', '#8E44AD', '#D35400', '#C0392B', '#7F8C8D', '#2C3E50']
bars = ax.bar(families, counts, color=colors_fam)
ax.set_ylabel('Number of genes', fontsize=12)
ax.set_title('Key gene families related to fruit quality & aroma in Zuili plum',
             fontsize=13, weight='bold')
for b, v in zip(bars, counts):
    ax.text(b.get_x() + b.get_width()/2, v, f'{v}', ha='center', va='bottom', fontsize=10)
ax.grid(axis='y', alpha=0.3)
plt.xticks(rotation=20, ha='right', fontsize=10)
plt.tight_layout()
out3 = os.path.join(img_dir, 'Fig20_key_gene_families.png')
plt.savefig(out3, dpi=200, bbox_inches='tight')
plt.close()
print(f'Saved: {out3}')

print('Done')