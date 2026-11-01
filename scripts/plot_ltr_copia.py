#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""绘制 LTR/Copia 在染色体上的详细分布图（与 LTR/Gypsy 对比）"""
import os
from collections import Counter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

# LTR/Copia 染色体分布（来自之前统计）
copia_chroms = {'Chr1': 5630, 'Chr2': 4273, 'Chr3': 3929, 'Chr4': 3608,
                'Chr5': 3578, 'Chr6': 3569, 'Chr7': 3556, 'Chr8': 2299}
gypsy_chroms = {'Chr1': 9466, 'Chr2': 8882, 'Chr3': 9272, 'Chr4': 7391,
                'Chr5': 7192, 'Chr6': 7458, 'Chr7': 6336, 'Chr8': 4312}

chroms = ['Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5', 'Chr6', 'Chr7', 'Chr8']
chr_len = {'Chr1': 53208954, 'Chr2': 38579706, 'Chr3': 36914270, 'Chr4': 35858089,
           'Chr5': 35300802, 'Chr6': 32489943, 'Chr7': 32115487, 'Chr8': 24114048}

# ============= 图 A: LTR/Copia vs LTR/Gypsy 数量与密度对比 =============
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

x = np.arange(len(chroms))
width = 0.35

# 左图: 数量
axes[0].bar(x - width/2, [gypsy_chroms[c] for c in chroms], width, label='LTR/Gypsy', color='#E74C3C', alpha=0.85)
axes[0].bar(x + width/2, [copia_chroms[c] for c in chroms], width, label='LTR/Copia', color='#3498DB', alpha=0.85)
axes[0].set_xlabel('Chromosome', fontsize=11)
axes[0].set_ylabel('Number of elements', fontsize=11)
axes[0].set_title('LTR/Gypsy vs LTR/Copia element counts', fontsize=12, weight='bold')
axes[0].set_xticks(x)
axes[0].set_xticklabels(chroms, fontsize=10)
axes[0].legend(fontsize=10)
axes[0].grid(axis='y', alpha=0.3)
for i, (g_v, c_v) in enumerate(zip([gypsy_chroms[c] for c in chroms], [copia_chroms[c] for c in chroms])):
    axes[0].text(i - width/2, g_v + 100, str(g_v), ha='center', va='bottom', fontsize=8)
    axes[0].text(i + width/2, c_v + 100, str(c_v), ha='center', va='bottom', fontsize=8)

# 右图: 密度（每Mb）
gypsy_dens = [gypsy_chroms[c] / (chr_len[c] / 1e6) for c in chroms]
copia_dens = [copia_chroms[c] / (chr_len[c] / 1e6) for c in chroms]
axes[1].bar(x - width/2, gypsy_dens, width, label='LTR/Gypsy', color='#E74C3C', alpha=0.85)
axes[1].bar(x + width/2, copia_dens, width, label='LTR/Copia', color='#3498DB', alpha=0.85)
axes[1].set_xlabel('Chromosome', fontsize=11)
axes[1].set_ylabel('Density (elements/Mb)', fontsize=11)
axes[1].set_title('LTR/Gypsy vs LTR/Copia density', fontsize=12, weight='bold')
axes[1].set_xticks(x)
axes[1].set_xticklabels(chroms, fontsize=10)
axes[1].legend(fontsize=10)
axes[1].grid(axis='y', alpha=0.3)
for i, (g_d, c_d) in enumerate(zip(gypsy_dens, copia_dens)):
    axes[1].text(i - width/2, g_d + 2, f'{g_d:.0f}', ha='center', va='bottom', fontsize=8)
    axes[1].text(i + width/2, c_d + 2, f'{c_d:.0f}', ha='center', va='bottom', fontsize=8)

plt.suptitle('Chromosomal distribution of LTR/Gypsy and LTR/Copia in Zuili plum',
             fontsize=13, weight='bold', y=1.02)
plt.tight_layout()
out1 = r"D:\xie's files\Xie's Data\李子数据\李子基因组、转录组测序信息\槜李-名媛\论文图表_v10\Fig27_LTR_Copia_dist.png"
plt.savefig(out1, dpi=120, bbox_inches='tight')
plt.close()
print(f'Saved Fig27: {out1}')

# ============= 图 B: LTR/Gypsy vs Copia 的 Gypsy/Copia 比值 =============
fig, ax = plt.subplots(figsize=(10, 5.5))
ratios = [gypsy_chroms[c] / max(copia_chroms[c], 1) for c in chroms]
colors = plt.cm.viridis(np.linspace(0, 0.8, len(chroms)))
bars = ax.bar(chroms, ratios, color=colors, alpha=0.85)

for i, (b, r) in enumerate(zip(bars, ratios)):
    ax.text(b.get_x() + b.get_width()/2, r + 0.05, f'{r:.2f}:1',
            ha='center', va='bottom', fontsize=10, weight='bold')

# 加平均线
avg = sum(gypsy_chroms[c] for c in chroms) / sum(copia_chroms[c] for c in chroms)
ax.axhline(y=avg, color='red', linestyle='--', alpha=0.6,
           label=f'Genome average: {avg:.2f}:1')

ax.set_xlabel('Chromosome', fontsize=12)
ax.set_ylabel('LTR/Gypsy : LTR/Copia ratio', fontsize=12)
ax.set_title('Gypsy/Copia ratio across chromosomes in Zuili plum genome',
             fontsize=13, weight='bold')
ax.legend(fontsize=11, loc='upper right')
ax.grid(axis='y', alpha=0.3)
ax.set_ylim(0, max(ratios) * 1.2)

plt.tight_layout()
out2 = r"D:\xie's files\Xie's Data\李子数据\李子基因组、转录组测序信息\槜李-名媛\论文图表_v10\Fig28_LTR_GC_ratio.png"
plt.savefig(out2, dpi=120, bbox_inches='tight')
plt.close()
print(f'Saved Fig28: {out2}')

# ============= 图 C: LTR/Copia 的 PercDiv 分布（与 Gypsy 对比） =============
fig, ax = plt.subplots(figsize=(10, 5.5))
bin_labels = ['0–5%', '5–10%', '10–15%', '15–20%', '20–25%', '25–30%', '30–50%', '>50%']
gypsy_divs = [19268, 8952, 7739, 8791, 7921, 5447, 2563, 0]
copia_divs = [3218, 4360, 5802, 6346, 6009, 3727, 1116, 0]

x = np.arange(len(bin_labels))
width = 0.35
ax.bar(x - width/2, gypsy_divs, width, label='LTR/Gypsy', color='#E74C3C', alpha=0.85)
ax.bar(x + width/2, copia_divs, width, label='LTR/Copia', color='#3498DB', alpha=0.85)

# 标注关键数值
for i, (g, c) in enumerate(zip(gypsy_divs, copia_divs)):
    if g > 5000:
        ax.text(i - width/2, g + 200, str(g), ha='center', va='bottom', fontsize=9, weight='bold')
    if c > 4000:
        ax.text(i + width/2, c + 200, str(c), ha='center', va='bottom', fontsize=9, weight='bold')

ax.set_xlabel('PercDiv (% divergence from consensus)', fontsize=12)
ax.set_ylabel('Number of elements', fontsize=12)
ax.set_title('LTR/Gypsy and LTR/Copia divergence distribution in Zuili plum',
             fontsize=13, weight='bold')
ax.set_xticks(x)
ax.set_xticklabels(bin_labels, rotation=30, ha='right')
ax.legend(fontsize=11)
ax.grid(axis='y', alpha=0.3)

# 注释关键区间
ax.annotate('Copia earlier peak\n(15-20% PercDiv,\n~6-7 Mya)',
            xy=(3, copia_divs[3]), xytext=(4.5, max(copia_divs)*0.7),
            arrowprops=dict(arrowstyle='->', color='blue'),
            fontsize=10, color='blue')
ax.annotate('Gypsy younger peak\n(0-5% PercDiv,\nrecent burst)',
            xy=(0, gypsy_divs[0]), xytext=(1, max(gypsy_divs)*0.85),
            arrowprops=dict(arrowstyle='->', color='red'),
            fontsize=10, color='red')

plt.tight_layout()
out3 = r"D:\xie's files\Xie's Data\李子数据\李子基因组、转录组测序信息\槜李-名媛\论文图表_v10\Fig29_LTR_divergence_compare.png"
plt.savefig(out3, dpi=120, bbox_inches='tight')
plt.close()
print(f'Saved Fig29: {out3}')

print('\n所有LTR/Copia相关图生成完成')