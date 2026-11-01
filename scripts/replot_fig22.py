#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""以更小尺寸重新绘制 Fig22"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

chroms = ['Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5', 'Chr6', 'Chr7', 'Chr8']
gh28 = [11, 14, 4, 14, 4, 7, 14, 2]
ce8 = [16, 11, 5, 16, 2, 2, 6, 2]
pl1  = [6, 3, 1, 2, 1, 1, 2, 3]
gh17 = [5, 32, 4, 5, 1, 6, 3, 2]

# 缩小画布
fig, ax = plt.subplots(figsize=(8, 5))
x = np.arange(len(chroms))
width = 0.2

ax.bar(x - 1.5*width, gh28, width, label='GH28 (PG)', color='#E74C3C', alpha=0.85)
ax.bar(x - 0.5*width, ce8, width, label='CE8 (PME)', color='#F39C12', alpha=0.85)
ax.bar(x + 0.5*width, pl1, width, label='PL1 (Pectate lyase)', color='#9B59B6', alpha=0.85)
ax.bar(x + 1.5*width, gh17, width, label='GH17 (β-1,3-glucanase)', color='#1ABC9C', alpha=0.85)

ax.annotate('Chr2: 5 consecutive GH28\n(MY006181-185)\nat 13.87-13.98 Mb',
            xy=(1, gh28[1]), xytext=(1, 36),
            arrowprops=dict(arrowstyle='->', color='red', lw=1),
            ha='center', fontsize=8, color='red')
ax.annotate(f'Chr2: GH17 burst\n({gh17[1]} genes, 55%)',
            xy=(1, gh17[1]), xytext=(3.2, 36),
            arrowprops=dict(arrowstyle='->', color='green', lw=1),
            ha='center', fontsize=8, color='green')

# 标注数值
for i, vals in enumerate([gh28, ce8, pl1, gh17]):
    offset = (i - 1.5) * width
    for j, v in enumerate(vals):
        if v >= 5:
            ax.text(j + offset, v + 0.5, str(v), ha='center', va='bottom', fontsize=8)

ax.set_xlabel('Chromosome', fontsize=11)
ax.set_ylabel('Number of genes', fontsize=11)
ax.set_title('Distribution of key cell-wall hydrolase gene families across chromosomes',
             fontsize=12, weight='bold')
ax.set_xticks(x)
ax.set_xticklabels(chroms, fontsize=10)
ax.legend(fontsize=8, loc='upper right')
ax.grid(axis='y', alpha=0.3)
ax.set_ylim(0, 40)

# 强制小尺寸 + dpi=100
out = r"D:\xie's files\Xie's Data\李子数据\李子基因组、转录组测序信息\槜李-名媛\论文图表_v10\Fig22_cell_wall_chromosome_dist.png"
plt.tight_layout()
plt.savefig(out, dpi=100, bbox_inches='tight', pil_kwargs={'optimize': True})
plt.close()
print(f'Saved: {out}')
print(f'  Size: {os.path.getsize(out)/1024:.1f} KB')

from PIL import Image
img = Image.open(out)
print(f'  Dimensions: {img.size}')
# 强制限制不超过 1000x800
if img.size[0] > 1000 or img.size[1] > 1000:
    img.thumbnail((1000, 1000))
    img.save(out, optimize=True)
    print(f'  Resized to: {img.size}')