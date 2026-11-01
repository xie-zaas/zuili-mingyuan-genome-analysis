#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""绘制染色体上特殊区段分布图"""
import os
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter, defaultdict

# 中文字体处理
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.family'] = 'sans-serif'

base = r"D:\xie's files\Xie's Data\李子数据\李子基因组、转录组测序信息\槜李-名媛"
csv_path = os.path.join(base, r"论文图表_v10\染色体特殊区段位置信息表.csv")
img_dir = os.path.join(base, r"论文图表_v10")

# 读 CSV
data = []
with open(csv_path, 'r', encoding='utf-8-sig') as f:
    reader = csv.DictReader(f)
    for row in reader:
        data.append(row)

# 染色体顺序
chroms = ['Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5', 'Chr6', 'Chr7', 'Chr8']
chr_len = {'Chr1': 53208954, 'Chr2': 38579706, 'Chr3': 36914270, 'Chr4': 35858089,
           'Chr5': 35300802, 'Chr6': 32489943, 'Chr7': 32115487, 'Chr8': 24114048}

# 收集数据
def collect(feature_name):
    """返回 {chrom: count}"""
    res = {}
    for row in data:
        if row['特征名称'] == feature_name:
            try:
                n = int(row['数量'].replace(',', ''))
                res[row['染色体']] = n
            except:
                pass
    return res

ltr_gypsy = collect('LTR/Gypsy')
ltr_copia = collect('LTR/Copia')
ltr_cassandra = collect('LTR/Cassandra')
ltr_recent = collect('LTR/Gypsy 近期 (PercDiv<5%)')
gh28 = collect('GH28 (PG)')
ce8 = collect('CE8 (PME)')
pl1 = collect('PL1 (Pectate lyase)')
gh17 = collect('GH17 (β-1,3-glucanase)')
dna_pif = collect('DNA/PIF-Harbinger')
dna_mule = collect('DNA/MULE-MuDR')
dna_cmc = collect('DNA/CMC-EnSpm')
dna_hat = collect('DNA/hAT-Ac')
line_l1 = collect('LINE/L1')
rc_helitron = collect('RC/Helitron')

# ============= 图1: LTR 各类亚家族染色体分布 =============
fig, ax = plt.subplots(figsize=(12, 6))

x = np.arange(len(chroms))
width = 0.35

gypsy_vals = [ltr_gypsy.get(c, 0) for c in chroms]
copia_vals = [ltr_copia.get(c, 0) for c in chroms]

ax.bar(x - width/2, gypsy_vals, width, label='LTR/Gypsy', color='#E74C3C', alpha=0.85)
ax.bar(x + width/2, copia_vals, width, label='LTR/Copia', color='#3498DB', alpha=0.85)

# 标注 Chr3 LTR/Cassandra
cassandra_vals = [ltr_cassandra.get(c, 0) for c in chroms]
for i, (c, v) in enumerate(zip(chroms, cassandra_vals)):
    if v > 0:
        ax.text(i, max(gypsy_vals[i], copia_vals[i]) + 500, f'Cassandra={v}',
                ha='center', va='bottom', fontsize=8,
                color='darkorange', weight='bold')

ax.set_xlabel('Chromosome', fontsize=12)
ax.set_ylabel('Number of LTR-RT elements', fontsize=12)
ax.set_title('Distribution of LTR/Gypsy and LTR/Copia across chromosomes in Zuili plum', fontsize=13, weight='bold')
ax.set_xticks(x)
ax.set_xticklabels(chroms)
ax.legend(fontsize=11)
ax.grid(axis='y', alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(img_dir, 'Fig21_LTR_chromosome_dist.png'), dpi=120, bbox_inches='tight')
plt.close()
print('Saved Fig21')

# ============= 图2: 关键细胞壁酶染色体分布 =============
fig, ax = plt.subplots(figsize=(12, 6))
gh28_v = [gh28.get(c, 0) for c in chroms]
ce8_v = [ce8.get(c, 0) for c in chroms]
pl1_v = [pl1.get(c, 0) for c in chroms]
gh17_v = [gh17.get(c, 0) for c in chroms]

width = 0.2
ax.bar(x - 1.5*width, gh28_v, width, label='GH28 (PG)', color='#E74C3C', alpha=0.85)
ax.bar(x - 0.5*width, ce8_v, width, label='CE8 (PME)', color='#F39C12', alpha=0.85)
ax.bar(x + 0.5*width, pl1_v, width, label='PL1 (Pectate lyase)', color='#9B59B6', alpha=0.85)
ax.bar(x + 1.5*width, gh17_v, width, label='GH17 (β-1,3-glucanase)', color='#1ABC9C', alpha=0.85)

# 标注 GH28 在 Chr2 13.87–13.98 Mb 区间存在 5 个连续基因
ax.annotate('Chr2: 5 consecutive GH28\n(MY006181-185) at 13.87-13.98 Mb',
            xy=(1, gh28_v[1]), xytext=(1, 35),
            arrowprops=dict(arrowstyle='->', color='red'),
            ha='center', fontsize=9, color='red', weight='bold')
# 标注 GH17 在 Chr2 富集
ax.annotate(f'Chr2: GH17 burst\n({gh17_v[1]} genes, 55%)',
            xy=(1, gh17_v[1]), xytext=(2.5, 38),
            arrowprops=dict(arrowstyle='->', color='green'),
            ha='center', fontsize=9, color='green', weight='bold')

ax.set_xlabel('Chromosome', fontsize=12)
ax.set_ylabel('Number of genes', fontsize=12)
ax.set_title('Distribution of key cell-wall hydrolase gene families across chromosomes', fontsize=13, weight='bold')
ax.set_xticks(x)
ax.set_xticklabels(chroms)
ax.legend(fontsize=10, loc='upper right')
ax.grid(axis='y', alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(img_dir, 'Fig22_cell_wall_chromosome_dist.png'), dpi=120, bbox_inches='tight')
plt.close()
print('Saved Fig22')

# ============= 图3: KEGG 关键基因家族染色体分布 =============
kegg_features = ['LOX', 'TYRAAT (BAHD)', 'ADH', 'TPS', 'SUS', 'SPS', 'CHS', 'F3H']
kegg_data = {k: collect(k) for k in kegg_features}

fig, ax = plt.subplots(figsize=(14, 7))
colors_kegg = ['#E74C3C', '#9B59B6', '#F39C12', '#3498DB', '#1ABC9C', '#16A085', '#2C3E50', '#8E44AD']
n_kegg = len(kegg_features)
width = 0.10
for i, (k, color) in enumerate(zip(kegg_features, colors_kegg)):
    offset = (i - n_kegg/2 + 0.5) * width
    vals = [kegg_data[k].get(c, 0) for c in chroms]
    bars = ax.bar(x + offset, vals, width, label=k, color=color, alpha=0.85)

ax.set_xlabel('Chromosome', fontsize=12)
ax.set_ylabel('Number of genes', fontsize=12)
ax.set_title('Distribution of KEGG key gene families across chromosomes', fontsize=13, weight='bold')
ax.set_xticks(x)
ax.set_xticklabels(chroms)
ax.legend(fontsize=9, ncol=2, loc='upper right')
ax.grid(axis='y', alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(img_dir, 'Fig23_KEGG_chromosome_dist.png'), dpi=120, bbox_inches='tight')
plt.close()
print('Saved Fig23')

# ============= 图4: 染色体条带图（多特征热图）=============
# 对每个染色体，画各特征的归一化"丰度"
fig, ax = plt.subplots(figsize=(14, 8))

# 数据矩阵: 特征 x 染色体
features_for_heatmap = [
    ('LTR/Gypsy', ltr_gypsy, '#E74C3C'),
    ('LTR/Copia', ltr_copia, '#3498DB'),
    ('LTR/Cassandra', ltr_cassandra, '#F39C12'),
    ('DNA/PIF-Harbinger', dna_pif, '#9B59B6'),
    ('DNA/MULE-MuDR', dna_mule, '#1ABC9C'),
    ('DNA/CMC-EnSpm', dna_cmc, '#16A085'),
    ('DNA/hAT-Ac', dna_hat, '#2980B9'),
    ('LINE/L1', line_l1, '#34495E'),
    ('RC/Helitron', rc_helitron, '#E67E22'),
    ('GH28 (PG)', gh28, '#C0392B'),
    ('CE8 (PME)', ce8, '#D35400'),
    ('PL1', pl1, '#8E44AD'),
    ('GH17', gh17, '#27AE60'),
    ('LTR/Gypsy recent', ltr_recent, '#FF6B6B'),
]

# 归一化（按行归一化至 0-1）
matrix = []
labels = []
colors_bar = []
for label, d, c in features_for_heatmap:
    vals = [d.get(ch, 0) for ch in chroms]
    max_v = max(vals) if max(vals) > 0 else 1
    norm_vals = [v/max_v for v in vals]
    matrix.append(norm_vals)
    labels.append(label)
    colors_bar.append(c)

matrix = np.array(matrix)

# 横向堆叠柱状图（feature 在 y，染色体分布归一化）
y_pos = np.arange(len(labels))  # 14 个特征
bar_height = 0.7
matrix_t = matrix.T  # 转置: shape (n_chroms=8, n_features=14)

ax = fig.add_subplot(111)
left = np.zeros(len(labels))  # 每个特征的 left 起点，长度=14
for i, ch in enumerate(chroms):
    vals = matrix_t[i]  # 该染色体上各特征的归一化值，长度=14
    ax.barh(y_pos, vals, height=bar_height, left=left,
            color=plt.cm.tab20(i / 8), label=ch, edgecolor='white', linewidth=0.5)
    # 标注最大值
    max_v_idx = np.argmax(vals)
    if vals[max_v_idx] > 0.5:
        ax.text(left[max_v_idx] + vals[max_v_idx]/2, y_pos[max_v_idx],
                ch, ha='center', va='center', fontsize=8, weight='bold', color='white')
    left = left + vals

ax.set_yticks(y_pos)
ax.set_yticklabels(labels, fontsize=10)
ax.set_xlabel('Normalized proportion (each row scaled to max=1.0)', fontsize=11)
ax.set_title('Chromosomal distribution of key genomic features in Zuili plum', fontsize=13, weight='bold')
ax.set_xlim(0, 8.2)
ax.grid(axis='x', alpha=0.3)
ax.legend(title='Chromosome', loc='lower right', ncol=4, fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(img_dir, 'Fig24_chromosome_heatmap.png'), dpi=120, bbox_inches='tight')
plt.close()
print('Saved Fig24')

# ============= 图5: 染色体基因密度 + LTR 密度综合图 =============
fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

# 上：基因密度
gene_density = {}
for row in data:
    if row['特征名称'] == '蛋白编码基因':
        ch = row['染色体']
        if ch not in chr_len:
            continue
        n = int(row['数量'].replace(',', ''))
        gene_density[ch] = n / (chr_len[ch] / 1e6)

axes[0].bar(chroms, [gene_density.get(c, 0) for c in chroms], color='#3498DB', alpha=0.8)
axes[0].set_ylabel('Gene density (gene/Mb)', fontsize=11)
axes[0].set_title('Gene density and LTR/Gypsy density across chromosomes in Zuili plum', fontsize=13, weight='bold')
axes[0].grid(axis='y', alpha=0.3)
for i, c in enumerate(chroms):
    axes[0].text(i, gene_density.get(c, 0), f'{gene_density.get(c, 0):.1f}',
                 ha='center', va='bottom', fontsize=9)

# 下：LTR 密度（含 Gypsy、Copia、Cassandra、近期爆发）
ltr_density_gypsy = [ltr_gypsy.get(c, 0) / (chr_len[c] / 1e6) for c in chroms]
ltr_density_recent = [ltr_recent.get(c, 0) / (chr_len[c] / 1e6) for c in chroms]
cass_density = [ltr_cassandra.get(c, 0) / (chr_len[c] / 1e6) for c in chroms]

axes[1].bar(chroms, ltr_density_gypsy, color='#E74C3C', alpha=0.7, label='LTR/Gypsy')
axes[1].bar(chroms, cass_density, color='#F39C12', alpha=0.7, label='LTR/Cassandra')
axes[1].bar(chroms, ltr_density_recent, color='#3498DB', alpha=0.7, label='LTR/Gypsy recent (PercDiv<5%)')
axes[1].set_ylabel('LTR density (elements/Mb)', fontsize=11)
axes[1].set_xlabel('Chromosome', fontsize=11)
axes[1].legend(loc='upper right', fontsize=10)
axes[1].grid(axis='y', alpha=0.3)

# 标注 Chr3 Cassandra 异常
axes[1].annotate(f'Chr3 Cassandra\npeak: {cass_density[2]:.1f}/Mb',
                 xy=(2, cass_density[2]), xytext=(2, max(cass_density)*1.3),
                 arrowprops=dict(arrowstyle='->', color='darkorange'),
                 ha='center', fontsize=9, color='darkorange', weight='bold')

plt.tight_layout()
plt.savefig(os.path.join(img_dir, 'Fig25_gene_LTR_density.png'), dpi=120, bbox_inches='tight')
plt.close()
print('Saved Fig25')

# ============= 图6: GH28 串联簇染色体位置示意 =============
# 显示每个簇在染色体上的物理位置
fig, ax = plt.subplots(figsize=(14, 6))

# GH28 簇数据
clusters_data = [
    ('Chr1', 3498315, 3521674, 3, 'MY000424-426'),
    ('Chr1', 30024300, 30159881, 3, 'MY002360/2378-2379'),
    ('Chr2', 13867719, 13981946, 5, 'MY006181-185'),
    ('Chr3', 33528053, 33532903, 3, 'MY011268-1269'),
    ('Chr4', 12581327, 12684919, 7, 'MY012333-237/343'),
    ('Chr7', 24747584, 24828544, 7, 'MY022445-247/249/451-452'),
]

colors_cl = ['#E74C3C', '#C0392B', '#9B59B6', '#3498DB', '#1ABC9C', '#F39C12']

# 按染色体分行
chrom_y = {ch: i for i, ch in enumerate(chroms)}
ax.set_yticks(range(len(chroms)))
ax.set_yticklabels(chroms)
ax.set_ylim(-0.5, len(chroms) - 0.5)

# 画染色体主干
for ch in chroms:
    y = chrom_y[ch]
    ax.barh(y, chr_len[ch] / 1e6, height=0.15, color='lightgray', edgecolor='gray', alpha=0.6)
    # 标 centromere 估算位置（约 1/3 - 1/2 处）
    ax.axvline(chr_len[ch] / 1e6 * 0.45, ymin=(y-0.075)/len(chroms),
               ymax=(y+0.075)/len(chroms), color='gray', linestyle='--', linewidth=0.5)

# 画 GH28 簇
for (ch, st, en, n, gene_ids), color in zip(clusters_data, colors_cl):
    y = chrom_y[ch]
    # 簇跨度（Mb）
    st_mb = st / 1e6
    en_mb = en / 1e6
    span_mb = en_mb - st_mb
    # 簇条
    ax.barh(y, span_mb, left=st_mb, height=0.25, color=color, edgecolor='black', linewidth=0.5)
    # 标注
    ax.text((st_mb + en_mb) / 2, y + 0.18, f'{n} genes\n{gene_ids}',
            ha='center', va='bottom', fontsize=8, color=color, weight='bold')

ax.set_xlabel('Position on chromosome (Mb)', fontsize=12)
ax.set_title('Physical positions of GH28 tandem clusters in Zuili plum genome', fontsize=13, weight='bold')
ax.set_xlim(0, max(chr_len.values()) / 1e6 + 5)
ax.grid(axis='x', alpha=0.3)

# 图例
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='lightgray', edgecolor='gray', label='Chromosome backbone'),
    plt.Line2D([0], [0], color='gray', linestyle='--', label='Approx. centromere'),
    Patch(facecolor='#E74C3C', edgecolor='black', label='GH28 cluster (3-7 genes)'),
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=9)

plt.tight_layout()
plt.savefig(os.path.join(img_dir, 'Fig26_GH28_clusters_map.png'), dpi=120, bbox_inches='tight')
plt.close()
print('Saved Fig26')

print('\n所有图已生成完毕')