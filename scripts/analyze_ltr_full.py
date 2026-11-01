#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""完整分析 LTR-RT 详细分类、Gypsy/Copia 比例、插入时间估算"""
import os
import re
from collections import Counter, defaultdict

gff = r"D:\xie's files\Xie's Data\李子数据\李子基因组、转录组测序信息\槜李-名媛\Result\03.Genome_Component\Repeat\MY.repbase.out.gff"

# 统计
class_count = Counter()
class_length = defaultdict(int)
class_percdiv = defaultdict(list)  # 用于估算插入时间

# 简易处理大文件
def parse_attrs(attr_str):
    """解析 GFF 第9列的属性"""
    res = {}
    for kv in attr_str.split(';'):
        if not kv:
            continue
        if '=' in kv:
            k, v = kv.split('=', 1)
            res[k] = v
    return res

# LTR 分类子项
ltr_subclass_count = Counter()
ltr_subclass_length = defaultdict(int)
ltr_percdiv_intact = []  # 用于估算整段 LTR 的插入时间
gypsy_count = 0
copia_count = 0
gypsy_len = 0
copia_len = 0

with open(gff, 'r', encoding='utf-8', errors='ignore') as f:
    for line in f:
        if line.startswith('#'):
            continue
        if not line.strip():
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9:
            continue
        chrom = parts[0]
        cls = parts[2]  # 这里其实是 type（如 LTR/Gypsy）
        start = int(parts[3])
        end = int(parts[4])
        length = end - start + 1
        attrs = parse_attrs(parts[8])
        cls_class = attrs.get('Class', 'Unknown')
        perc_div_str = attrs.get('PercDiv', '0.0')
        try:
            perc_div = float(perc_div_str)
        except:
            perc_div = 0.0

        # 一级分类统计（用 Class 字段）
        class_count[cls_class] += 1
        class_length[cls_class] += length
        class_percdiv[cls_class].append(perc_div)

        # 二级 LTR 子分类
        if cls_class.startswith('LTR/'):
            subclass = cls_class[4:]  # 去掉 'LTR/'
            ltr_subclass_count[subclass] += 1
            ltr_subclass_length[subclass] += length

# 输出结果
genome_size = 294_811_541

print("=" * 70)
print("Class 分类（一级）")
print("=" * 70)
print(f"{'Class':<30s}{'N':>8s}{'Total Len':>15s}{'%Genome':>12s}{'Avg%Div':>10s}")
for cls, cnt in class_count.most_common():
    pct = class_length[cls] / genome_size * 100
    avg_div = sum(class_percdiv[cls]) / len(class_percdiv[cls]) if class_percdiv[cls] else 0
    print(f"{cls:<30s}{cnt:>8d}{class_length[cls]:>15d}{pct:>11.2f}%{avg_div:>9.2f}%")

print()
print("=" * 70)
print("LTR 子分类（详细）")
print("=" * 70)
ltr_total_count = sum(ltr_subclass_count.values())
ltr_total_len = sum(ltr_subclass_length.values())
print(f"{'SubClass':<25s}{'N':>8s}{'Total Len':>15s}{'%Genome':>12s}{'%of LTR':>10s}")
for sub, cnt in ltr_subclass_count.most_common():
    pct = class_length[f'LTR/{sub}'] / genome_size * 100
    pct_ltr = ltr_subclass_length[sub] / ltr_total_len * 100 if ltr_total_len else 0
    print(f"  LTR/{sub:<20s}{cnt:>8d}{ltr_subclass_length[sub]:>15d}{pct:>11.2f}%{pct_ltr:>9.2f}%")

print(f"\nLTR总数: {ltr_total_count}, LTR总长: {ltr_total_len} bp, 占基因组: {ltr_total_len/genome_size*100:.2f}%")

# Gypsy vs Copia
print()
print("=" * 70)
print("LTR/Gypsy vs LTR/Copia 比较")
print("=" * 70)
gypsy_n = ltr_subclass_count.get('Gypsy', 0)
copia_n = ltr_subclass_count.get('Copia', 0)
gypsy_len_v = ltr_subclass_length.get('Gypsy', 0)
copia_len_v = ltr_subclass_length.get('Copia', 0)
print(f"  LTR/Gypsy: {gypsy_n} 个, {gypsy_len_v} bp ({gypsy_len_v/genome_size*100:.2f}%)")
print(f"  LTR/Copia: {copia_n} 个, {copia_len_v} bp ({copia_len_v/genome_size*100:.2f}%)")
print(f"  Gypsy:Copia 数量比: {gypsy_n/max(copia_n,1):.2f}:1")
print(f"  Gypsy:Copia 长度比: {gypsy_len_v/max(copia_len_v,1):.2f}:1")
total_known = gypsy_len_v + copia_len_v
print(f"  Gypsy 占已知 Gypsy+Copia: {gypsy_len_v/total_known*100:.2f}%")
print(f"  Copia 占已知 Gypsy+Copia: {copia_len_v/total_known*100:.2f}%")

# 简单估算 Gypsy/Copia 的插入时间分布
print()
print("=" * 70)
print("LTR 插入时间分布（PercDiv 反映分歧度）")
print("=" * 70)
# 公式：T = PercDiv / (2 * mu * 100)
# mu 植物 LTR 突变率约 1.3e-8 ~ 1.0e-8 / 每年 (Wicker et al.)
# 采用 1.3e-9 substitutions/site/year
mu = 1.3e-9  # substitutions/site/year

# 按 PercDiv 区间统计
bins = [0, 5, 10, 15, 20, 25, 30, 50, 100]
bin_gypsy = [0]*len(bins)
bin_copia = [0]*len(bins)

with open(gff, 'r', encoding='utf-8', errors='ignore') as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9:
            continue
        attrs = parse_attrs(parts[8])
        cls_class = attrs.get('Class', '')
        try:
            perc_div = float(attrs.get('PercDiv', '0'))
        except:
            perc_div = 0
        if cls_class == 'LTR/Gypsy':
            for i in range(len(bins)-1, -1, -1):
                if perc_div >= bins[i]:
                    bin_gypsy[i] += 1
                    break
        elif cls_class == 'LTR/Copia':
            for i in range(len(bins)-1, -1, -1):
                if perc_div >= bins[i]:
                    bin_copia[i] += 1
                    break

print(f"{'PercDiv%':<12s}{'Gypsy':>10s}{'Copia':>10s}{'Gypsy时代(My)':>18s}")
for i in range(len(bins)-1):
    t_myr = bins[i+1] / (2 * mu * 100) / 1e6  # 换算成百万年
    print(f"  {bins[i]}–{bins[i+1]}%  {bin_gypsy[i]:>10d}{bin_copia[i]:>10d}{t_myr:>14.2f} Mya")

# 计算平均分歧度
gypsy_divs = []
copia_divs = []
with open(gff, 'r', encoding='utf-8', errors='ignore') as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9:
            continue
        attrs = parse_attrs(parts[8])
        cls_class = attrs.get('Class', '')
        try:
            perc_div = float(attrs.get('PercDiv', '0'))
        except:
            perc_div = 0
        if cls_class == 'LTR/Gypsy':
            gypsy_divs.append(perc_div)
        elif cls_class == 'LTR/Copia':
            copia_divs.append(perc_div)

if gypsy_divs:
    print(f"\nGypsy 平均 PercDiv: {sum(gypsy_divs)/len(gypsy_divs):.2f}%, 推算平均插入年龄: {sum(gypsy_divs)/len(gypsy_divs)/(2*mu*100)/1e6:.2f} Mya")
if copia_divs:
    print(f"Copia 平均 PercDiv: {sum(copia_divs)/len(copia_divs):.2f}%, 推算平均插入年龄: {sum(copia_divs)/len(copia_divs)/(2*mu*100)/1e6:.2f} Mya")