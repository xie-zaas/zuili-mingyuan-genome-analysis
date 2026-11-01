#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""从CAZy详细注释统计 GH28、TMT2、LOX、AAT 等基因家族"""
import os

base = r"D:\xie's files\Xie's Data\李子数据\李子基因组、转录组测序信息\槜李-名媛\Result\04.Genome_Annotation\MY\CAZy\MY.anno.xls"

# 读全部内容
with open(base, 'r', encoding='utf-8', errors='ignore') as f:
    lines = f.readlines()

# 跳过表头
header = lines[0].strip().split('\t')
# 列: #Query, Length, Subject, Length, Q-start, Q-end, S-start, S-end, Evalue, Coverage, cazy-class, cazy-activities

# 统计 cazy-class
from collections import Counter, defaultdict
class_count = Counter()
class_genes = defaultdict(set)
class_subjects = defaultdict(Counter)
specific_targets = {
    'GH28': 'Polygalacturonase / 多聚半乳糖醛酸酶',
    'PL': 'Pectate lyase',
    'CE8': 'Pectin methylesterase / 果胶甲酯酶',
    'PME': 'Pectin methylesterase',
    'GT': 'Glycosyltransferase',
    'GH': 'Glycoside Hydrolase',
    'GH17': 'beta-1,3-glucanase / 葡聚糖酶',
    'LOX': 'Lipoxygenase (脂氧合酶)',
    'AAT': 'Alcohol acyltransferase (醇酰基转移酶)',
    'TMT': 'Tonoplast monosaccharide transporter (液泡糖转运)',
    'SUT': 'Sucrose transporter',
    'SWEET': 'Sugar transporter',
    'INV': 'Invertase / 转化酶',
    'SUS': 'Sucrose synthase / 蔗糖合酶',
    'SPS': 'Sucrose phosphate synthase',
    'ACO': 'ACC oxidase',
    'ACS': 'ACC synthase',
}

specific_gene_count = Counter()
specific_gene_genes = defaultdict(set)

for line in lines[1:]:
    parts = line.strip().split('\t')
    if len(parts) < 12:
        continue
    query = parts[0]
    subject = parts[2]
    cazy_class = parts[10]
    activity = parts[11]
    class_count[cazy_class] += 1
    class_genes[cazy_class].add(query)
    class_subjects[cazy_class][subject] += 1

    # 关键字匹配
    subj_upper = subject.upper()
    for key, desc in specific_targets.items():
        if key.upper() in subj_upper or (key.upper() in activity.upper() if activity else False):
            specific_gene_count[key] += 1
            specific_gene_genes[key].add(query)
            break

print("=" * 70)
print("CAZy 分类统计")
print("=" * 70)
for c, n in class_count.most_common():
    print(f"  {c:8s} 模块数: {n:5d}  对应基因数: {len(class_genes[c]):5d}")

print()
print("=" * 70)
print("重点关注的基因家族")
print("=" * 70)
for key, desc in specific_targets.items():
    n = specific_gene_count.get(key, 0)
    genes = specific_gene_genes.get(key, set())
    print(f"  {key:8s} ({desc:30s}): 模块数 {n:4d}  基因数 {len(genes):4d}")

print()
print("=" * 70)
print("CAZy 内 Subject 详细分类（每个家族类型取前10）")
print("=" * 70)
for c in ['GH', 'PL', 'CE', 'GT', 'AA', 'CBM']:
    print(f"\n  {c} 家族细分:")
    for subj, n in class_subjects[c].most_common(15):
        print(f"    {subj:8s}: {n:4d}")

# 关键GH家族列表
print()
print("=" * 70)
print("GH28 多聚半乳糖醛酸酶（PG）相关条目")
print("=" * 70)
gh28_genes = []
gh28_subjects = Counter()
for line in lines[1:]:
    parts = line.strip().split('\t')
    if len(parts) < 12:
        continue
    if 'GH28' in parts[2] or 'GH28' in parts[11]:
        gh28_genes.append((parts[0], parts[2], parts[11][:80]))
        gh28_subjects[parts[2]] += 1
print(f"  共 {len(gh28_genes)} 个 GH28 注释条目（按基因组装）")
for g, s, a in gh28_genes[:30]:
    print(f"    {g}\t{s}\t{a}")
print(f"\n  涉及 Subject 类型:")
for s, n in gh28_subjects.most_common():
    print(f"    {s}: {n}")

# PME / 果胶甲酯酶相关
print()
print("=" * 70)
print("CE8 果胶甲酯酶（PME）")
print("=" * 70)
ce8_genes = set()
ce8_subjects = Counter()
for line in lines[1:]:
    parts = line.strip().split('\t')
    if len(parts) < 12:
        continue
    if 'CE8' in parts[2] or 'pectin methylesterase' in (parts[11] if len(parts) > 11 else '').lower():
        ce8_genes.add(parts[0])
        ce8_subjects[parts[2]] += 1
print(f"  PME/CE8 相关基因: {len(ce8_genes)} 个")
for s, n in ce8_subjects.most_common():
    print(f"    {s}: {n}")

# PL 多糖裂解酶
print()
print("=" * 70)
print("PL 多糖裂解酶（果胶酸裂解酶、纤维素裂解酶等）")
print("=" * 70)
pl_genes = set()
pl_subjects = Counter()
for line in lines[1:]:
    parts = line.strip().split('\t')
    if len(parts) < 12:
        continue
    if parts[10] == 'PL':
        pl_genes.add(parts[0])
        pl_subjects[parts[2]] += 1
print(f"  PL 类基因: {len(pl_genes)} 个")
for s, n in pl_subjects.most_common():
    print(f"    {s}: {n}")