#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""查找醇酰基转移酶 (AAT, EC 2.3.1.84) - 真正的水果香气合成酶"""
import os
from collections import defaultdict, Counter

cat_file = r"D:\xie's files\Xie's Data\李子数据\李子基因组、转录组测序信息\槜李-名媛\Result\04.Genome_Annotation\MY\KEGG\kegg.category.xls"

# 查所有 acetyltransferase / acyltransferase 相关条目
print("=" * 70)
print("所有酰基转移酶 (涉及香气酯类合成)")
print("=" * 70)
genes = defaultdict(set)
with open(cat_file, 'r', encoding='utf-8', errors='ignore') as f:
    next(f)
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 7:
            continue
        gene, entry_id, cat1, cat2, cat3, pathway, desc = parts
        dlow = desc.lower()
        # AAT EC 2.3.1.84: alcohol O-acetyltransferase
        if 'acetyltransferase' in dlow or 'acyltransferase' in dlow:
            if 'alcohol' in dlow or 'O-acetyl' in dlow or 'O-acyl' in dlow or 'BAHD' in dlow or 'acetyl-CoA' in dlow:
                genes[desc[:80]].add(gene)

for desc, gs in sorted(genes.items(), key=lambda x: -len(x[1])):
    print(f"  {desc}: {len(gs)} 个基因")
    if len(gs) <= 8:
        for g in gs:
            print(f"    {g}")

# 香气物质代谢相关
print()
print("=" * 70)
print("香气物质前体合成酶（萜类、酯类、醛类）")
print("=" * 70)
keywords = {
    'AAT 真正醇酰基转移酶 (EC 2.3.1.84)': ['alcohol O-acetyltransferase', 'alcohol acetyltransferase', 'AAT;'],
    'BAHD 酰基转移酶': ['BAHD'],
    'TPS 萜类合成酶': ['terpene synthase', 'TPS'],
    '苯丙烷类 / 苯甲醇苯甲醛相关': ['phenylpropanoid', 'benzaldehyde', 'benzyl alcohol'],
    'ADH 乙醇脱氢酶': ['alcohol dehydrogenase', 'ADH;'],
    'AOS 丙二烯氧化物合酶': ['allene oxide synthase', 'AOS;'],
    'HPL 脂氢过氧化物裂解酶': ['hydroperoxide lyase', 'HPL;'],
    '萜类主链合成': ['farnesyl', 'geranyl', 'mevalonate'],
}

with open(cat_file, 'r', encoding='utf-8', errors='ignore') as f:
    next(f)
    rows = list(f)

for label, kws in keywords.items():
    found_genes = set()
    found_descs = Counter()
    for line in rows:
        parts = line.strip().split('\t')
        if len(parts) < 7:
            continue
        gene, entry_id, cat1, cat2, cat3, pathway, desc = parts
        dlow = desc.lower()
        for kw in kws:
            if kw.lower() in dlow:
                found_genes.add(gene)
                found_descs[desc[:80]] += 1
                break
    print(f"\n  {label}: {len(found_genes)} 个基因")
    for d, n in found_descs.most_common(8):
        print(f"    {n}x {d}")