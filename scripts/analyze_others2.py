#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""从 KEGG category 提取 LOX/AAT/TMT/SWEET/SUT 等家族"""
import os
from collections import Counter, defaultdict
import re

cat_file = r"D:\xie's files\Xie's Data\李子数据\李子基因组、转录组测序信息\槜李-名媛\Result\04.Genome_Annotation\MY\KEGG\kegg.category.xls"

# 关键 KO 与功能描述关键词
target_keywords = {
    'LOX': ['lipoxygenase', 'LOX', 'linoleate'],
    'AAT': ['alcohol acyltransferase', 'AAT', 'BADH', 'BAHD acyltransferase', 'benzyl alcohol O-acetyltransferase'],
    'TMT': ['tonoplast monosaccharide transporter', 'TMT', 'MST'],
    'SWEET': ['SWEET', 'bidirectional sugar transporter'],
    'SUT': ['sucrose transporter', 'SUT', 'SUC'],
    'INV': ['invertase', 'INV', 'beta-fructofuranosidase'],
    'SUS': ['sucrose synthase', 'SUS'],
    'SPS': ['sucrose-phosphate synthase', 'SPS'],
    'ACO': ['ACC oxidase', 'ACO'],
    'ACS': ['ACC synthase', 'ACS'],
    'ALDH': ['aldehyde dehydrogenase', 'ALDH'],
    'F3H': ['flavanone 3-hydroxylase', 'F3H'],
    'CHS': ['chalcone synthase', 'CHS'],
    'ANR': ['anthocyanidin reductase', 'ANR'],
    'LAR': ['leucoanthocyanidin reductase', 'LAR'],
    'UFGT': ['UDP-glucose:flavonoid', 'UFGT', 'anthocyanidin 3-O-glucosyltransferase'],
    'MYB': [],  # 转录因子，不在 KEGG
    'bHLH': [],
    'WRKY': [],
    'NBS-LRR': [],  # 抗病基因
}

# 基因-功能映射
gene_func = defaultdict(set)

with open(cat_file, 'r', encoding='utf-8', errors='ignore') as f:
    next(f)  # 跳过表头
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 7:
            continue
        gene, entry_id, cat1, cat2, cat3, pathway, desc = parts
        desc_lower = desc.lower()

        for key, kws in target_keywords.items():
            if not kws:
                continue
            for kw in kws:
                if kw.lower() in desc_lower:
                    gene_func[key].add(gene)
                    break

print("=" * 70)
print("KEGG 注释中关键基因家族统计")
print("=" * 70)
for key in target_keywords:
    n = len(gene_func.get(key, set()))
    print(f"  {key:<12s}: {n:>4d} 个基因")

# 输出详细
print()
print("=" * 70)
print("LOX (脂氧合酶) 基因列表（涉及香气物质合成）")
print("=" * 70)
with open(cat_file, 'r', encoding='utf-8', errors='ignore') as f:
    next(f)
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 7:
            continue
        gene, entry_id, cat1, cat2, cat3, pathway, desc = parts
        if 'lipoxygenase' in desc.lower() or 'LOX' in desc:
            print(f"  {gene}  {entry_id}  {desc[:100]}")

print()
print("=" * 70)
print("AAT (醇酰基转移酶/BAHD) 基因列表（涉及香气酯类合成）")
print("=" * 70)
with open(cat_file, 'r', encoding='utf-8', errors='ignore') as f:
    next(f)
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 7:
            continue
        gene, entry_id, cat1, cat2, cat3, pathway, desc = parts
        if any(k.lower() in desc.lower() for k in ['alcohol acyltransferase', 'AAT', 'BAHD']):
            print(f"  {gene}  {entry_id}  {desc[:100]}")

print()
print("=" * 70)
print("TMT (液泡单糖转运蛋白) 基因")
print("=" * 70)
with open(cat_file, 'r', encoding='utf-8', errors='ignore') as f:
    next(f)
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 7:
            continue
        gene, entry_id, cat1, cat2, cat3, pathway, desc = parts
        if any(k.lower() in desc.lower() for k in ['tonoplast monosaccharide transporter', 'TMT']):
            print(f"  {gene}  {entry_id}  {desc[:120]}")

print()
print("=" * 70)
print("SWEET 基因（糖外排转运蛋白）")
print("=" * 70)
with open(cat_file, 'r', encoding='utf-8', errors='ignore') as f:
    next(f)
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 7:
            continue
        gene, entry_id, cat1, cat2, cat3, pathway, desc = parts
        if 'SWEET' in desc:
            print(f"  {gene}  {entry_id}  {desc[:120]}")

print()
print("=" * 70)
print("SUT/SUC（蔗糖转运蛋白）")
print("=" * 70)
with open(cat_file, 'r', encoding='utf-8', errors='ignore') as f:
    next(f)
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 7:
            continue
        gene, entry_id, cat1, cat2, cat3, pathway, desc = parts
        if 'sucrose transporter' in desc.lower() or 'SUT' in desc or 'SUC' in desc:
            print(f"  {gene}  {entry_id}  {desc[:120]}")

print()
print("=" * 70)
print("INV（转化酶）")
print("=" * 70)
with open(cat_file, 'r', encoding='utf-8', errors='ignore') as f:
    next(f)
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 7:
            continue
        gene, entry_id, cat1, cat2, cat3, pathway, desc = parts
        if 'invertase' in desc.lower():
            print(f"  {gene}  {entry_id}  {desc[:120]}")