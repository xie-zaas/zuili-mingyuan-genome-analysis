#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""从repeat注释GFF中提取LTR-RT详细分类（Gypsy/Copia/CACTA等）"""
import os
from collections import Counter, defaultdict

gff = r"D:\xie's files\Xie's Data\李子数据\李子基因组、转录组测序信息\槜李-名媛\Result\03.Genome_Component\Repeat\MY.repbase.out.gff"

# RepeatMasker .out.gff 格式（实际是 .out.txt）通常：
# 列1: 染色体  列2: 源  列3: 类型(class/family)  列4: start  列5: end  列6: score  列7: strand  列8: phase  列9: attributes

# 先看看前30行
print("=== GFF前30行 ===")
with open(gff, 'r', encoding='utf-8', errors='ignore') as f:
    head_lines = [next(f) for _ in range(30)]
for i, l in enumerate(head_lines):
    print(f"  {i+1}: {l.rstrip()[:200]}")