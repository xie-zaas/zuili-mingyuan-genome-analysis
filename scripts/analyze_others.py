#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""从KEGG注释中查找 LOX/AAT/TMT/SWEET/SUT 等基因家族"""
import os
from collections import Counter, defaultdict

ko_file = r"D:\xie's files\Xie's Data\李子数据\李子基因组、转录组测序信息\槜李-名媛\Result\04.Genome_Annotation\MY\KEGG\ko_table.xls"
cat_file = r"D:\xie's files\Xie's Data\李子数据\李子基因组、转录组测序信息\槜李-名媛\Result\04.Genome_Annotation\MY\KEGG\kegg.category.xls"

# 查看 ko_table 列
print("=== ko_table.xls 前 5 行 ===")
with open(ko_file, 'r', encoding='utf-8', errors='ignore') as f:
    for i, l in enumerate(f):
        if i > 5:
            break
        print(f"  {l.rstrip()[:300]}")

print()
print("=== kegg.category.xls 前 5 行 ===")
with open(cat_file, 'r', encoding='utf-8', errors='ignore') as f:
    for i, l in enumerate(f):
        if i > 8:
            break
        print(f"  {l.rstrip()[:300]}")