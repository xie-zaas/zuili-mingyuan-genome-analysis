#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""读取Result文件夹中的xls文件"""
import pandas as pd
import os

base_dir = r"D:\xie's files\Xie's Data\李子数据\李子基因组、转录组测序信息\槜李-名媛\Result"

files = [
    ('Repeat', r'03.Genome_Component\Repeat\MY.repbase.stat.xls'),
    ('TRF',    r'03.Genome_Component\Repeat\MY.trf.stat.xls'),
    ('Gene',   r'03.Genome_Component\Gene\MY.gene.stat.xls'),
    ('ncRNA',  r'03.Genome_Component\ncRNA\MY.ncRNA.stat.xls'),
    ('CAZy_full',  r'04.Genome_Annotation\MY\CAZy\MY.anno.st.xls'),
    ('CAZy_class', r'04.Genome_Annotation\MY\CAZy\MY.anno.xls'),
    ('Kmer',   r'02.Assembly\K-mer.xls'),
    ('Assembly', r'02.Assembly\assemblyresult.xls'),
    ('Chr.stat', r'02.Assembly\MY.Chr.stat.xls'),
    ('annotation', r'04.Genome_Annotation\annotation.stat.xls'),
    ('all_anno', r'04.Genome_Annotation\MY\all_anno_st.xls'),
]

for name, rel in files:
    p = os.path.join(base_dir, rel)
    if not os.path.exists(p):
        print(f'\n===== {name}: NOT EXIST =====')
        continue
    print(f'\n===== {name} =====')
    print(f'  File: {rel}')
    try:
        # 旧版 .xls 用 xlrd
        df = pd.read_excel(p, engine='xlrd')
        print(f'  Shape: {df.shape}')
        print(f'  Columns: {list(df.columns)}')
        # 全部打印（小文件）
        with pd.option_context('display.max_rows', None, 'display.max_columns', None, 'display.width', 200, 'display.max_colwidth', 50):
            print(df.to_string(index=False))
    except Exception as e:
        print(f'  Error: {e}')
        # 尝试用纯文本读取
        with open(p, 'rb') as f:
            data = f.read()
        print(f'  Raw bytes preview (first 500):')
        print(data[:500])