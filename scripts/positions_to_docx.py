#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""将染色体特殊区段位置信息表 + 6 张分布图保存为 Word 文档（投稿补充材料）"""
import os
from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
import csv


def set_cell_font(run, font_name='宋体', size=10):
    run.font.name = font_name
    run._element.rPr.rFonts.set(qn('w:eastAsia'), font_name)
    run.font.size = Pt(size)


def add_table_from_csv(doc, csv_path, max_rows=None):
    """从 CSV 创建Word表格"""
    with open(csv_path, 'r', encoding='utf-8-sig') as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = list(reader)

    if max_rows:
        rows = rows[:max_rows]

    n_rows = len(rows) + 1
    n_cols = len(header)
    table = doc.add_table(rows=n_rows, cols=n_cols, style='Table Grid')
    table.alignment = WD_ALIGN_PARAGRAPH.CENTER

    # 表头
    for j, h in enumerate(header):
        cell = table.rows[0].cells[j]
        cell.text = ''
        p = cell.paragraphs[0]
        run = p.add_run(h)
        set_cell_font(run, '宋体', 9)
        run.bold = True
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        # 表头背景色
        shading = OxmlElement('w:shd')
        shading.set(qn('w:fill'), 'B4C7E7')
        cell._tc.get_or_add_tcPr().append(shading)

    # 数据行
    for i, row_data in enumerate(rows):
        for j, val in enumerate(row_data):
            if j < n_cols:
                cell = table.rows[i+1].cells[j]
                cell.text = ''
                p = cell.paragraphs[0]
                run = p.add_run(str(val))
                set_cell_font(run, '宋体', 9)
                p.alignment = WD_ALIGN_PARAGRAPH.LEFT if j == 1 else WD_ALIGN_PARAGRAPH.CENTER


def add_md_table(doc, md_text):
    """将 Markdown 表格插入 Word"""
    lines = md_text.strip().split('\n')
    table_lines = [l for l in lines if l.strip().startswith('|')]
    if not table_lines:
        return
    # 跳过分隔行
    clean = []
    for l in table_lines:
        if l.strip().startswith('|---') or '---' in l and '|' in l:
            continue
        cells = [c.strip() for c in l.strip().strip('|').split('|')]
        clean.append(cells)

    n_rows = len(clean)
    n_cols = max(len(r) for r in clean)
    table = doc.add_table(rows=n_rows, cols=n_cols, style='Table Grid')
    table.alignment = WD_ALIGN_PARAGRAPH.CENTER

    for i, row_data in enumerate(clean):
        for j in range(n_cols):
            cell = table.rows[i].cells[j]
            cell.text = ''
            p = cell.paragraphs[0]
            val = row_data[j] if j < len(row_data) else ''
            run = p.add_run(val)
            set_cell_font(run, '宋体', 9)
            if i == 0:
                run.bold = True
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER
                shading = OxmlElement('w:shd')
                shading.set(qn('w:fill'), 'B4C7E7')
                cell._tc.get_or_add_tcPr().append(shading)
            else:
                p.alignment = WD_ALIGN_PARAGRAPH.LEFT if j == 1 else WD_ALIGN_PARAGRAPH.CENTER


def add_md_content(doc, md_text):
    """把 Markdown 内容转 Word 段落（标题、粗体、表格）"""
    lines = md_text.split('\n')
    i = 0
    table_buffer = []

    def flush_table():
        nonlocal table_buffer
        if table_buffer:
            add_md_table(doc, '\n'.join(table_buffer))
            table_buffer = []

    while i < len(lines):
        line = lines[i].rstrip('\n')

        if not line.strip():
            flush_table()
            i += 1
            continue

        if line.startswith('# '):
            flush_table()
            p = doc.add_paragraph()
            run = p.add_run(line[2:].strip())
            run.font.name = '黑体'
            run._element.rPr.rFonts.set(qn('w:eastAsia'), '黑体')
            run.font.size = Pt(18)
            run.bold = True
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            i += 1
            continue
        if line.startswith('## '):
            flush_table()
            p = doc.add_paragraph()
            run = p.add_run(line[3:].strip())
            run.font.name = '黑体'
            run._element.rPr.rFonts.set(qn('w:eastAsia'), '黑体')
            run.font.size = Pt(14)
            run.bold = True
            i += 1
            continue
        if line.startswith('### '):
            flush_table()
            p = doc.add_paragraph()
            run = p.add_run(line[4:].strip())
            run.font.name = '黑体'
            run._element.rPr.rFonts.set(qn('w:eastAsia'), '黑体')
            run.font.size = Pt(12)
            run.bold = True
            i += 1
            continue

        if line.strip().startswith('|'):
            table_buffer.append(line)
            i += 1
            continue
        else:
            flush_table()

        if line.strip() == '---':
            i += 1
            continue

        # 普通段落 (处理粗体)
        text = line.strip()
        if text:
            p = doc.add_paragraph()
            parts = text.split('**')
            for idx, part in enumerate(parts):
                if not part:
                    continue
                run = p.add_run(part)
                set_cell_font(run, '宋体', 10.5)
                if idx % 2 == 1:  # 粗体部分
                    run.bold = True
            p.paragraph_format.first_line_indent = Cm(0.74)

        i += 1

    flush_table()


def add_image(doc, img_path, width_cm=14):
    if not os.path.exists(img_path):
        return
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = p.add_run()
    run.add_picture(img_path, width=Cm(width_cm))


def add_page_break(doc):
    from docx.enum.text import WD_BREAK
    p = doc.add_paragraph()
    run = p.add_run()
    run.add_break(WD_BREAK.PAGE)


# ============================================
# 主程序
# ============================================
base = r"D:\xie's files\Xie's Data\李子数据\李子基因组、转录组测序信息\槜李-名媛"
img_dir = os.path.join(base, r"论文图表_v10")
out_dir = base

# 创建文档
doc = Document()

# 页面设置
section = doc.sections[0]
section.top_margin = Cm(2.0)
section.bottom_margin = Cm(2.0)
section.left_margin = Cm(2.0)
section.right_margin = Cm(2.0)

# 默认样式
style = doc.styles['Normal']
style.font.name = 'Times New Roman'
style.element.rPr.rFonts.set(qn('w:eastAsia'), '宋体')
style.font.size = Pt(10.5)

# ============= 标题 =============
p = doc.add_paragraph()
run = p.add_run('槜李\'名媛\'染色体特殊区段位置信息表')
run.font.name = '黑体'
run._element.rPr.rFonts.set(qn('w:eastAsia'), '黑体')
run.font.size = Pt(20)
run.bold = True
p.alignment = WD_ALIGN_PARAGRAPH.CENTER

p = doc.add_paragraph()
run = p.add_run('（投稿补充材料 / Supplementary Material）')
run.font.size = Pt(12)
run.italic = True
p.alignment = WD_ALIGN_PARAGRAPH.CENTER

p = doc.add_paragraph()
run = p.add_run('作者：谢小波    浙江省农业科学院 园艺研究所，杭州 310021\n邮箱：xiexb@zaas.ac.cn 或 xie.xiaobo@qq.com')
set_cell_font(run, '宋体', 11)
p.alignment = WD_ALIGN_PARAGRAPH.CENTER

doc.add_paragraph()  # 空行

# ============= 表格 =============
# 先添加 Markdown 表格内容
md_path = os.path.join(img_dir, '染色体特殊区段位置信息表.md')
with open(md_path, 'r', encoding='utf-8') as f:
    md_content = f.read()

# 删除第一个一级标题（已添加）
md_lines = md_content.split('\n')
filtered = []
skip_title = True
for line in md_lines:
    if skip_title and line.startswith('# '):
        skip_title = False
        continue
    filtered.append(line)
md_content = '\n'.join(filtered)

add_md_content(doc, md_content)

# ============= 完整 CSV 数据表 =============
add_page_break(doc)

p = doc.add_paragraph()
run = p.add_run('附录：染色体特殊区段位置信息完整数据表（CSV 内容）')
run.font.name = '黑体'
run._element.rPr.rFonts.set(qn('w:eastAsia'), '黑体')
run.font.size = Pt(16)
run.bold = True
p.alignment = WD_ALIGN_PARAGRAPH.LEFT

p = doc.add_paragraph()
run = p.add_run('（共 ' + str(sum(1 for _ in open(os.path.join(img_dir, "染色体特殊区段位置信息表.csv"), 'r', encoding='utf-8-sig')) - 1) + ' 条记录）')
run.italic = True
set_cell_font(run, '宋体', 10)

csv_path = os.path.join(img_dir, '染色体特殊区段位置信息表.csv')
add_table_from_csv(doc, csv_path)

# ============= 染色体分布图 =============
for i in range(21, 27):
    add_page_break(doc)
    fig_path = os.path.join(img_dir, f'Fig{i}_*.png')
    import glob
    matches = glob.glob(fig_path)
    if not matches:
        continue
    img_path = matches[0]

    titles_cn = {
        21: '图1 LTR/Gypsy 与 LTR/Copia 在8条染色体上的分布',
        22: '图2 关键细胞壁水解酶基因家族染色体分布',
        23: '图3 KEGG 关键基因家族（香气/糖代谢）染色体分布',
        24: '图4 染色体多特征分布归一化堆叠图',
        25: '图5 基因密度与 LTR/Gypsy 密度染色体比较',
        26: '图6 GH28 多聚半乳糖醛酸酶基因串联簇染色体物理位置',
    }
    p = doc.add_paragraph()
    run = p.add_run(titles_cn.get(i, f'图{i}'))
    run.font.name = '黑体'
    run._element.rPr.rFonts.set(qn('w:eastAsia'), '黑体')
    run.font.size = Pt(14)
    run.bold = True
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER

    add_image(doc, img_path, width_cm=15)

# ============= 文档说明 =============
add_page_break(doc)
p = doc.add_paragraph()
run = p.add_run('数据来源与处理说明')
run.font.name = '黑体'
run._element.rPr.rFonts.set(qn('w:eastAsia'), '黑体')
run.font.size = Pt(14)
run.bold = True

texts = [
    '1. LTR 注释 GFF：Result/03.Genome_Component/Repeat/MY.repbase.out.gff（包含 Class 字段和 PercDiv 字段）。',
    '2. 基因 GFF：Result/03.Genome_Component/Gene/MY.gff（提取每个基因的染色体与起止位置）。',
    '3. CAZy 详细注释：Result/04.Genome_Annotation/MY/CAZy/MY.anno.xls（按基因 ID 反向映射染色体位置）。',
    '4. KEGG 注释：Result/04.Genome_Annotation/MY/KEGG/kegg.category.xls（按 Entry 和 Description 匹配关键基因家族）。',
    '5. 染色体长度：Result/02.Assembly/MY.Chr.stat.xls。',
    '6. LTR 插入时间推算公式：T = PercDiv / (2 × μ × 100)，其中 μ = 1.3 × 10⁻⁹ substitutions/site/year (Wicker et al., 2003)。',
    '7. GH28 串联簇识别标准：染色体上位置相邻 ≤ 200 kb 且 ≥3 个 GH28 基因。',
    '8. KEGG 中 TMT2 (液泡糖转运蛋白) 和 SWEET 因 KEGG 数据库注释不全未直接命中，需后续通过 InterProScan/Pfam 注释补充；本文报告的 12 个 LOX、24 个 LOX、24 个 ADH、37 个 SUT 等均通过 KEGG 注释直接确认。',
    '',
    '主要发现：',
    '• GH28（多聚半乳糖醛酸酶）共 70 个，在 Chr2、Chr4、Chr7 上各含 14 个；Chr2 13.87–13.98 Mb 区间存在 5 个连续 GH28 基因（MY006181–MY006185），提示近期串联扩张。',
    '• CE8（果胶甲酯酶 PME）共 60 个（KEGG 直接命中 61 个），与 GH28 协同构成"果胶降解三联体"，是槜李果肉"化浆可吸食"的细胞学基础。',
    '• LTR/Gypsy 共 60 681 个，占基因组 25.10%；其中近期活跃（PercDiv<5%）元件 19 056 个，主要分布于 Chr3（3 411 个）和 Chr2（2 993 个），与 2 500 年无性繁殖历史吻合。',
    '• LTR/Cassandra 仅 1 786 个，但整体 PercDiv 仅 5.01%（平均插入年龄 ~1.93 Mya），属于近期活跃爆发的 LTR 家族；Chr3 上 Cassandra 数量（968）异常突出，提示 Chr3 为该家族扩张热点。',
    '• KEGG 关键家族中，LOX（脂氧合酶，12 个）主要分布于 Chr3（5 个）和 Chr4（3 个），与 C6/C9 醛类香气合成相关；ADH（24 个）在 Chr5 极度富集（14 个，70%）；TPS（萜类合成酶，17 个）以 Chr6（6 个）为主。',
    '',
    '本附表与论文正文 §3.6 性状关联分析互为补充，可作为审稿人复核相关结论的原始数据来源。',
]
for t in texts:
    p = doc.add_paragraph()
    run = p.add_run(t)
    set_cell_font(run, '宋体', 10.5)
    p.paragraph_format.first_line_indent = Cm(0.74)

# 保存
out_path = os.path.join(out_dir, '槜李_名媛_染色体特殊区段位置信息表.docx')
doc.save(out_path)
print(f'Saved: {out_path}')
print(f'  Size: {os.path.getsize(out_path)/1024:.1f} KB')