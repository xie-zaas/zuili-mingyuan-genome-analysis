# Zuili plum cv. 'Mingyuan' genome analysis code

Analysis and plotting scripts for the chromosome-level reference genome of the
Zuili plum (*Prunus salicina* Lindl. var. *Zuili* cv. 'Mingyuan').

## Citation

If you use these scripts, please cite the accompanying manuscript. The code is
archived on Zenodo; the concept DOI is given in the manuscript's Code availability
statement.

## What is here

| Script | Purpose |
|---|---|
| `analyze_aat2.py` | Amino-acid-transposition / repeat-family composition summaries |
| `analyze_family.py` | Repeat-family composition tables |
| `analyze_ltr.py` | LTR retrotransposon PercDiv-based insertion-age estimation |
| `analyze_ltr_full.py` | Full LTR analysis including subfamily assignment |
| `analyze_others.py` | Auxiliary repeat-class analyses |
| `analyze_others2.py` | Auxiliary repeat-class analyses (part 2) |
| `analyze_positions.py` | Physical positions of GH28 and other feature clusters |
| `plot_ltr.py` | LTR PercDiv distribution and insertion-age figures |
| `plot_ltr_copia.py` | LTR/Copia-specific distribution figure |
| `plot_positions.py` | Cluster-position map figure (Fig. 10) |
| `replot_fig22.py` | Regeneration of the cell-wall enzyme family distribution figure (Fig. 9) |
| `read_xls.py` | Parser for the vendor-supplied Excel-derived tables |
| `positions_to_docx.py` | Exports the cluster position table to Word |

## Environment

- Python 3.13
- Key packages: `pandas`, `numpy`, `matplotlib`, `python-docx`
- Reference R version for the enrichment analyses: R 4.1.0 with `clusterProfiler` 3.18.1

Install the core dependencies with:

```bash
pip install pandas numpy matplotlib python-docx
```

## Inputs

The scripts read the intermediate summary tables produced by the assembly and
annotation pipeline. Those tables are deposited with the supplementary data at
Figshare (see the accompanying manuscript for the DOI) and are also archived in the
data deposit described in the manuscript.

## Software used for the primary analyses

Assembly and annotation were performed with NextDenovo v2.5.0, Racon v1.4.3,
Medaka v1.4.4, Pilon v1.23, 3D-DNA v180922, Juicebox v1.11.08, RepeatMasker v4.1.2,
LTR_FINDER v1.0.7, MITE-Hunter v1.0, RepeatScout v1.0.6, TRF v4.0.9, AUGUSTUS v3.5.0,
GeMoMa v1.7.1, Trinity v2.15.2, PASA v1.3.3, EVidenceModeler v1.1.1, BLAST+ v2.10.0,
InterProScan v5.34, tRNAscan-SE v2.0.7, RNAmmer v1.2, INFERNAL v1.1.3,
OrthoFinder v2.4.0, MAFFT v7.471, RAxML v8.2.12, MCScanX and BUSCO v5.3.2.

The comparative-genomics step (OrthoFinder, MAFFT/RAxML phylogeny, MCScanX
collinearity) was carried out by Shanghai Ling'en Biotechnology Co., Ltd.
(上海凌恩生物科技有限公司); the scripts used for that step are not distributed here.

## License

MIT — see `LICENSE`.
