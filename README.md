<img align="left" src="./imgs/logo.png">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<img src="./imgs/jhu.png">

# ASCOT

[Jonathan P. Ling](https://scholar.google.com/citations?user=dGBD72YAAAAJ), [Chris Wilks](https://github.com/ChristopherWilks), [Rone Charles](https://github.com/ch4rr0), [Ben Langmead](http://www.langmead-lab.org/), [Seth Blackshaw](http://neuroscience.jhu.edu/research/faculty/7)

[ASCOT](http://ascot.cs.jhu.edu) quantifies alternative splicing and gene expression across tens of thousands of bulk and single-cell RNA-Seq datasets in human and mouse. This repository contains the scripts used to generate the PSI/NAUC tables for this resource.

ASCOT uses [annotation-free methods](http://www.biorxiv.org/) to detect exon percent spliced-in (PSI) values within [Snaptron](http://snaptron.cs.jhu.edu/), a rapidly queryable database of splice junction counts derived from public RNA-seq data. Gene expression levels are calculated using a normalized "area-under-curve" (NAUC) metric as described in [recount2](http://google.com).

Please refer to our [bioRxiv preprint](http://www.biorxiv.org/) and the accompanying [ASCOT website](http://ascot.cs.jhu.edu). All data tables are available for [download](http://ascot.cs.jhu.edu/data).

Comments and suggestions are always welcome: [ascotfeedback@gmail.com](ascotfeedback@gmail.com)

#### Usage instructions (we recommend a system with at least 30Gb ram):
```
git clone https://github.com/jpling/ascot.git
cd ./ascot/software/snaptron
make
```

#### To generate the splicing PSI data tables:
```
>Mouse cell types and tissues from bulk RNA-Seq (MESA) - ## datasets, ## exon query
    python3 ascot_psi.py --i ./exons/mesa_exons.tsv --a mesaall --c mesalinked --o mesa_psi.tsv

>Mouse single-cell RNA-Seq data (CellTower) - ## datasets, ## exon query
    python3 ascot_psi.py --i ./exons/ctms_exons.tsv --a ctmsall --c ctmslinked --min 10 --f 1 --o ctms_psi.tsv

>Human GTEx tissues + eye (GTEx) - ## datasets, ## exon query
    python3 ascot_psi.py --i ./exons/gtexeye_exons_1.tsv --a gtexeyeall --c gtexeyelinked --o gtexeye_psi_part1.tsv
    python3 ascot_psi.py --i ./exons/gtexeye_exons_2.tsv --a gtexeyeall --c gtexeyelinked --o gtexeye_psi_part2.tsv
    python3 ascot_psi.py --i ./exons/gtexeye_exons_3.tsv --a gtexeyeall --c gtexeyelinked --o gtexeye_psi_part3.tsv

>ENCODE shRNA-Seq knockdown data (ENCODE) - # datasets, ## exon query
    python3 ascot_psi.py --i ./exons/encodegtex_exons_1.tsv --a encodegtexall --c encodegtexlinked --o encodegtex_psi_part1.tsv
    python3 ascot_psi.py --i ./exons/encodegtex_exons_2.tsv --a encodegtexall --c encodegtexlinked --o encodegtex_psi_part2.tsv
```

#### To generate the gene expression NAUC data tables:
```
>Mouse cell types and tissues from bulk RNA-Seq (MESA)
    python3 ascot_nauc.py --a mesaall --c mesalinked --o mesa_nauc.tsv
    
>Mouse single-cell RNA-Seq data (CellTower)
    python3 ascot_nauc.py --a ctmsall --c ctmslinked --o ctms_nauc.tsv
    
>Human GTEx tissues + eye (GTEx)
    python3 ascot_nauc.py --a gtexeyeall --c gtexeyelinked --o gtexeye_nauc.tsv
    
>ENCODE shRNA-Seq knockdown data (ENCODE)
    python3 ascot_nauc.py --a encodegtexall --c encodegtexlinked --o encodegtex_nauc.tsv
```
