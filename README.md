# Filip1l DNA methylation in a murine cSCC model

Code for the manuscript *"Clinically relevant aberrant Filip1l DNA methylation detected in a murine model of cutaneous squamous cell carcinoma"* (Roth, Coussement et al., EBioMedicine, 2021, [doi:10.1016/j.ebiom.2021.103383](https://doi.org/10.1016/j.ebiom.2021.103383)).

We profiled cutaneous squamous cell carcinoma (cSCC) of SKH-1 hairless mice exposed to solar-simulated UV radiation with reduced representation bisulfite sequencing (RRBS) and oxidative RRBS (oxRRBS): 7 laser-capture microdissected tumours with matched control skin of the same mouse, plus 2 independent healthy skin samples of non-irradiated mice. One tumour and one control were sequenced in both pools as technical replicates. Differential methylation between tumour and control was tested per CpG (limma on M-values, adjusted for mouse) and per region (BiSeq). The mouse results were then compared with public human Infinium MethylationEPIC data of healthy skin, actinic keratosis (AK) and cSCC via human-mouse orthologues, and FILIP1L methylation was tested in high- vs low-risk human cSCC.

## Repository contents

| File | Purpose |
|---|---|
| `Import_RRBS.R` | Reads the Bismark `.cov` files (RRBS `B*` and oxRRBS `O*` samples) of both sequencing pools into one BiSeq `BSraw` object with the sample annotation; saves `allRRBS.rdata` |
| `Analysis_RRBS_mouse.R` | Main mouse analysis: replicate clustering, merging of the technical replicates, coverage filtering, global methylation and coverage tests, hydroxymethylation (RRBS − oxRRBS), per-CpG limma analysis, promoter/genic/intergenic annotation (biomaRt, Ensembl January 2019 archive), DMR analysis (BiSeq) |
| `Analysis_EPIC_human_Lyko.R` | Human EPIC analysis of the Rodríguez-Paredes et al. (Lyko lab, Nature Communications 2018) data and comparison with the mouse results: limma per contrast (healthy skin, AK, cSCC, stem-cell-like and keratinocyte-like subtypes), orthologue mapping, overlap of differentially methylated genes, genes of interest, keratin cluster, high- vs low-risk lesions (E-MTAB-8542) |
| `AdditionalFunctions.R` | Functions for average logFC and gene symbols per gene, orthologue mapping, overlap of the mouse and human gene sets, bedGraph tracks and DMR analysis (sourced by `Analysis_EPIC_human_Lyko.R`) |
| `Summary_RRBS_mouse.R` | Post-processing of `output_rrbs.txt`: summary bar plots, mapping statistics table, all CpGs in Filip1l |
| `ProcessBigWig.bash` | Public ChIP-seq tracks (CTCF, H3K27Ac, H3K27me3; GSM3259343-GSM3259351): bigWig → bedGraph, liftOver mm9 → mm10, merge of overlapping intervals, merge of replicates → one bigWig per mark |
| `mergeRegions.R` | Merges overlapping bedGraph intervals after liftOver (called by `ProcessBigWig.bash`) |
| `UCSCtracks_histonmod.R` | Downloads public histone modification peaks (H3K27Ac, H3K4me3, H3K27me3; GSM3634100-GSM3634105) and writes UCSC narrowPeak tracks in mm10 and mm9 (liftOver) |
| `conda_environment.yml` | Conda environment with the R packages and tools (see Installation for the versions) |
| `environment_example.yml` | Template for the local paths (copy to `environment.yml`) |

## Installation

```bash
git clone https://github.com/ljcousse/Filip1l_cSCC_mouse_methylation.git
cd Filip1l_cSCC_mouse_methylation

conda env create -f conda_environment.yml
conda activate cscc_mouse
```

The environment contains:
- Bioconductor: `BiSeq`, `GenomicRanges`, `GenomicFeatures`, `rtracklayer`, `limma`, `biomaRt`, `wateRmelon`, `methylumi`, `ChAMPdata`;
- CRAN: `yaml`, `gplots`, `ggplot2`, `xlsx`;
- UCSC tools `bigWigToBedGraph`, `liftOver`, `bedSort`, `bedGraphToBigWig`, `bigWigMerge`, `fetchChromSizes`, plus `pigz` and `wget` (for the genome browser tracks).

The pre-processing is not part of this repository: reads (NextSeq 500, 76 bp single-end) were trimmed with Trim Galore and mapped to GRCm38 with Bismark (v0.16.3, Bowtie2 mode, default settings). The analysis starts from the Bismark `.cov` files.

## Configuration

Copy the template and fill in your paths:

```bash
cp environment_example.yml environment.yml
```

`environment.yml` contains:
- `dir_pool20`: Bismark `.cov` files of sequencing pool 20;
- `dir_pool16`: Bismark `.cov` files of sequencing pool 16;
- `dir_work`: working directory, all input files are expected here and all output is written here.

Run the scripts from the repository directory (they read `environment.yml` from there and then change to `dir_work`). Keep `environment.yml` local (it is listed in `.gitignore`).

## Running the analysis

1. **Import:** `Rscript Import_RRBS.R` (Bismark `.cov` files → `BSraw` object `allRRBS.rdata` in `dir_work`). 36 libraries: 18 RRBS and 18 oxRRBS.
2. **Mouse analysis:** `Rscript Analysis_RRBS_mouse.R`. Technical replicates are clustered (top 10,000 most variable CpGs) and then merged by summing the read counts. Only RRBS CpGs covered by at least 8 reads in all samples and with at least 6 methylated reads in total are kept; the oxRRBS data are restricted to the same CpGs. The RRBS − oxRRBS difference (hydroxymethylation) was mainly noise, so the further analysis uses the RRBS data only. The DMR analysis uses 20 cores (`mc.cores = 20`); change this for smaller machines. Needs `Mouse_CGI.txt` and downloads the Ensembl gene annotation once (`genes.Rda`).
3. **Comparison with human data:** `Rscript Analysis_EPIC_human_Lyko.R` (EPIC `.idat` files → M-values → limma per contrast → overlap with `output_rrbs.txt` via orthologues). `AdditionalFunctions.R` uses 25 cores (`mc.cores = 25`). Plots go to `dir_work/Plots/` and `dir_work/Risk_analysis/`.
4. **Post-processing:** `Rscript Summary_RRBS_mouse.R` (`output_rrbs.txt` → summary plots, `table_mapstat.txt`, `Filip1l_RRBS_mouse.txt`).
5. **Genome browser tracks:** `bash ProcessBigWig.bash` (merged wild-type ChIP-seq bigWigs in `dir_work/chipseq_tracks/` → `CTCF.bw`, `H3K27Ac.bw`, `H3K27Me3.bw` in mm10; only files with `WT` and `merged.bw` in their name are processed). `Rscript UCSCtracks_histonmod.R` (GEO peak files → `*_mm10_narrowPeak.txt`, `*_mm9_narrowPeak.txt`). Both scripts download the UCSC liftOver chain files they need.

### Main outputs

| Output | Content |
|---|---|
| `output_rrbs.txt` | limma results per CpG (tumour vs control) with gene annotation, mean and SD per group, and CGI context (`cgi_sort`: island, shore (≤ 2 kb), shelf (2-4 kb), open sea) |
| `output_rrbs_suppl.txt` | `output_rrbs.txt` with all overlapping genes and gene symbols per CpG |
| `DMR_analysis.txt` | Differentially methylated regions (BiSeq) with annotation |
| `Betas_rrbs.txt`, `Annot_rrbs.txt` | Methylation degree per CpG and sample, and the CpG positions |
| `results_lyko_Ker.txt`, `results_lyko_Stem.txt` | Human EPIC results, healthy skin vs keratinocyte-like and stem-cell-like tumours |
| `CpG_in_GOI.txt` | Number of CpGs in the genes of interest (mouse and human) |
| `Methylation_plot_paper.tiff`, `Methylation_plot_paper_small.tiff` | Global methylation, tumour vs control (all CpGs, and CpGs after filtering) |
| `RRBS_tumour.BedGraph`, `RRBS_control.BedGraph` | UCSC tracks with the mean methylation per CpG in tumour and control; significant CpGs (FDR < 0.10) as positive values, the other CpGs as negative values |
| `CTCF.bw`, `H3K27Ac.bw`, `H3K27Me3.bw`, `*_narrowPeak.txt` | Public ChIP-seq tracks and histone peaks for the UCSC genome browser |
| `Clustering.jpg`, `Risk_analysis_FILIP1L_sign.jpg` | Replicate clustering, FILIP1L in high- vs low-risk lesions |

### Statistical design (short)

- **Global methylation (mouse):** mean methylation per sample, ANOVA `y ~ tissue + mouse` and `y ~ tissue`, and a paired t-test of the 7 tumour/control pairs.
- **Per CpG (mouse):** limma on M-values (`log2((β + 0.01)/(1 − β + 0.01))`), `~ tissue + mouse`, coefficient tissue, Benjamini-Hochberg correction. The independent healthy skin samples are controls with their own mouse level.
- **DMRs (mouse):** BiSeq, CpG clusters with at least 15 CpGs, maximum distance 200 bp and coverage in 75% of the samples; beta regression (probit link) `~ tissue + mouse`; cluster FDR 0.1, CpG FDR 0.05. The `FDR` column of `DMR_analysis.txt` is Benjamini-Hochberg on the median p-value per DMR (the published version of the code used the stricter Holm correction here).
- **Human EPIC:** limma on M-values per contrast (healthy skin vs AK and/or cSCC, AK vs cSCC, stem-cell-like vs keratinocyte-like, high- vs low-risk). Differentially methylated genes: adjusted p < 0.05 (human) and < 0.10 (mouse). Overlap of the gene sets tested with chi-squared tests, direction of change with a binomial test.

Promoters are defined as −2,000 to +500 bp around the transcription start site; genic CpGs are the remaining CpGs within a gene.

## Data

The mouse RRBS and oxRRBS data are available in NCBI GEO under accession [GSE151657](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151657).

The public data used:
- human EPIC data of healthy skin, AK and cSCC: ArrayExpress [E-MTAB-5738](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-5738);
- human EPIC data of high- and low-risk lesions: ArrayExpress [E-MTAB-8542](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-8542);
- mouse ChIP-seq tracks: GEO GSM3259343-GSM3259351 and GSM3634100-GSM3634105.

The analysis expects these input files in `dir_work`:
- `allRRBS.rdata` (output of `Import_RRBS.R`);
- `Mouse_CGI.txt` (CpG island positions, UCSC Genome Browser CpG islands track for mm10; columns `chrom`, `chromStart`, `chromEnd`);
- `idats/`, `idats_bis/` (E-MTAB-5738 `.idat` files, split in two folders because the files have different lengths) and `E-MTAB-5738.sdrf.txt`;
- `idats_Sandoval/`, `E-MTAB-8542.sdrf.txt` and `Sample_info_Sandoval.txt` (lesion type per sample);
- `CpG_panel.csv` (DNA methylation prognosis signature of 32 CpGs predicting bad prognosis in cSCC, from the supplementary data of Hervás-Marín et al. (E-MTAB-8542). Expected columns: `IlmnID`, coefficient, direction of change (`UP`/`DOWN`), `CHR`, `MAPINFO`, then the Illumina annotation including `UCSC_RefGene_Name`);
- `Pooling info_OXRRBSAlbena.xlsx` and `amountCpGs.txt` (mapping statistics, for `Summary_RRBS_mouse.R`);
- `chipseq_tracks/` with the merged bigWig files of GSM3259343-GSM3259351.

## License

This code is released under the MIT License (see `LICENSE`).

## Citation

If you use this code, please cite:

Roth K\*, Coussement L\*, Knatko EV, Higgins M, Steyaert S, Proby CM, De Meyer T, Dinkova-Kostova AT. Clinically relevant aberrant Filip1l DNA methylation detected in a murine model of cutaneous squamous cell carcinoma. EBioMedicine. 2021;67:103383. https://doi.org/10.1016/j.ebiom.2021.103383.

\* shared first authorship
