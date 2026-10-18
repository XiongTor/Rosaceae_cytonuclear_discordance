# Phylogenomic Dissection of Cytonuclear Discordance and Reticulate Evolution in Rosaceae

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Phylogenomics](https://img.shields.io/badge/Pipeline-Phylogenomics-blue.svg)]()
[![R-packages](https://img.shields.io/badge/R-%3E%3D4.2.0-green.svg)]()
[![Python](https://img.shields.io/badge/Python-%3E%3D3.9-blue.svg)]()

This repository contains the complete analytical pipelines, custom scripts, and figure-generation code for the study:  
**"Phylogenomic Dissection of Cytonuclear Discordance, Incomplete Lineage Sorting, and Reticulate Evolution in Rosaceae"**  
*(《蔷薇科核质冲突及其成因探究》)*

---

## 📖 Overview & Scientific Scope

Phylogenomic relationships among the three major subfamilies of Rosaceae (Rosoideae, Amygdaloideae, Dryadoideae) have long been controversial. By integrating **three independent genomic datasets**—complete/partial chloroplast genomes, target-enrichment nuclear genes (**Angiosperms353**), and transcriptome/whole-genome orthologs (**OrthoFinder**)—this study systematically dissects the macroevolutionary patterns and underlying drivers of phylogenetic discordance across Rosaceae:

1. **Deep Cytonuclear Discordance**: Nuclear data strongly support Rosoideae as basal, whereas plastid genomes place Amygdaloideae as basal.
2. **Topological & Node-level Quantification**: DensiTree clouds, multidimensional scaling (MDS), Robinson–Foulds distances, and node metrics (Quartet Sampling QC/QD/QI, Internode Certainty IC/ICA, and concordance vectors $\psi_1-\psi_4$) identify Amygdaloideae, particularly the tribe **Maleae**, as the conflict epicenter.
3. **Sequence Property Basis**: Systematic comparison of 15 sequence/tree properties and property-matching experiments demonstrate that evolutionary dynamics differences between organellar and nuclear genomes contribute substantially to cytonuclear divergence.
4. **Biological Drivers (ILS vs. IH)**: Multi-angle approaches—coalescent simulations under the multispecies coalescent (MSC), MSCquartets hypothesis testing, QuIBL internal branch-length modeling, Dsuite ABBA–BABA $D$-statistics/$f_4$-ratio, and SNaQ phylogenetic networks—demonstrate that **ILS** is the primary biological driver, accompanied by **episodic, localized introgression** predominantly within Maleae.
5. **Unified Relative Importance Framework**: Partitioning GTEE (gene-tree estimation error), ILS, and IH via hierarchical variance decomposition (`relaimpo`).
6. **Functional Heterogeneity of Introgressed Genes**: Sliding-window $|f_{dM}|$ scans, EggNOG annotation, and clusterProfiler functional enrichment across six high-confidence trios reveal polygenic adaptations to stress responses and metabolic homeostasis.

---

## 🗂 Repository Structure

```text
final_github_code/
├── Analysis/
│   ├── hybpiper_analy/                           # Module 1: HybPiper target assembly & quality control
│   │   ├── trimmomatic.sh                        # Read trimming and adapter removal
│   │   ├── hybpiper_step.sh                      # HybPiper assembly of Angiosperms353 loci
│   │   ├── hybpiper_paftol_combine.sh            # Target gene concatenation
│   │   ├── make_tree_step.sh                     # Gene tree generation pipeline
│   │   └── plot_angiosperms353_recovery_boxplot.R# Boxplot of locus recovery across sequencing strategies
│   │
│   ├── bulidtree/                                # Module 2: Species tree & gene tree construction
│   │   ├── make_tree.sh                          # IQ-TREE supermatrix & ASTRAL coalescent trees
│   │   ├── make_condence_tree.sh                 # Consensus tree generation
│   │   ├── Sortadate.sh                          # SortaDate clock-like locus filtering
│   │   ├── treeshrink_pipeline.sh                # TreeShrink long-branch outlier filtering
│   │   ├── genetree_reroot.sh                    # Gene tree outgroup rooting
│   │   ├── tree_reroot_mono.sh                   # Monophyly enforcement & rerooting
│   │   ├── hybsuite_use.sh                       # HybSuite pipeline runner
│   │   └── plot_tree.R                           # Phylogenetic tree visualization with ggtree
│   │
│   ├── ILS/
│   │   ├── cytonuclear_display/                  # Module 3: Visual conflict & topological distances
│   │   │   ├── densitree.R                       # DensiTree visual superposition of gene trees
│   │   │   ├── RF_path_MDS_new.R                 # Robinson-Foulds & Path distance MDS projection
│   │   │   ├── QS_test.sh & QS_visual.R          # Quartet Sampling analysis & visualization
│   │   │   └── phyparts.sh & phyparts_pies.sh    # PhyParts pie-chart bipartition conflicts
│   │   │
│   │   ├── concordance_node_label/               # Module 4: Multimetric node concordance
│   │   │   ├── QS_test.sh                        # Quartet Concordance (QC, QD, QI) computation
│   │   │   ├── ICA_IC_value.sh & ICA_plot.R      # RAxML Internode Certainty (IC, ICA)
│   │   │   ├── concordance_factor_gcf_scf_qcf.sh # IQ-TREE gene/site concordance factors (gCF/sCF)
│   │   │   └── phytop.sh                         # ASTRAL quartet concordance vectors (ψ1-ψ4)
│   │   │
│   │   ├── factor_alignment_tree/                # Module 5: 15 Sequence properties & matching
│   │   │   ├── factor.sh & factor_parallel.sh    # High-throughput property extraction
│   │   │   └── factor_plot/                      # 14 R scripts for PCA, correlation & matching:
│   │   │       ├── 01_data_loading.R             # Data ingestion and normalization
│   │   │       ├── 02_density_plots.R            # Property distribution comparison
│   │   │       ├── 03_correlation_analysis.R     # Property correlation matrix
│   │   │       ├── 04_pca_analysis.R             # PCA of nuclear vs plastid property space
│   │   │       ├── 05_rf_distance_analysis.R     # Baseline RF distance evaluation
│   │   │       ├── 06_chloroplast_matching.R     # Property-matched subsetting experiment
│   │   │       └── 07_final_comparison.R         # Multi-metric topological divergence re-assessment
│   │   │
│   │   ├── simulate_theta/                       # Module 6: Coalescent simulation under MSC
│   │   │   ├── get_theta.R                       # Branch-specific mutation parameter θ estimation
│   │   │   ├── phybase_simulate.R                # Phybase neutral MSC gene tree simulation
│   │   │   ├── count_four_aray.R & _lm.R         # MSCquartets topology frequency & linear regression
│   │   │   └── get_number.py & plot.R            # Goodness-of-fit (R²) extraction and plotting
│   │   │
│   │   ├── MSCquartets/                          # Module 7: MSCquartets hypothesis testing
│   │   │   ├── MSCquartets.R                     # Statistical hypothesis tests (T1 and T3 models)
│   │   │   ├── get_sp_ILS_IH_contribution.R      # Species-pair ILS vs IH quantification
│   │   │   └── get_gene_ILS_IH_contribution.R    # Locus-level ILS vs IH quantification
│   │   │
│   │   ├── QuIBL/                                # Module 8: QuIBL internal branch-length modeling
│   │   │   ├── run_quibl.sh                      # QuIBL execution across triplets
│   │   │   └── QuIBL_plot.R                      # ΔBIC filtering, introgression heatmaps
│   │   │
│   │   ├── D-statistics/                         # Module 9: Dsuite genome-wide introgression
│   │   │   ├── msa2vcf_biallelic.py              # MSA to biallelic SNP VCF conversion
│   │   │   ├── D-statistics.sh & run_Dsuite*.py  # D-statistic & f4-ratio with SNP thinning
│   │   │   └── plot_box.R                        # Subfamily-level introgression distribution boxplots
│   │   │
│   │   ├── PhyloNetwork/                         # Module 10: SNaQ phylogenetic networks
│   │   │   ├── runSNaQ_hn.jl                     # Julia PhyloNetworks SNaQ network inference (h=0..6)
│   │   │   ├── PhyloNetworks.R                   # Network parsing and inheritance probability γ
│   │   │   └── plot_final_loglikelihood.R        # Log-likelihood slope plot for best h selection
│   │   │
│   │   └── ILS_GTEE_IH/                          # Module 11: Unified relative importance partitioning
│   │       ├── new_use_in_globel/                # Full pipeline & relaimpo LMG decomposition
│   │       │   ├── relaimpo.R                    # Relative contribution of GTEE vs ILS vs IH
│   │       │   ├── extract_node_values_theta.py  # Branch length & θ feature extraction
│   │       │   └── triplet_mapper.py             # Empirical vs simulated triplet mapper
│   │       └── new_use_in_GTEE_part/             # Seq-Gen simulation for estimation error
│   │           └── generate_seqgen_commands.py   # Seq-Gen simulation command generator
│   │
│   └── Introgression_Functional_Enrichment/      # Module 12: Introgressed genes & GO/KEGG enrichment
│       ├── 01_filter_high_confidence_trios_f4_quibl.R  # Dual-filtering (top 10% f4-ratio + QuIBL)
│       ├── 02_eggnog_orgdb_clusterprofiler_enrichment.R # EggNOG annotation & OrgDb construction
│       ├── 03_plot_GO_multi_threshold_barplots.R       # GO enrichment barplots (top 3% & top 5% |fdM|)
│       ├── 04_analyze_GO_venn_six_groups_exact_intersections.R # Exact 6-trio Venn intersections
│       └── 05_batch_assemble_trio_figures.py           # Multi-panel vector figure assembly
│
├── Tools/                                        # Core bioinformatics utility scripts
│   ├── cophylo.R                                 # Face-to-face cophylogenetic tanglegrams
│   ├── count_RF_path_dist.R                      # Robinson-Foulds and path distance matrices
│   ├── count_OG_per_species.sh                   # OrthoFinder gene count summarizer
│   ├── msa2vcf_biallelic.py                      # Multi-alignment to VCF converter
│   ├── kaks.sh                                   # Ka/Ks calculator wrapper
│   ├── mono.R & monophy.R                        # Monophyly testing
│   └── remove_lt_100bp_length_moregap.sh         # Alignment length and gap filter
│
├── picture/                                      # Multi-panel figure assembly & format conversion
│   ├── concat_images.py & combine_pdf.py         # Vector & raster figure stitching
│   └── svg_to_pdf.py & pdf_to_png.py             # Format converters
│
├── imag/                                         # Project workflow diagrams
│   └── readme_flowwork/                          # Methodology flowcharts
│
├── code_curation_report.md                       # Comprehensive code audit and curation report
└── Gene_Introgression_Workflow.pdf               # Gene introgression workflow overview
```

---

## 💻 Software & Environmental Prerequisites

### 1. R Environment (R $\ge$ 4.2.0)
Install required packages via CRAN and Bioconductor:
```R
# CRAN packages
install.packages(c("ape", "phangorn", "phytools", "ggplot2", "dplyr", "tidyr", 
                   "scales", "pheatmap", "corrplot", "relaimpo", "pdfplumber"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("ggtree", "treeio", "clusterProfiler", "AnnotationForge", 
                       "enrichplot", "qvalue"))

# GitHub packages
devtools::install_github("nstenz/quiblR")
```

### 2. Python Environment (Python $\ge$ 3.9)
```bash
pip install biopython pandas numpy matplotlib seaborn scipy reportlab PyPDF2 pdf2image
```

### 3. Julia Environment (for SNaQ phylogenetic networks)
```julia
using Pkg
Pkg.add(["PhyloNetworks", "PhyloPlots", "CSV", "DataFrames"])
```

### 4. Standalone Bioinformatics Tools
Ensure the following command-line tools are installed and accessible in your `$PATH`:
- **IQ-TREE 3** (`iqtree2` or `iqtree3`, v3.0+)
- **ASTRAL-III / ASTRAL-MP / ASTRAL-pro**
- **RAxML / RAxML-NG** (v8.2.12+)
- **Dsuite** (v0.4+)
- **QuIBL** (Python-based branch length modeling)
- **HybPiper** (v2.0+)
- **TreeShrink** (v1.3+)
- **SortaDate** (Python)
- **Quartet Sampling** (`quartet_sampling.py`, Pease et al. 2018)
- **Seq-Gen** (v1.3.4+)
- **EggNOG-mapper** (v2.1.13+)

---

## 🚀 Step-by-Step Execution Guide

### Step 1: Target Enrichment & Phylogenomic Assembly
```bash
# 1. Clean reads
bash Analysis/hybpiper_analy/trimmomatic.sh

# 2. Extract Angiosperms353 loci with HybPiper
bash Analysis/hybpiper_analy/hybpiper_step.sh

# 3. Detect long-branch outliers
bash Analysis/bulidtree/treeshrink_pipeline.sh

# 4. Construct ML supermatrix & ASTRAL coalescent species trees
bash Analysis/bulidtree/make_tree.sh
```

### Step 2: Cytonuclear Discordance & Multimetric Node Quantification
```bash
# 1. Visualize gene tree discordance with DensiTree
Rscript Analysis/ILS/cytonuclear_display/densitree.R

# 2. Project gene tree topology space via MDS
Rscript Analysis/ILS/cytonuclear_display/RF_path_MDS_new.R

# 3. Compute Quartet Concordance (QC, QD, QI)
bash Analysis/ILS/concordance_node_label/QS_test.sh

# 4. Compute Internode Certainty (IC, ICA) with RAxML
bash Analysis/ILS/concordance_node_label/ICA_IC_value.sh

# 5. Compute Concordance Vectors (ψ1-ψ4)
bash Analysis/ILS/concordance_node_label/phytop.sh
```

### Step 3: Sequence & Tree Properties and Property-Matching Experiment
```bash
# 1. Extract 15 properties across nuclear and plastid alignments
bash Analysis/ILS/factor_alignment_tree/factor_parallel.sh

# 2. Run PCA, correlation, and property-matching experiment
Rscript Analysis/ILS/factor_alignment_tree/factor_plot/09-run_all.R
```

### Step 4: Dissecting Biological Drivers (ILS vs. IH)
```bash
# 1. Estimate θ and perform neutral MSC coalescent simulations
Rscript Analysis/ILS/simulate_theta/get_theta.R
Rscript Analysis/ILS/simulate_theta/phybase_simulate.R
Rscript Analysis/ILS/simulate_theta/count_four_aray_lm.R

# 2. MSCquartets T1 and T3 hypothesis tests
Rscript Analysis/ILS/MSCquartets/MSCquartets.R

# 3. QuIBL branch-length modeling
bash Analysis/ILS/QuIBL/run_quibl.sh
Rscript Analysis/ILS/QuIBL/QuIBL_plot.R

# 4. Dsuite genome-wide ABBA-BABA D-statistics
bash Analysis/ILS/D-statistics/D-statistics.sh
python Analysis/ILS/D-statistics/run_Dsuite_xt.py

# 5. SNaQ phylogenetic networks in Julia
julia Analysis/ILS/PhyloNetwork/runSNaQ_hn.jl
Rscript Analysis/ILS/PhyloNetwork/plot_final_loglikelihood.R
```

### Step 5: Unified Hierarchical Partitioning (GTEE vs. ILS vs. IH)
```bash
# 1. Simulate sequences with Seq-Gen to quantify estimation error
python Analysis/ILS/ILS_GTEE_IH/new_use_in_GTEE_part/generate_seqgen_commands.py
bash Analysis/ILS/ILS_GTEE_IH/new_use_in_GTEE_part/seqgen_commands.sh

# 2. Run relaimpo LMG relative contribution partitioning
Rscript Analysis/ILS/ILS_GTEE_IH/new_use_in_globel/relaimpo.R
```

### Step 6: Candidate Introgressed Gene Identification & Functional Enrichment
```bash
# 1. Filter candidate trios based on top decile f4-ratio and QuIBL
Rscript Analysis/Introgression_Functional_Enrichment/01_filter_high_confidence_trios_f4_quibl.R

# 2. Build custom OrgDb and run clusterProfiler GO/KEGG enrichment
Rscript Analysis/Introgression_Functional_Enrichment/02_eggnog_orgdb_clusterprofiler_enrichment.R

# 3. Plot multi-threshold barplots (top 3% & top 5% |fdM|)
Rscript Analysis/Introgression_Functional_Enrichment/03_plot_GO_multi_threshold_barplots.R

# 4. Calculate exact 6-trio Venn intersections
Rscript Analysis/Introgression_Functional_Enrichment/04_analyze_GO_venn_six_groups_exact_intersections.R
```

---

## 📊 Manuscript Figure & Table to Code Mapping Matrix

| Manuscript Item | Description | Generating / Analytical Scripts in this Repository |
| :--- | :--- | :--- |
| **Figure 1** | Cytonuclear tanglegram, DensiTree clouds, and basal topological conflicts | `Tools/cophylo.R`, `Analysis/ILS/cytonuclear_display/densitree.R` |
| **Figure 2** | Genus-level species tree, node concordance (QC, QD, QI, IC, ICA, $\psi_1-\psi_4$) | `Analysis/ILS/concordance_node_label/QS_visual.R`, `ICA_plot.R`, `phytop.sh` |
| **Figure 3** | 15 Sequence properties, correlation heatmap, PCA, and property-matching experiment | `Analysis/ILS/factor_alignment_tree/factor_plot/03_correlation_analysis.R`, `04_pca_analysis.R`, `06_chloroplast_matching.R` |
| **Figure 4** | Coalescent simulation fit ($R^2$), QuIBL heatmaps, and Dsuite $D$-statistic distributions | `Analysis/ILS/simulate_theta/count_four_aray_lm.R`, `Analysis/ILS/QuIBL/QuIBL_plot.R`, `Analysis/ILS/D-statistics/plot_box.R` |
| **Figs. S1–S6** | Data quality control, Angiosperms353 recovery, and paralog pruning | `Analysis/hybpiper_analy/plot_angiosperms353_recovery_boxplot.R`, `Tools/count_OG_per_species.sh` |
| **Figs. S7–S11** | Robinson-Foulds and Path distance distributions across taxonomic levels | `Tools/count_RF_path_dist.R`, `Analysis/ILS/cytonuclear_display/RF_path_MDS_new.R` |
| **Figs. S12–S16** | Concatenation (IQ-TREE) and coalescent (ASTRAL) nuclear/plastid species trees | `Analysis/bulidtree/make_tree.sh`, `plot_tree.R` |
| **Figs. S17–S18** | Comprehensive node concordance metrics (QC, QD, QI, IC, ICA, gCF, sCF, $\psi$) | `Analysis/ILS/concordance_node_label/QS_test.sh`, `ICA_IC_value.sh`, `phytop.sh` |
| **Figs. S21–S23** | Sequence property density curves, feature correlation matrices, and distance metrics | `Analysis/ILS/factor_alignment_tree/factor_plot/02_density_plots.R`, `07_final_comparison.R` |
| **Figs. S24–S30** | Coalescent simulation $R^2$ regression plots across $\theta$ gradients (family, subfamilies) | `Analysis/ILS/simulate_theta/count_four_aray_lm.R`, `plot.R` |
| **Figs. S31–S38** | MSCquartets hypothesis test ratios (T1 and T3 models) across $\alpha$ and $\beta$ thresholds | `Analysis/ILS/MSCquartets/MSCquartets.R`, `get_sp_ILS_IH_contribution.R` |
| **Figs. S41–S45** | QuIBL introgression proportion heatmaps for AGS353 and OrthoFinder | `Analysis/ILS/QuIBL/QuIBL_plot.R` |
| **Figs. S46–S48** | Dsuite ABBA-BABA $D$-statistic and $f_4$-ratio boxplots and spatial patterns | `Analysis/ILS/D-statistics/plot_box.R` |
| **Figs. S49–S52** | SNaQ semi-directed phylogenetic networks ($h=0$ to $h=6$) and log-likelihood curves | `Analysis/ILS/PhyloNetwork/runSNaQ_hn.jl`, `plot_final_loglikelihood.R` |
| **Figs. S53–S58** | Relative importance hierarchical variance decomposition (GTEE vs ILS vs IH) | `Analysis/ILS/ILS_GTEE_IH/new_use_in_globel/relaimpo.R` |
| **Figs. S59–S64** | Multi-panel introgression scans ($|f_{dM}|$) and GO functional enrichment for 6 trios | `Analysis/Introgression_Functional_Enrichment/03_plot_GO_multi_threshold_barplots.R`, `05_batch_assemble_trio_figures.py` |
| **Figs. S69–S70** | Nuclear gene tree DensiTree clouds and multidimensional scaling (MDS) | `Analysis/ILS/cytonuclear_display/densitree.R`, `RF_path_MDS_new.R` |
| **Tables S7–S8** | QuIBL discordant quartet distribution and model fit statistics | `Analysis/ILS/QuIBL/QuIBL_plot.R` |
| **Table S22** | Approximately Unbiased (AU) topology test statistics | `Analysis/bulidtree/make_tree.sh` |
| **Table S25** | High-confidence introgressed gene candidate functional annotations and $P$-values | `Analysis/Introgression_Functional_Enrichment/02_eggnog_orgdb_clusterprofiler_enrichment.R` |

---

## 📄 License & Attribution

This project is licensed under the MIT License - see the LICENSE file for details.

If you find this repository or associated methodologies useful in your research, please cite our manuscript:
```bibtex
@article{xiong2026rosaceae,
  title={Phylogenomic Dissection of Cytonuclear Discordance, Incomplete Lineage Sorting, and Reticulate Evolution in Rosaceae},
  author={Xiong, Tao and et al.},
  journal={In review},
  year={2026}
}
```

For questions, bug reports, or collaborative inquiries, please contact:  
**Tao Xiong** (xiongtao@... / GitHub: [@XiongTor](https://github.com/XiongTor))