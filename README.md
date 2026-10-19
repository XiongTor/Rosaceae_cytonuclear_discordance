# Phylogenomic Dissection of Cytonuclear Discordance and Reticulate Evolution in Rosaceae

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Pipeline: Phylogenomics](https://img.shields.io/badge/Pipeline-Phylogenomics-blue.svg)]()
[![R: >=4.2.0](https://img.shields.io/badge/R-%3E%3D4.2.0-green.svg)]()
[![Python: >=3.9](https://img.shields.io/badge/Python-%3E%3D3.9-blue.svg)]()
[![Julia: >=1.9](https://img.shields.io/badge/Julia-%3E%3D1.9-purple.svg)]()

This repository contains the complete analytical pipelines, custom scripts, and figure-generation workflows for the study:  
**"Decomposing Phylogenetic  Discordance in Rosaceae: The Roles of Incomplete Lineage Sorting, Introgression, and Gene-Tree Error"**  

<p align="center">
  <img src="workflow/Rosaceae_Cytonuclear_Workflow.jpg" alt="Rosaceae Cytonuclear Analytical Workflow" width="100%">
</p>

---

## 1. Scientific Overview and Abstract

Phylogenetic relationships among the major lineages of Rosaceae have long remained contentious, characterized by pronounced conflict between plastid and nuclear phylogenies, as well as extensive topological discordance among nuclear gene trees. Here, we comprehensively investigate the nature, magnitude, and underlying drivers of phylogenetic discordance across Rosaceae by sampling 79 genera—representing approximately 75% of the naturally occurring generic diversity across all three subfamilies (Rosoideae, Amygdaloideae, Dryadoideae) and all tribes—along with three outgroup taxa.

We integrated three complementary genomic datasets:
1. **Chloroplast Dataset**
2. **Angiosperms353 Dataset**
3. **OrthoFinder Dataset**

### Key Findings
1. **Profound Cytonuclear Incongruence**: The chloroplast concatenated tree places Amygdaloideae as early-diverging/basal with maximum bootstrap support, whereas both concatenated and coalescent-based species trees derived from the two independent nuclear datasets unambiguously place Rosoideae as basal.
2. **Discordance Epicenter Mapping**: Genome-wide DensiTree clouds, multidimensional scaling (MDS) of Robinson-Foulds and path distances, and internode conflict metrics (Quartet Sampling QC/QD/QI, Internode Certainty IC/ICA, and the four-dimensional concordance vector $\psi_1-\psi_4$) identify Amygdaloideae—and particularly the tribe **Maleae**—as the primary epicenter of topological incongruence.
3. **Sequence Property Divergence**: Systematically quantifying 15 sequence and tree properties via AMAS and PhyKit, followed by principal component analysis (PCA) and property-matched subsampling experiments, demonstrates that intrinsic evolutionary heterogeneity (e.g., GC content, evolutionary rate, and proportion of variable sites) between nuclear and organellar genomes contributes significantly to cytonuclear divergence.
4. **Dissecting Biological Drivers (ILS vs. IH)**: Multi-angle investigations—neutral multispecies coalescent (MSC) simulations, MSCquartets hypothesis testing under T1 and T3 models, QuIBL internal branch-length mixture modeling, genome-wide Dsuite ABBA-BABA, and SNaQ semi-directed phylogenetic networks—reveal that **Incomplete Lineage Sorting (ILS)** is the pervasive family-wide biological driver of gene-tree discordance, whereas **introgressive hybridization (IH)** is largely localized and concentrated within Maleae.
5. **Unified Hierarchical Importance Partitioning**: Regression of gene concordance factors (gCF) on estimated ILS ($\theta$), gene-tree estimation error (GTEE via Seq-Gen simulations), and introgression intensity (Reticulation Index, RI) shows that ILS and GTEE represent the dominant factors driving internode gene-tree conflict.
6. **Functional Adaptations of Introgressed Loci**: Whole-genome sliding-window $|f_{dM}|$ scans, GMAP chromosomal mapping onto the *Malus domestica* reference, and clusterProfiler GO enrichment across high-confidence Maleae trios reveal that introgressed genes are significantly enriched in pathways related to environmental stress tolerance, metabolic homeostasis, and reproductive defense.

---

## 2. Repository Structure

The repository is structured into 11 dedicated, modular scientific analysis suites:

```text
Rosaceae_cytonuclear_discordance/
├── concordance_node_label/                               # Multimetric internode concordance factors
│   ├── concordance_factor_gcf_scf_qcf.sh                 # IQ-TREE 3 gene (gCF) and site (sCF) concordance factors
│   ├── gcf_scf.sh                                        # Summary and parsing of gCF and sCF output tables
│   ├── phytop.sh                                         # ASTRAL quartet concordance vector calculation (ψ1-ψ4)
│   ├── ICA_IC_value.sh                                   # RAxML Internode Certainty (IC and ICA) computation
│   ├── ICA_plot.R                                        # Plotting of Internode Certainty distributions
│   └── QS_test.sh                                        # Quartet Sampling calculation across species-tree nodes
│
├── cytonuclear_display/                                  # Visual conflict and topological distance projection
│   ├── densitree.R                                       # DensiTree visual cloud overlay of gene trees
│   ├── RF_path_MDS_new.R & RF_path_MDS.R                 # Robinson-Foulds and Path distance MDS 2D projection
│   ├── QS_visual.R & QS_test.sh                          # Quartet Sampling QC/QD/QI visualization
│   └── phyparts.sh & phyparts_pies.sh                    # PhyParts bipartition conflict pie-chart calculation
│
├── factor_alignment_tree/                                # Sequence & tree properties and matching experiments
│   ├── factor.sh & factor_parallel.sh                    # High-throughput property extraction via AMAS & PhyKit
│   ├── make_tree_factor.sh                               # Supermatrix building and ML tree inference for matched sets
│   ├── random_nuclear.sh                                 # Random nuclear gene sampling as background controls
│   └── factor_plot/                                      # Modular R pipeline for statistical analysis & plotting:
│       ├── 00-utils_functions.R                          # Shared utility functions and theme definitions
│       ├── 01_data_loading.R                             # Dataset ingestion and normalization
│       ├── 02_density_plots.R                            # Density distribution comparisons between nuclear & plastid
│       ├── 03_correlation_analysis.R                     # 15-property pairwise correlation matrices
│       ├── 04_pca_analysis.R                             # Principal Component Analysis (PCA)
│       ├── 05_rf_distance_analysis.R                     # Baseline topological distance evaluation
│       ├── 06_chloroplast_matching.R                     # Property-matched subsetting experiment
│       ├── 07_final_comparison.R                         # Topological divergence reassessment
│       ├── 08_density_based_rf_analysis.R                # Density-stratified RF distance analysis
│       ├── 09-run_all.R                                  # Master runner for property analysis
│       ├── 10_chloroplast_similarity_analysis.R          # Multidimensional distance quantification
│       ├── factor_plot.R                                 # Monolithic plotting script
│       ├── README.md                                     # Submodule documentation
│       └── README_RF_ANALYSIS.md                         # Documentation for topological distance comparisons
│
├── simulate_theta/                                       # Coalescent simulation under the MSC model
│   ├── get_theta.R                                       # Internal node mutation parameter θ calculation
│   ├── phybase_simulate.R                                # Phybase coalescent simulation across θ gradient (3,000 trees)
│   ├── count_four_aray.R & count_four_aray_lm.R          # Quartet frequency extraction & linear regression (R²)
│   ├── get_number.py & plot.R                            # Regression R² summary and scatter plotting
│   ├── rename.sh                                         # File management for simulation replicates
│   └── RAEDME.md                                         # Submodule documentation
│
├── MSCquartets/                                          # MSCquartets hypothesis testing
│   ├── MSCquartets.R                                     # Hypothesis testing under T1 (tree) and T3 (coalescent) models
│   ├── get_sp_ILS_IH_contribution.R                      # Species-pair level ILS vs IH quantification
│   ├── get_gene_ILS_IH_contribution.R                    # Locus-level ILS vs IH scoring
│   ├── analy_df_num.R                                    # Hypothesis test classification summarizer
│   └── results_classify.sh                               # Automation script for quartet partitioning
│
├── Get_ILS_IH_gene/                                      # Score visualization for ILS and IH loci
│   └── get_and_plot_IH_gene.R                            # Boxplots and barplots of locus-level ILS/IH scores
│
├── QuIBL/                                                # QuIBL branch-length mixture modeling
│   ├── run_quibl.sh                                      # Triplet branch-length distribution estimation
│   └── QuIBL_plot.R                                      # ΔBIC classification and introgression heatmaps
│
├── D-statistics/                                         # Genome-wide ABBA-BABA introgression tests
│   ├── msa2vcf_biallelic.py                              # Extraction of biallelic SNPs from MSA to VCF format
│   ├── D-statistics.sh                                   # Dsuite Dtrios, Fbranch, and visualization runner
│   ├── run_Dsuite.py & run_Dsuite_xt.py                  # Automated batch execution of Dsuite
│   ├── run_Dsuite_no_snp_thinning.py                     # Dsuite execution without SNP thinning
│   ├── run_Dsutie_xt_snp_thinning.py                     # Dsuite execution with physical distance SNP thinning
│   ├── make_runbash.sh & run.sh                          # Shell script generator and parallel execution runner
│   ├── plot_box.R                                        # Subfamily and inter-subfamily D-statistic boxplots
│   ├── sets.txt                                          # Species-to-clade assignment definitions
│   └── plot_order.txt                                    # Species tree topology ordering for Dsuite heatmaps
│
├── PhyloNetwork/                                         # SNaQ semi-directed phylogenetic network inference
│   ├── runSNaQ_hn.jl                                     # Julia PhyloNetworks SNaQ optimization (h = 0 to 6)
│   ├── PhyloNetworks.R                                   # Network topology and inheritance probability γ extraction
│   ├── plot_final_loglikelihood.R                        # Log-pseudolikelihood slope plot for optimal h selection
│   └── plot_result_julia.R                               # Network visualization via PhyloPlots
│
├── ILS_GTEE_IH/                                          # Integrated variance decomposition (GTEE vs ILS vs IH)
│   ├── ILS_GTEE_IH.SH                                    # Master pipeline runner
│   ├── get_node_inf.R                                    # Node ID and reticulation index extraction
│   ├── get_theta_value.R                                 # Per-node mutation parameter θ estimation
│   ├── README_ILS_IH_GTEE.md                             # Pipeline design and module documentation
│   ├── new_use_in_globel/                                # Global analysis components:
│   │   ├── relaimpo.R                                    # Hierarchical relative importance (LMG) with LOESS 95% CI
│   │   ├── extract_node_values_theta.py                  # Node feature table construction (with θ)
│   │   ├── extract_node_values_without_theta.py           # Node feature table construction (without θ)
│   │   ├── MSC_geneTr_simulator.R                        # Gene tree simulator under the coalescent
│   │   ├── simulator_tripletCounter_tripletMapper.sh     # Triplet enumeration pipeline
│   │   ├── triplet_mapper.py                             # Empirical vs simulated triplet mapper
│   │   ├── triple_frequency_counter.py                   # Triplet topology counter
│   │   ├── find_unbalanced_triplets.py                   # Identification of significantly unbalanced triplets
│   │   ├── clean_supp.R                                  # Tree support cleanup utility
│   │   └── make_treelength_gt_0.R                        # Branch length zero-elimination utility
│   └── new_use_in_GTEE_part/                             # GTEE quantification via Seq-Gen:
│       ├── generate_seqgen_commands.py                   # Seq-Gen simulation command generator
│       ├── seqgen_commands.sh                            # Batch execution of Seq-Gen simulations
│       ├── clip_missing_taxa.sh                          # Alignment pruning for empirical taxon occupancy
│       └── transform_astral_q1_result.py                 # Summary of simulated bipartition recovery
│
├── Introgression_Functional_Enrichment/                  # Introgressed gene detection & GO enrichment
│   ├── 01_filter_high_confidence_trios_f4_quibl.R        # Dual filtering (top decile f4-ratio & QuIBL mixprop2)
│   ├── 02_eggnog_orgdb_clusterprofiler_enrichment.R      # EggNOG parsing, custom OrgDb, and clusterProfiler
│   ├── 03_plot_GO_multi_threshold_barplots.R             # Multi-threshold (top 3% & top 5% |fdM|) GO barplots
│   ├── 04_analyze_GO_venn_six_groups_exact_intersections.R # Exact 6-trio Venn intersections
│   ├── 05_batch_assemble_trio_figures.py                 # Publication figure assembly (scatter, tree, GO)
│   └── README_Introgression_GO.md                        # Workflow documentation
│
├── workflow/                                             # Analytical pipeline overview and flowchart
│   └── Rosaceae_Cytonuclear_Workflow.jpg                 # Full analytical workflow diagram
│
└── code_curation_report.md                               # Comprehensive code audit and curation report
```

---

## 3. Computational Environment and Prerequisites

### 3.1 R Environment (R >= 4.2.0)
Install required CRAN, Bioconductor, and GitHub packages:
```R
# CRAN dependencies
install.packages(c(
  "ape", "phangorn", "phytools", "ggplot2", "dplyr", "tidyr", 
  "scales", "pheatmap", "corrplot", "relaimpo", "gridExtra", 
  "foreach", "doParallel", "reshape2", "psych", "TreeDist", 
  "Quartet", "dispRity", "MonoPhy"
))

# Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c(
  "ggtree", "treeio", "clusterProfiler", "AnnotationForge", 
  "enrichplot", "qvalue"
))

# GitHub dependencies
devtools::install_github("nstenz/quiblR")
```

### 3.2 Python Environment (Python >= 3.9)
```bash
pip install biopython pandas numpy scipy matplotlib seaborn reportlab PyPDF2 pdf2image
```

### 3.3 Julia Environment (for SNaQ Network Inference, Julia >= 1.9)
```julia
using Pkg
Pkg.add(["PhyloNetworks", "PhyloPlots", "CSV", "DataFrames"])
```

### 3.4 Main Bioinformatics Software
Ensure the following tools are compiled and available in system `$PATH`:
- **IQ-TREE 3** (v3.0.0 or later; Minh et al. 2020)
- **ASTRAL-Pro3 / ASTRAL-III** (Zhang et al. 2020)
- **RAxML / RAxML-NG** (Stamatakis 2014)
- **Dsuite** (v0.4 or later; Malinsky et al. 2021)
- **QuIBL** (Edelman et al. 2019)
- **Quartet Sampling** (Pease et al. 2018)
- **Seq-Gen** (v1.3.4; Rambaut and Grass 1997)
- **PhyKit & AMAS** (Steenwyk et al. 2021; Borowiec 2016)
- **EggNOG-mapper** (v2.1.13; Cantalapiedra et al. 2021)

---

## 4. Citation and Contact

If you find this codebase, analytical workflows, or custom scripts useful in your research, please cite our paper:

```bibtex
@article{xiong2026rosaceae,
  title={Phylogenomic Dissection of Cytonuclear Discordance and Reticulate Evolution in Rosaceae},
  author={Xiong, Tao and et al.},
  journal={In review},
  year={2026}
}
```

For questions regarding script execution, methodology, or data availability, please contact:
- **Tao Xiong** (xiongtao@... / GitHub: [@XiongTor](https://github.com/XiongTor))
- **Miao Sun** (miaosun@mail.hzau.edu.cn)