# Phylogenomic Dissection of Cytonuclear Discordance and Reticulate Evolution in Rosaceae

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Pipeline: Phylogenomics](https://img.shields.io/badge/Pipeline-Phylogenomics-blue.svg)]()
[![R: >=4.2.0](https://img.shields.io/badge/R-%3E%3D4.2.0-green.svg)]()
[![Python: >=3.9](https://img.shields.io/badge/Python-%3E%3D3.9-blue.svg)]()
[![Julia: >=1.9](https://img.shields.io/badge/Julia-%3E%3D1.9-purple.svg)]()

This repository contains the complete analytical pipelines, custom scripts, and figure-generation workflows for the study:  
**"Phylogenomic Dissection of Cytonuclear Discordance and Reticulate Evolution in Rosaceae"**  
*(《蔷薇科核质冲突及其成因探究》)*

---

## 1. Scientific Overview and Abstract

Phylogenetic relationships among the major lineages of Rosaceae have long remained contentious, characterized by pronounced conflict between plastid and nuclear phylogenies, as well as extensive topological discordance among nuclear gene trees. Here, we comprehensively investigate the nature, magnitude, and underlying drivers of phylogenetic discordance across Rosaceae by sampling 79 genera—representing approximately 75% of the naturally occurring generic diversity across all three subfamilies (Rosoideae, Amygdaloideae, Dryadoideae) and all tribes—along with three outgroup taxa.

We integrated three complementary genomic datasets:
1. **Chloroplast Dataset**: Complete and partial plastid genomes assembled from 35 reference plastomes.
2. **Angiosperms353 Dataset**: Target-enrichment single-copy nuclear orthologs (353 loci).
3. **OrthoFinder Dataset**: Broad-scale transcriptome and whole-genome orthologous groups (2,987 single-copy orthologous genes derived from six reference genomes across diverse tribes).

### Key Findings
1. **Profound Cytonuclear Incongruence**: The chloroplast concatenated tree places Amygdaloideae as early-diverging/basal with maximum bootstrap support, whereas both concatenated and coalescent-based species trees derived from the two independent nuclear datasets unambiguously place Rosoideae as basal.
2. **Discordance Epicenter Mapping**: Genome-wide DensiTree clouds, multidimensional scaling (MDS) of Robinson-Foulds and path distances, and internode conflict metrics (Quartet Sampling QC/QD/QI, Internode Certainty IC/ICA, and the four-dimensional concordance vector $\psi_1-\psi_4$) identify Amygdaloideae—and particularly the tribe **Maleae**—as the primary epicenter of topological incongruence.
3. **Sequence Property Divergence**: Systematically quantifying 15 sequence and tree properties via AMAS and PhyKit, followed by principal component analysis (PCA) and property-matched subsampling experiments, demonstrates that intrinsic evolutionary heterogeneity (e.g., GC content, evolutionary rate, and proportion of variable sites) between nuclear and organellar genomes contributes significantly to cytonuclear divergence.
4. **Dissecting Biological Drivers (ILS vs. IH)**: Multi-angle investigations—neutral multispecies coalescent (MSC) simulations, MSCquartets hypothesis testing under T1 and T3 models, QuIBL internal branch-length mixture modeling, genome-wide Dsuite ABBA-BABA $D$-statistics/$f_4$-ratio/$F_{\mathrm{branch}}$, and SNaQ semi-directed phylogenetic networks—reveal that **Incomplete Lineage Sorting (ILS)** is the pervasive family-wide biological driver of gene-tree discordance, whereas **introgressive hybridization (IH)** is largely localized and concentrated within Maleae.
5. **Unified Hierarchical Importance Partitioning**: Regression of gene concordance factors (gCF) on estimated ILS ($\theta$), gene-tree estimation error (GTEE via Seq-Gen simulations), and introgression intensity (Reticulation Index, RI) shows that ILS and GTEE represent the dominant factors driving internode gene-tree conflict.
6. **Functional Adaptations of Introgressed Loci**: Whole-genome sliding-window $|f_{dM}|$ scans, GMAP chromosomal mapping onto the *Malus domestica* reference, and clusterProfiler GO enrichment across high-confidence Maleae trios reveal that introgressed genes are significantly enriched in pathways related to environmental stress tolerance, metabolic homeostasis, and reproductive defense.

---

## 2. Repository Structure

```text
Rosaceae_cytonuclear_discordance/
├── Analysis/
│   ├── hybpiper_analy/                                   # Target capture assembly & quality assessment
│   │   ├── trimmomatic.sh                                # Raw WGS/transcriptome read trimming & QC
│   │   ├── hybpiper_step.sh                              # HybPiper assembly of target nuclear orthologs
│   │   ├── hybpiper_paftol_combine.sh                    # Extraction and merging of target sequences
│   │   ├── make_tree_step.sh                             # MAFFT alignment, trimAl trimming, and preliminary gene trees
│   │   └── plot_angiosperms353_recovery_boxplot.R        # Boxplot of gene recovery across sequencing strategies
│   │
│   ├── bulidtree/                                        # Phylogenomic species tree & gene tree inference
│   │   ├── make_tree.sh                                  # IQ-TREE ML supermatrix & ASTRAL-Pro3 coalescent trees
│   │   ├── make_condence_tree.sh                         # Condensation of species trees to tribe and subfamily levels
│   │   ├── Sortadate.sh                                  # SortaDate clock-like locus filtering (--order 3,1,2)
│   │   ├── treeshrink_pipeline.sh                        # Outgroup-aware TreeShrink long-branch filtering (-q 0.2 -b 20)
│   │   ├── genetree_reroot.sh                            # Gene-tree outgroup rerooting via nw_labels and pxrr
│   │   ├── tree_reroot_mono.sh                           # Outgroup rerooting and tip name standardization
│   │   ├── hybsuite_use.sh                               # HybSuite paralog pruning workflow execution
│   │   └── plot_tree.R                                   # Publication-ready phylogenetic tree plotting (ggtree)
│   │
│   ├── ILS/
│   │   ├── cytonuclear_display/                          # Visual conflict and topological distance projection
│   │   │   ├── densitree.R                               # DensiTree visual cloud overlay of gene trees
│   │   │   ├── RF_path_MDS_new.R & RF_path_MDS.R         # Robinson-Foulds and Path distance MDS 2D projection
│   │   │   ├── QS_test.sh & QS_visual.R                  # Quartet Sampling QC/QD/QI visualization
│   │   │   └── phyparts.sh & phyparts_pies.sh            # PhyParts bipartition conflict pie-chart calculation
│   │   │
│   │   ├── concordance_node_label/                       # Multimetric internode concordance factors
│   │   │   ├── QS_test.sh                                # Quartet Sampling calculation across species-tree nodes
│   │   │   ├── ICA_IC_value.sh & ICA_plot.R              # RAxML Internode Certainty (IC and ICA) computation
│   │   │   ├── concordance_factor_gcf_scf_qcf.sh         # IQ-TREE 3 gene (gCF) and site (sCF) concordance factors
│   │   │   ├── gcf_scf.sh                                # Parsing gCF and sCF output tables
│   │   │   └── phytop.sh                                 # ASTRAL quartet concordance vector calculation (ψ1-ψ4)
│   │   │
│   │   ├── factor_alignment_tree/                        # Sequence & tree properties and matching experiments
│   │   │   ├── factor.sh & factor_parallel.sh            # High-throughput property extraction via AMAS & PhyKit
│   │   │   ├── make_tree_factor.sh                       # Supermatrix building and ML tree inference for matched sets
│   │   │   ├── random_nuclear.sh                         # Random nuclear gene sampling as background controls
│   │   │   └── factor_plot/                              # Modular R pipeline for statistical analysis & plotting:
│   │   │       ├── 00-utils_functions.R                  # Shared utility functions and theme definitions
│   │   │       ├── 01_data_loading.R                     # Dataset ingestion and normalization
│   │   │       ├── 02_density_plots.R                    # Density distribution comparisons between nuclear & plastid
│   │   │       ├── 03_correlation_analysis.R             # 15-property pairwise correlation matrices
│   │   │       ├── 04_pca_analysis.R                     # Principal Component Analysis (PCA)
│   │   │       ├── 05_rf_distance_analysis.R             # Baseline topological distance evaluation
│   │   │       ├── 06_chloroplast_matching.R             # Property-matched subsetting experiment
│   │   │       ├── 07_final_comparison.R                 # Topological divergence reassessment
│   │   │       ├── 08_density_based_rf_analysis.R        # Density-stratified RF distance analysis
│   │   │       ├── 09-run_all.R                          # Master runner for property analysis
│   │   │       ├── 10_chloroplast_similarity_analysis.R  # Multidimensional distance quantification
│   │   │       ├── factor_plot.R                         # Monolithic plotting script
│   │   │       ├── README.md                             # Submodule documentation
│   │   │       └── README_RF_ANALYSIS.md                 # Documentation for topological distance comparisons
│   │   │
│   │   ├── simulate_theta/                               # Coalescent simulation under the MSC model
│   │   │   ├── get_theta.R                               # Internal node mutation parameter θ calculation
│   │   │   ├── phybase_simulate.R                        # Phybase coalescent simulation across θ gradient (3,000 trees)
│   │   │   ├── count_four_aray.R & _lm.R                 # Quartet frequency extraction & linear regression (R²)
│   │   │   ├── get_number.py & plot.R                    # Regression R² summary and scatter plotting
│   │   │   ├── rename.sh                                 # File management for simulation replicates
│   │   │   └── RAEDME.md                                 # Submodule documentation
│   │   │
│   │   ├── MSCquartets/                                  # MSCquartets hypothesis testing
│   │   │   ├── MSCquartets.R                             # Hypothesis testing under T1 (tree) and T3 (coalescent) models
│   │   │   ├── get_sp_ILS_IH_contribution.R              # Species-pair level ILS vs IH quantification
│   │   │   ├── get_gene_ILS_IH_contribution.R            # Locus-level ILS vs IH scoring
│   │   │   ├── analy_df_num.R                            # Hypothesis test classification summarizer
│   │   │   └── results_classify.sh                       # Automation script for quartet partitioning
│   │   │
│   │   ├── Get_ILS_IH_gene/                              # Score visualization for ILS and IH loci
│   │   │   └── get_and_plot_IH_gene.R                    # Boxplots and barplots of locus-level ILS/IH scores
│   │   │
│   │   ├── QuIBL/                                        # QuIBL branch-length mixture modeling
│   │   │   ├── run_quibl.sh                              # Triplet branch-length distribution estimation
│   │   │   └── QuIBL_plot.R                              # ΔBIC classification and introgression heatmaps
│   │   │
│   │   ├── D-statistics/                                 # Genome-wide ABBA-BABA introgression tests
│   │   │   ├── msa2vcf_biallelic.py                      # Extraction of biallelic SNPs from MSA to VCF format
│   │   │   ├── D-statistics.sh                           # Dsuite Dtrios, Fbranch, and visualization runner
│   │   │   ├── run_Dsuite.py & run_Dsuite_xt.py          # Automated batch execution of Dsuite
│   │   │   ├── run_Dsuite_no_snp_thinning.py             # Dsuite execution without SNP thinning
│   │   │   ├── run_Dsutie_xt_snp_thinning.py             # Dsuite execution with physical distance SNP thinning
│   │   │   ├── make_runbash.sh & run.sh                  # Shell script generator and parallel execution runner
│   │   │   ├── plot_box.R                                # Subfamily and inter-subfamily D-statistic boxplots
│   │   │   ├── sets.txt                                  # Species-to-clade assignment definitions
│   │   │   └── plot_order.txt                            # Species tree topology ordering for Dsuite heatmaps
│   │   │
│   │   ├── PhyloNetwork/                                 # SNaQ semi-directed phylogenetic network inference
│   │   │   ├── runSNaQ_hn.jl                             # Julia PhyloNetworks SNaQ optimization (h = 0 to 6)
│   │   │   ├── PhyloNetworks.R                           # Network topology and inheritance probability γ extraction
│   │   │   ├── plot_final_loglikelihood.R                # Log-pseudolikelihood slope plot for optimal h selection
│   │   │   └── plot_result_julia.R                       # Network visualization via PhyloPlots
│   │   │
│   │   ├── ILS_GTEE_IH/                                  # Integrated variance decomposition (GTEE vs ILS vs IH)
│   │   │   ├── ILS_GTEE_IH.SH                            # Master pipeline runner
│   │   │   ├── get_node_inf.R                            # Node ID and reticulation index extraction
│   │   │   ├── get_theta_value.R                         # Per-node mutation parameter θ estimation
│   │   │   ├── README_ILS_IH_GTEE.md                     # Pipeline design and module documentation
│   │   │   ├── new_use_in_globel/                        # Global analysis components:
│   │   │   │   ├── relaimpo.R                            # Hierarchical relative importance (LMG) with LOESS 95% CI
│   │   │   │   ├── extract_node_values_theta.py          # Node feature table construction (with θ)
│   │   │   │   ├── extract_node_values_without_theta.py   # Node feature table construction (without θ)
│   │   │   │   ├── MSC_geneTr_simulator.R                # Gene tree simulator under the coalescent
│   │   │   │   ├── simulator_tripletCounter_tripletMapper.sh # Triplet enumeration pipeline
│   │   │   │   ├── triplet_mapper.py                     # Empirical vs simulated triplet mapper
│   │   │   │   ├── triple_frequency_counter.py           # Triplet topology counter
│   │   │   │   ├── find_unbalanced_triplets.py           # Identification of significantly unbalanced triplets
│   │   │   │   ├── clean_supp.R                          # Tree support cleanup utility
│   │   │   │   └── make_treelength_gt_0.R                # Branch length zero-elimination utility
│   │   │   └── new_use_in_GTEE_part/                     # GTEE quantification via Seq-Gen:
│   │   │       ├── generate_seqgen_commands.py           # Seq-Gen simulation command generator
│   │   │       ├── seqgen_commands.sh                    # Batch execution of Seq-Gen simulations
│   │   │       ├── clip_missing_taxa.sh                  # Alignment pruning for empirical taxon occupancy
│   │   │       └── transform_astral_q1_result.py         # Summary of simulated bipartition recovery
│   │   │
│   │   ├── tree_addoutg/                                 # Outgroup standardization for gene trees
│   │   │   ├── add_outgroup_gene_specific_change.py      # Gene-specific outgroup grafting based on sister taxa
│   │   │   └── keep_one_outg.sh                          # Standardization to a single representative outgroup
│   │   │
│   │   └── count_sp_rich_value.sh                        # Post-trimming gene occupancy summarizer per species
│   │
│   └── Introgression_Functional_Enrichment/              # Introgressed gene detection & GO enrichment
│       ├── 01_filter_high_confidence_trios_f4_quibl.R    # Dual filtering (top decile f4-ratio & QuIBL mixprop2)
│       ├── 02_eggnog_orgdb_clusterprofiler_enrichment.R  # EggNOG parsing, custom OrgDb, and clusterProfiler
│       ├── 03_plot_GO_multi_threshold_barplots.R         # Multi-threshold (top 3% & top 5% |fdM|) GO barplots
│       ├── 04_analyze_GO_venn_six_groups_exact_intersections.R # Exact 6-trio Venn intersections
│       ├── 05_batch_assemble_trio_figures.py             # Publication figure assembly (scatter, tree, GO)
│       └── README_Introgression_GO.md                    # Workflow documentation
│
├── Tools/                                                # General-purpose phylogenomic utilities
│   ├── cophylo.R                                         # Face-to-face cophylogenetic tanglegrams (ape/phytools)
│   ├── count_RF_path_two_trees.R                         # Pairwise RF and Path distance between two trees
│   ├── count_RF_path_dist.R                              # Distribution of RF and Path distances (sp vs genes)
│   ├── count_OG_per_species.sh                           # Summarizer of orthogroups per species from OrthoFinder
│   ├── count_length_gap.sh                               # Sequence length and gap proportion calculator
│   ├── Get_support_value.R                               # Extraction of ASTRAL support values with node IDs
│   ├── make_hybpiper_reference.sh                        # Reference fasta formatter for HybPiper
│   ├── mono.R                                            # Clade monophyly inspection using customized ape::zoom
│   ├── monophy.R                                         # Monophyly testing and outlier detection via MonoPhy
│   └── remove_lt_100bp_length_moregap.sh                 # Alignment length (<100 bp) and excessive gap filter
│
├── picture/                                              # Vector and raster image processing tools
│   ├── pdf_to_png.py                                     # High-DPI PDF to PNG batch converter
│   ├── concat_images.py                                  # Grid and multi-panel image concatenation
│   ├── combine_pdf.py                                    # Multi-page PDF merger
│   ├── combine_png.py                                    # Image grid stitcher
│   ├── svg_to_pdf.py                                     # Vector SVG to PDF converter
│   └── svg_to_png.py                                     # Vector SVG to PNG rasterizer
│
└── imag/                                                 # Methodological schematics and flowcharts
    └── readme_flowwork/
        └── 蔷薇科网状进化及不完全谱系筛选研究.svg         # Complete methodological flowchart
```

---

## 3. Computational Environment and Dependencies

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
*Note: `pdf2image` requires `poppler` binaries accessible in system PATH.*

### 3.3 Julia Environment (for SNaQ Network Inference, Julia >= 1.9)
```julia
using Pkg
Pkg.add(["PhyloNetworks", "PhyloPlots", "CSV", "DataFrames"])
```

### 3.4 Standalone Bioinformatics Software
Ensure the following tools are compiled and available in system `$PATH`:
- **IQ-TREE 3** (v3.0.0 or later; Minh et al. 2020)
- **ASTRAL-Pro3 / ASTRAL-III** (Zhang et al. 2020)
- **RAxML / RAxML-NG** (Stamatakis 2014)
- **Dsuite** (v0.4 or later; Malinsky et al. 2021)
- **QuIBL** (Edelman et al. 2019)
- **HybPiper** (v2.0 or later; Johnson et al. 2016)
- **HybSuite** (v1.1.0; Liu et al. 2026)
- **TreeShrink** (v1.3.5; Mai and Mirarab 2018)
- **SortaDate** (Smith et al. 2018)
- **Quartet Sampling** (Pease et al. 2018)
- **Seq-Gen** (v1.3.4; Rambaut and Grass 1997)
- **phyx** (Brown et al. 2017)
- **Newick Utilities** (Junier and Zdobnov 2010)
- **EggNOG-mapper** (v2.1.13; Cantalapiedra et al. 2021)
- **BWA, SAMtools, GATK v4.6+, VCFtools, PLINK v1.9, GMAP** (for upstream trio variant calling)

---

## 4. Step-by-Step Analytical Workflow

### Step 1: Quality Control, Assembly, and Sequence Curation
```bash
# 1. Adapter trimming and quality filtering of raw reads
bash Analysis/hybpiper_analy/trimmomatic.sh

# 2. HybPiper assembly of target nuclear loci
bash Analysis/hybpiper_analy/hybpiper_step.sh

# 3. Multiple sequence alignment (MAFFT) and automatic trimming (trimAl -auto)
bash Analysis/hybpiper_analy/make_tree_step.sh

# 4. Filter sequences shorter than 100 bp
bash Tools/remove_lt_100bp_length_moregap.sh

# 5. Outgroup-aware TreeShrink long-branch detection (-q 0.2 -b 20)
bash Analysis/bulidtree/treeshrink_pipeline.sh
```

### Step 2: Species Tree and Gene Tree Phylogenomics
```bash
# 1. Reroot gene trees using outgroup taxa
bash Analysis/bulidtree/genetree_reroot.sh

# 2. Maximum likelihood supermatrix concatenation (IQ-TREE 3) and coalescent (ASTRAL-Pro3) trees
bash Analysis/bulidtree/make_tree.sh

# 3. Visualize species trees with ggtree
Rscript Analysis/bulidtree/plot_tree.R
```

### Step 3: Cytonuclear Discordance and Internode Conflict Quantification
```bash
# 1. Face-to-face cophylogenetic tanglegrams (chloroplast vs. nuclear)
Rscript Tools/cophylo.R

# 2. Pairwise and distribution-level RF and Path topological distances
Rscript Tools/count_RF_path_two_trees.R <tree1.tre> <tree2.tre>
Rscript Tools/count_RF_path_dist.R <species_tree.tre> <gene_trees.tre>

# 3. DensiTree visual cloud overlay
Rscript Analysis/ILS/cytonuclear_display/densitree.R

# 4. Multidimensional scaling (MDS) projection of topological space
Rscript Analysis/ILS/cytonuclear_display/RF_path_MDS_new.R

# 5. Calculate gene and site concordance factors (gCF and sCF)
bash Analysis/ILS/concordance_node_label/concordance_factor_gcf_scf_qcf.sh

# 6. Internode Certainty (IC / ICA) via RAxML
bash Analysis/ILS/concordance_node_label/ICA_IC_value.sh
Rscript Analysis/ILS/concordance_node_label/ICA_plot.R

# 7. Quartet Sampling (QC, QD, QI) across species-tree nodes
bash Analysis/ILS/concordance_node_label/QS_test.sh
Rscript Analysis/ILS/cytonuclear_display/QS_visual.R

# 8. Four-dimensional concordance vector calculation (ψ1-ψ4)
bash Analysis/ILS/concordance_node_label/phytop.sh
```

### Step 4: Sequence and Tree Property-Matching Analysis
```bash
# 1. High-throughput extraction of 15 sequence and tree properties
bash Analysis/ILS/factor_alignment_tree/factor_parallel.sh <alignment_dir> <genetree_dir>

# 2. Execute correlation analysis, PCA, and property-matched sampling
Rscript Analysis/ILS/factor_alignment_tree/factor_plot/09-run_all.R
```

### Step 5: Distinguishing Incomplete Lineage Sorting and Introgression
```bash
# 1. Estimate mutation parameter θ and execute neutral coalescent simulations in Phybase
Rscript Analysis/ILS/simulate_theta/get_theta.R <coalescent_tree.tre> <mutation_tree.tre>
Rscript Analysis/ILS/simulate_theta/phybase_simulate.R
Rscript Analysis/ILS/simulate_theta/count_four_aray_lm.R

# 2. MSCquartets hypothesis testing under T1 (species tree) and T3 (coalescent) models
Rscript Analysis/ILS/MSCquartets/MSCquartets.R
Rscript Analysis/ILS/MSCquartets/get_sp_ILS_IH_contribution.R
Rscript Analysis/ILS/MSCquartets/get_gene_ILS_IH_contribution.R <tree_dir> <output_dir>

# 3. QuIBL branch-length mixture distribution modeling
bash Analysis/ILS/tree_addoutg/keep_one_outg.sh
bash Analysis/ILS/QuIBL/run_quibl.sh
Rscript Analysis/ILS/QuIBL/QuIBL_plot.R

# 4. Whole-genome ABBA-BABA D-statistics and f4-ratio via Dsuite
python Analysis/ILS/D-statistics/msa2vcf_biallelic.py <alignment.fasta> <output.vcf> --min_cov 0.8
bash Analysis/ILS/D-statistics/D-statistics.sh
Rscript Analysis/ILS/D-statistics/plot_box.R

# 5. Semi-directed phylogenetic network inference with SNaQ (h = 0 to 6)
julia Analysis/ILS/PhyloNetwork/runSNaQ_hn.jl
Rscript Analysis/ILS/PhyloNetwork/plot_final_loglikelihood.R
```

### Step 6: Unified Hierarchical Partitioning (GTEE vs. ILS vs. IH)
```bash
# 1. Simulate sequences with Seq-Gen along species-tree topology to quantify GTEE
python Analysis/ILS/ILS_GTEE_IH/new_use_in_GTEE_part/generate_seqgen_commands.py
bash Analysis/ILS/ILS_GTEE_IH/new_use_in_GTEE_part/seqgen_commands.sh

# 2. Enumerate rooted triplets in empirical and coalescent-simulated gene trees
bash Analysis/ILS/ILS_GTEE_IH/new_use_in_globel/simulator_tripletCounter_tripletMapper.sh

# 3. Execute relative importance (LMG) variance partitioning with LOESS 95% CI
Rscript Analysis/ILS/ILS_GTEE_IH/new_use_in_globel/relaimpo.R <input_metrics.csv>
```

### Step 7: Candidate Introgressed Gene Identification and Functional Enrichment
```bash
# 1. Dual-filtering of candidate trios based on top decile f4-ratio and QuIBL mixprop2
Rscript Analysis/Introgression_Functional_Enrichment/01_filter_high_confidence_trios_f4_quibl.R

# 2. Build custom OrgDb from EggNOG annotations and execute clusterProfiler GO enrichment
Rscript Analysis/Introgression_Functional_Enrichment/02_eggnog_orgdb_clusterprofiler_enrichment.R

# 3. Plot multi-threshold barplots (top 3% and top 5% |fdM|)
Rscript Analysis/Introgression_Functional_Enrichment/03_plot_GO_multi_threshold_barplots.R

# 4. Calculate exact 6-trio non-overlapping Venn intersections
Rscript Analysis/Introgression_Functional_Enrichment/04_analyze_GO_venn_six_groups_exact_intersections.R

# 5. Assemble publication figures combining scatter plots, trees, and GO functional profiles
python Analysis/Introgression_Functional_Enrichment/05_batch_assemble_trio_figures.py
```

---

## 5. Manuscript Figures and Tables to Code Mapping

| Manuscript Item | Description / Content | Underlying Script(s) in this Repository |
| :--- | :--- | :--- |
| **Figure 1** | Cytonuclear tanglegram, DensiTree clouds, and deep topological conflict | `Tools/cophylo.R`, `Analysis/ILS/cytonuclear_display/densitree.R` |
| **Figure 2** | Genus-level species tree with node concordance (QC, QD, QI, IC, ICA, $\psi_1-\psi_4$) | `Analysis/ILS/concordance_node_label/QS_visual.R`, `ICA_plot.R`, `phytop.sh` |
| **Figure 3** | 15 sequence/tree properties, correlation heatmap, PCA, and property-matching experiment | `Analysis/ILS/factor_alignment_tree/factor_plot/03_correlation_analysis.R`, `04_pca_analysis.R`, `06_chloroplast_matching.R` |
| **Figure 4** | Neutral coalescent simulation fit ($R^2$), QuIBL introgression heatmaps, and Dsuite $D$ distributions | `Analysis/ILS/simulate_theta/count_four_aray_lm.R`, `Analysis/ILS/QuIBL/QuIBL_plot.R`, `Analysis/ILS/D-statistics/plot_box.R` |
| **Figs. S1–S4** | Target capture evaluation, Angiosperms353 recovery, and paralog pruning strategies | `Analysis/hybpiper_analy/plot_angiosperms353_recovery_boxplot.R`, `Analysis/bulidtree/hybsuite_use.sh` |
| **Figs. S5–S6** | Branch support collapsing thresholds (<10%, <20%) and species-tree stability | `Analysis/bulidtree/make_tree.sh` |
| **Figs. S7–S11** | Robinson-Foulds and Path distance distributions across taxonomic levels | `Tools/count_RF_path_dist.R`, `Tools/count_RF_path_two_trees.R`, `Analysis/ILS/cytonuclear_display/RF_path_MDS_new.R` |
| **Figs. S12–S16** | Concatenation (IQ-TREE 3) and coalescent (ASTRAL-Pro3) nuclear and plastid phylogenies | `Analysis/bulidtree/make_tree.sh`, `Analysis/bulidtree/plot_tree.R` |
| **Figs. S17–S18** | Comprehensive internode concordance factors (QC, QD, QI, IC, ICA, gCF, sCF, $\psi$) | `Analysis/ILS/concordance_node_label/concordance_factor_gcf_scf_qcf.sh`, `ICA_IC_value.sh`, `QS_test.sh`, `phytop.sh` |
| **Figs. S21–S23** | Sequence property density curves, feature correlation matrices, and distance metrics | `Analysis/ILS/factor_alignment_tree/factor_plot/02_density_plots.R`, `07_final_comparison.R` |
| **Figs. S24–S30** | Coalescent simulation $R^2$ regression plots across $\theta$ gradients (family, subfamilies) | `Analysis/ILS/simulate_theta/count_four_aray_lm.R`, `Analysis/ILS/simulate_theta/plot.R` |
| **Figs. S31–S38** | MSCquartets hypothesis test ratios across $\alpha$ and $\beta$ thresholds, and locus-level scores | `Analysis/ILS/MSCquartets/MSCquartets.R`, `get_sp_ILS_IH_contribution.R`, `Analysis/ILS/Get_ILS_IH_gene/get_and_plot_IH_gene.R` |
| **Figs. S41–S45** | QuIBL introgression proportion heatmaps for AGS353 and OrthoFinder | `Analysis/ILS/QuIBL/QuIBL_plot.R` |
| **Figs. S46–S48** | Dsuite ABBA-BABA $D$-statistic and $f_4$-ratio distributions and spatial patterns | `Analysis/ILS/D-statistics/plot_box.R`, `Analysis/ILS/D-statistics/D-statistics.sh` |
| **Figs. S49–S52** | SNaQ semi-directed phylogenetic networks ($h = 0$ to $6$) and pseudolikelihood slopes | `Analysis/ILS/PhyloNetwork/runSNaQ_hn.jl`, `plot_final_loglikelihood.R` |
| **Figs. S53–S58** | Relative importance hierarchical variance decomposition (GTEE vs. ILS vs. IH) with LOESS | `Analysis/ILS/ILS_GTEE_IH/new_use_in_globel/relaimpo.R` |
| **Figs. S59–S64** | Multi-panel introgression scans ($|f_{dM}|$) and GO functional enrichment for 6 trios | `Analysis/Introgression_Functional_Enrichment/03_plot_GO_multi_threshold_barplots.R`, `05_batch_assemble_trio_figures.py` |
| **Figs. S69–S70** | Nuclear gene-tree DensiTree clouds and multidimensional scaling (MDS) | `Analysis/ILS/cytonuclear_display/densitree.R`, `RF_path_MDS_new.R` |
| **Tables S7–S8** | QuIBL discordant quartet distribution and mixture model fit statistics | `Analysis/ILS/QuIBL/QuIBL_plot.R` |
| **Table S22** | Approximately Unbiased (AU) topology test statistics | `Analysis/bulidtree/make_tree.sh` |
| **Table S25** | High-confidence introgressed gene candidates, functional annotations, and $P$-values | `Analysis/Introgression_Functional_Enrichment/02_eggnog_orgdb_clusterprofiler_enrichment.R` |

---

## 6. Citation and Contact

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