# Phylogenomic Dissection of Cytonuclear Discordance and Reticulate Evolution in Rosaceae

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Pipeline: Phylogenomics](https://img.shields.io/badge/Pipeline-Phylogenomics-blue.svg)]()
[![R: >=4.2.0](https://img.shields.io/badge/R-%3E%3D4.2.0-green.svg)]()
[![Python: >=3.9](https://img.shields.io/badge/Python-%3E%3D3.9-blue.svg)]()
[![Julia: >=1.9](https://img.shields.io/badge/Julia-%3E%3D1.9-purple.svg)]()

This repository contains the complete analytical pipelines, custom scripts for the study:  
**"Decomposing Phylogenetic  Discordance in Rosaceae: The Roles of Incomplete Lineage Sorting, Introgression, and Gene-Tree Error"**  

<p align="center">
  <img src="workflow/Rosaceae_Cytonuclear_Workflow.jpg" alt="Rosaceae Cytonuclear Analytical Workflow" width="100%">
</p>

---



##  Computational Environment and Prerequisites

### R Environment (R >= 4.2.0)
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

### Python Environment (Python >= 3.9)
```bash
pip install biopython pandas numpy scipy matplotlib seaborn reportlab PyPDF2 pdf2image
```

### Julia Environment (for SNaQ Network Inference, Julia >= 1.9)
```julia
using Pkg
Pkg.add(["PhyloNetworks", "PhyloPlots", "CSV", "DataFrames"])
```

### Main Bioinformatics Software
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

## Repository Structure

The analytical codebase is modularly organized into 12 dedicated scientific analysis suites. Below is the structured catalog of pipelines, scripts, and documentation files, where each module's contents are organized in a dedicated table and all script names provide direct clickable links to the source files:

### 1. `concordance_node_label/` — Multimetric Internode Concordance Factors

| Script / File | Type | Description |
| :--- | :---: | :--- |
| [`concordance_factor_gcf_scf_qcf.sh`](concordance_node_label/concordance_factor_gcf_scf_qcf.sh) | Shell | IQ-TREE 3 gene (gCF) and site (sCF) concordance factors computation |
| [`gcf_scf.sh`](concordance_node_label/gcf_scf.sh) | Shell | Summary and parsing of gCF and sCF output tables |
| [`phytop.sh`](concordance_node_label/phytop.sh) | Shell | ASTRAL quartet concordance vector calculation ($\psi_1-\psi_4$) |
| [`ICA_IC_value.sh`](concordance_node_label/ICA_IC_value.sh) | Shell | RAxML Internode Certainty (IC and ICA) computation |
| [`ICA_plot.R`](concordance_node_label/ICA_plot.R) | R | Plotting of Internode Certainty distributions |
| [`QS_test.sh`](concordance_node_label/QS_test.sh) | Shell | Quartet Sampling calculation across species-tree nodes |

### 2. `cytonuclear_display/` — Visual Conflict and Topological Distance Projection

| Script / File | Type | Description |
| :--- | :---: | :--- |
| [`densitree.R`](cytonuclear_display/densitree.R) | R | DensiTree visual cloud overlay of gene trees |
| [`RF_path_MDS_new.R`](cytonuclear_display/RF_path_MDS_new.R) | R | Robinson-Foulds and Path distance MDS 2D projection (updated) |
| [`RF_path_MDS.R`](cytonuclear_display/RF_path_MDS.R) | R | Robinson-Foulds and Path distance MDS 2D projection |
| [`QS_visual.R`](cytonuclear_display/QS_visual.R) | R | Quartet Sampling QC/QD/QI visualization |
| [`QS_test.sh`](cytonuclear_display/QS_test.sh) | Shell | Quartet Sampling calculation across species-tree nodes |
| [`phyparts.sh`](cytonuclear_display/phyparts.sh) | Shell | PhyParts bipartition conflict analysis |
| [`phyparts_pies.sh`](cytonuclear_display/phyparts_pies.sh) | Shell | PhyParts bipartition conflict pie-chart calculation |

### 3. `factor_alignment_tree/` — Sequence & Tree Properties and Matching Experiments

| Script / File | Type | Description |
| :--- | :---: | :--- |
| [`factor.sh`](factor_alignment_tree/factor.sh) | Shell | High-throughput property extraction via AMAS & PhyKit |
| [`factor_parallel.sh`](factor_alignment_tree/factor_parallel.sh) | Shell | Parallel property extraction runner |
| [`make_tree_factor.sh`](factor_alignment_tree/make_tree_factor.sh) | Shell | Supermatrix building and ML tree inference for matched sets |
| [`random_nuclear.sh`](factor_alignment_tree/random_nuclear.sh) | Shell | Random nuclear gene sampling as background controls |
| [`factor_plot/00-utils_functions.R`](factor_alignment_tree/factor_plot/00-utils_functions.R) | R | Shared utility functions and theme definitions |
| [`factor_plot/01_data_loading.R`](factor_alignment_tree/factor_plot/01_data_loading.R) | R | Dataset ingestion and normalization |
| [`factor_plot/02_density_plots.R`](factor_alignment_tree/factor_plot/02_density_plots.R) | R | Density distribution comparisons between nuclear & plastid |
| [`factor_plot/03_correlation_analysis.R`](factor_alignment_tree/factor_plot/03_correlation_analysis.R) | R | 15-property pairwise correlation matrices |
| [`factor_plot/04_pca_analysis.R`](factor_alignment_tree/factor_plot/04_pca_analysis.R) | R | Principal Component Analysis (PCA) |
| [`factor_plot/05_rf_distance_analysis.R`](factor_alignment_tree/factor_plot/05_rf_distance_analysis.R) | R | Baseline topological distance evaluation |
| [`factor_plot/06_chloroplast_matching.R`](factor_alignment_tree/factor_plot/06_chloroplast_matching.R) | R | Property-matched subsetting experiment |
| [`factor_plot/07_final_comparison.R`](factor_alignment_tree/factor_plot/07_final_comparison.R) | R | Topological divergence reassessment |
| [`factor_plot/08_density_based_rf_analysis.R`](factor_alignment_tree/factor_plot/08_density_based_rf_analysis.R) | R | Density-stratified RF distance analysis |
| [`factor_plot/09-run_all.R`](factor_alignment_tree/factor_plot/09-run_all.R) | R | Master runner for property analysis |
| [`factor_plot/10_chloroplast_similarity_analysis.R`](factor_alignment_tree/factor_plot/10_chloroplast_similarity_analysis.R) | R | Multidimensional distance quantification |
| [`factor_plot/factor_plot.R`](factor_alignment_tree/factor_plot/factor_plot.R) | R | Monolithic plotting script |
| [`factor_plot/README.md`](factor_alignment_tree/factor_plot/README.md) | Docs | Submodule documentation |
| [`factor_plot/README_RF_ANALYSIS.md`](factor_alignment_tree/factor_plot/README_RF_ANALYSIS.md) | Docs | Documentation for topological distance comparisons |

### 4. `simulate_theta/` — Coalescent Simulation Under the MSC Model

| Script / File | Type | Description |
| :--- | :---: | :--- |
| [`get_theta.R`](simulate_theta/get_theta.R) | R | Internal node mutation parameter $\theta$ calculation |
| [`phybase_simulate.R`](simulate_theta/phybase_simulate.R) | R | Phybase coalescent simulation across $\theta$ gradient (3,000 trees) |
| [`count_four_aray.R`](simulate_theta/count_four_aray.R) | R | Quartet frequency extraction |
| [`count_four_aray_lm.R`](simulate_theta/count_four_aray_lm.R) | R | Linear regression ($R^2$) between quartet frequencies and coalescent expectation |
| [`get_number.py`](simulate_theta/get_number.py) | Python | Regression $R^2$ summary and extraction |
| [`plot.R`](simulate_theta/plot.R) | R | Regression scatter plotting |
| [`rename.sh`](simulate_theta/rename.sh) | Shell | File management for simulation replicates |
| [`RAEDME.md`](simulate_theta/RAEDME.md) | Docs | Submodule documentation |

### 5. `MSCquartets/` — MSCquartets Hypothesis Testing

| Script / File | Type | Description |
| :--- | :---: | :--- |
| [`MSCquartets.R`](MSCquartets/MSCquartets.R) | R | Hypothesis testing under T1 (tree) and T3 (coalescent) models |
| [`get_sp_ILS_IH_contribution.R`](MSCquartets/get_sp_ILS_IH_contribution.R) | R | Species-pair level ILS vs IH quantification |
| [`get_gene_ILS_IH_contribution.R`](MSCquartets/get_gene_ILS_IH_contribution.R) | R | Locus-level ILS vs IH scoring |
| [`analy_df_num.R`](MSCquartets/analy_df_num.R) | R | Hypothesis test classification summarizer |
| [`results_classify.sh`](MSCquartets/results_classify.sh) | Shell | Automation script for quartet partitioning |

### 6. `Get_ILS_IH_gene/` — Score Visualization for ILS and IH Loci

| Script / File | Type | Description |
| :--- | :---: | :--- |
| [`get_and_plot_IH_gene.R`](Get_ILS_IH_gene/get_and_plot_IH_gene.R) | R | Boxplots and barplots of locus-level ILS/IH scores |

### 7. `QuIBL/` — QuIBL Branch-Length Mixture Modeling

| Script / File | Type | Description |
| :--- | :---: | :--- |
| [`run_quibl.sh`](QuIBL/run_quibl.sh) | Shell | Triplet branch-length distribution estimation |
| [`QuIBL_plot.R`](QuIBL/QuIBL_plot.R) | R | $\Delta\text{BIC}$ classification and introgression heatmaps |

### 8. `D-statistics/` — Genome-Wide ABBA-BABA Introgression Tests

| Script / File | Type | Description |
| :--- | :---: | :--- |
| [`msa2vcf_biallelic.py`](D-statistics/msa2vcf_biallelic.py) | Python | Extraction of biallelic SNPs from MSA to VCF format |
| [`D-statistics.sh`](D-statistics/D-statistics.sh) | Shell | Dsuite Dtrios, Fbranch, and visualization runner |
| [`run_Dsuite.py`](D-statistics/run_Dsuite.py) | Python | Automated batch execution of Dsuite |
| [`run_Dsuite_xt.py`](D-statistics/run_Dsuite_xt.py) | Python | Automated batch execution of Dsuite (extended workflow) |
| [`run_Dsuite_no_snp_thinning.py`](D-statistics/run_Dsuite_no_snp_thinning.py) | Python | Dsuite execution without SNP thinning |
| [`run_Dsutie_xt_snp_thinning.py`](D-statistics/run_Dsutie_xt_snp_thinning.py) | Python | Dsuite execution with physical distance SNP thinning |
| [`make_runbash.sh`](D-statistics/make_runbash.sh) | Shell | Shell script generator for batch runs |
| [`run.sh`](D-statistics/run.sh) | Shell | Parallel execution runner |
| [`plot_box.R`](D-statistics/plot_box.R) | R | Subfamily and inter-subfamily D-statistic boxplots |
| [`sets.txt`](D-statistics/sets.txt) | Text | Species-to-clade assignment definitions |
| [`plot_order.txt`](D-statistics/plot_order.txt) | Text | Species tree topology ordering for Dsuite heatmaps |

### 9. `PhyloNetwork/` — SNaQ Semi-Directed Phylogenetic Network Inference

| Script / File | Type | Description |
| :--- | :---: | :--- |
| [`runSNaQ_hn.jl`](PhyloNetwork/runSNaQ_hn.jl) | Julia | Julia PhyloNetworks SNaQ optimization ($h = 0$ to $6$) |
| [`PhyloNetworks.R`](PhyloNetwork/PhyloNetworks.R) | R | Network topology and inheritance probability $\gamma$ extraction |
| [`plot_final_loglikelihood.R`](PhyloNetwork/plot_final_loglikelihood.R) | R | Log-pseudolikelihood slope plot for optimal $h$ selection |
| [`plot_result_julia.R`](PhyloNetwork/plot_result_julia.R) | R | Network visualization via PhyloPlots |

### 10. `ILS_GTEE_IH/` — Integrated Variance Decomposition (GTEE vs ILS vs IH)

| Script / File | Type | Description |
| :--- | :---: | :--- |
| [`ILS_GTEE_IH.SH`](ILS_GTEE_IH/ILS_GTEE_IH.SH) | Shell | Master pipeline runner |
| [`get_node_inf.R`](ILS_GTEE_IH/get_node_inf.R) | R | Node ID and reticulation index extraction |
| [`get_theta_value.R`](ILS_GTEE_IH/get_theta_value.R) | R | Per-node mutation parameter $\theta$ estimation |
| [`README_ILS_IH_GTEE.md`](ILS_GTEE_IH/README_ILS_IH_GTEE.md) | Docs | Pipeline design and module documentation |
| [`new_use_in_globel/relaimpo.R`](ILS_GTEE_IH/new_use_in_globel/relaimpo.R) | R | Hierarchical relative importance (LMG) with LOESS 95% CI |
| [`new_use_in_globel/extract_node_values_theta.py`](ILS_GTEE_IH/new_use_in_globel/extract_node_values_theta.py) | Python | Node feature table construction (with $\theta$) |
| [`new_use_in_globel/extract_node_values_without_theta.py`](ILS_GTEE_IH/new_use_in_globel/extract_node_values_without_theta.py) | Python | Node feature table construction (without $\theta$) |
| [`new_use_in_globel/MSC_geneTr_simulator.R`](ILS_GTEE_IH/new_use_in_globel/MSC_geneTr_simulator.R) | R | Gene tree simulator under the coalescent |
| [`new_use_in_globel/simulator_tripletCounter_tripletMapper.sh`](ILS_GTEE_IH/new_use_in_globel/simulator_tripletCounter_tripletMapper.sh) | Shell | Triplet enumeration pipeline |
| [`new_use_in_globel/triplet_mapper.py`](ILS_GTEE_IH/new_use_in_globel/triplet_mapper.py) | Python | Empirical vs simulated triplet mapper |
| [`new_use_in_globel/triple_frequency_counter.py`](ILS_GTEE_IH/new_use_in_globel/triple_frequency_counter.py) | Python | Triplet topology counter |
| [`new_use_in_globel/find_unbalanced_triplets.py`](ILS_GTEE_IH/new_use_in_globel/find_unbalanced_triplets.py) | Python | Identification of significantly unbalanced triplets |
| [`new_use_in_globel/clean_supp.R`](ILS_GTEE_IH/new_use_in_globel/clean_supp.R) | R | Tree support cleanup utility |
| [`new_use_in_globel/make_treelength_gt_0.R`](ILS_GTEE_IH/new_use_in_globel/make_treelength_gt_0.R) | R | Branch length zero-elimination utility |
| [`new_use_in_GTEE_part/generate_seqgen_commands.py`](ILS_GTEE_IH/new_use_in_GTEE_part/generate_seqgen_commands.py) | Python | Seq-Gen simulation command generator |
| [`new_use_in_GTEE_part/seqgen_commands.sh`](ILS_GTEE_IH/new_use_in_GTEE_part/seqgen_commands.sh) | Shell | Batch execution of Seq-Gen simulations |
| [`new_use_in_GTEE_part/clip_missing_taxa.sh`](ILS_GTEE_IH/new_use_in_GTEE_part/clip_missing_taxa.sh) | Shell | Alignment pruning for empirical taxon occupancy |
| [`new_use_in_GTEE_part/transform_astral_q1_result.py`](ILS_GTEE_IH/new_use_in_GTEE_part/transform_astral_q1_result.py) | Python | Summary of simulated bipartition recovery |

### 11. `Introgression_Functional_Enrichment/` — Introgressed Gene Detection & GO Enrichment

| Script / File | Type | Description |
| :--- | :---: | :--- |
| [`01_filter_high_confidence_trios_f4_quibl.R`](Introgression_Functional_Enrichment/01_filter_high_confidence_trios_f4_quibl.R) | R | Dual filtering (top decile $f_4$-ratio & QuIBL mixprop2) |
| [`02_eggnog_orgdb_clusterprofiler_enrichment.R`](Introgression_Functional_Enrichment/02_eggnog_orgdb_clusterprofiler_enrichment.R) | R | EggNOG parsing, custom OrgDb, and clusterProfiler |
| [`03_plot_GO_multi_threshold_barplots.R`](Introgression_Functional_Enrichment/03_plot_GO_multi_threshold_barplots.R) | R | Multi-threshold (top 3% & top 5% \|$f_{dM}$\|) GO barplots |
| [`04_analyze_GO_venn_six_groups_exact_intersections.R`](Introgression_Functional_Enrichment/04_analyze_GO_venn_six_groups_exact_intersections.R) | R | Exact 6-trio Venn intersections |
| [`05_batch_assemble_trio_figures.py`](Introgression_Functional_Enrichment/05_batch_assemble_trio_figures.py) | Python | Publication figure assembly (scatter, tree, GO) |
| [`README_Introgression_GO.md`](Introgression_Functional_Enrichment/README_Introgression_GO.md) | Docs | Workflow documentation |

### 12. `workflow/` — Analytical Pipeline Overview

| File | Type | Description |
| :--- | :---: | :--- |
| [`Rosaceae_Cytonuclear_Workflow.jpg`](workflow/Rosaceae_Cytonuclear_Workflow.jpg) | Image | Full high-resolution analytical workflow diagram |

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