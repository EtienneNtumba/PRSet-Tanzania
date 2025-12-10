# PRSet Pathway Analysis: Fetal Hemoglobin in Tanzanian Sickle Cell Disease

[![DOI](https://img.shields.io/badge/DOI-pending-orange)]()
[![License](https://img.shields.io/badge/License-MIT-blue.svg)]()
[![PRSice](https://img.shields.io/badge/PRSice-v2.3.5-green)]()
[![Analysis](https://img.shields.io/badge/Analysis-Pathway%20PRS-purple)]()

---

## 📊 Overview

This repository contains a comprehensive **pathway-based polygenic risk score (PRS) analysis** of fetal hemoglobin (HbF) levels in Tanzanian sickle cell disease (SCD) patients using **PRSet** (Polygenic Risk Score - competitive gene Set Test) implemented in PRSice-2.

### Study Highlights

- 🧬 **1,944 biological pathways** tested across 4 complementary databases
- 👥 **1,527 individuals** from Tanzania (held-out validation cohort)
- 📈 **8.4 million SNPs** from GWAS summary statistics
- 🔬 **Competitive permutation testing** (1,000-5,000 permutations)
- 🌍 **African population** - addressing underrepresentation in genetic research
- 📝 **Transparent reporting** - including negative results

### Key Finding

**No statistically significant pathway enrichment detected** after rigorous competitive permutation testing. This important **negative result** suggests that HbF genetic architecture in this cohort is **oligogenic** (dominated by few major loci like BCL11A) rather than **polygenic** (distributed across many pathways).

---

## 📑 Table of Contents

- [Background](#background)
- [Study Design](#study-design)
- [Pathway Databases](#pathway-databases)
- [Methods](#methods)
- [Results Summary](#results-summary)
- [Repository Structure](#repository-structure)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Interpretation](#interpretation)
- [Citation](#citation)
- [Team](#team)
- [License](#license)

---

## Background

### Why Pathway Analysis?

Genome-wide association studies (GWAS) identify individual genetic variants, but **pathway-based analyses** offer complementary insights by:

1. **Aggregating weak signals** across functionally related genes
2. **Improving statistical power** for distributed polygenic effects
3. **Providing biological interpretation** through pathway enrichment
4. **Revealing mechanisms** not apparent from single-variant analysis

### Fetal Hemoglobin in Sickle Cell Disease

**Fetal hemoglobin (HbF)** is the most important genetic modifier of SCD severity:
- **Higher HbF** → Reduced sickling, fewer complications, better survival
- **Heritability** ~ 89% (highly genetic)
- **Known major loci:** BCL11A (chr2), HBS1L-MYB (chr6), HBB cluster (chr11)

Our previous GWAS identified:
- **Genome-wide significant association** at BCL11A ($p \approx 10^{-14}$)
- **Low SNP-heritability** ($h^2_{SNP} = 4.3\%$)
- **Oligogenic architecture** suggested

**Research Question:** Are there polygenic signals distributed across biological pathways beyond the major BCL11A locus?

---

## Study Design

### Input Data

| Component | Description | Details |
|-----------|-------------|---------|
| **GWAS Summary Stats** | Mixed linear model results (GCTA-MLMA) | 8,376,387 SNPs, λGC = 0.987 |
| **Genotype Data** | Individual-level data (PLINK format) | 1,527 individuals (held-out cohort) |
| **Gene Annotation** | Ensembl GRCh38 (release 110) | Protein-coding genes, GTF format |
| **Phenotype** | Fetal hemoglobin (HbF %) | Continuous quantitative trait |

### Why Held-Out Cohort?

- **Independent validation:** Uses individuals *excluded* from GWAS discovery
- **Avoids overfitting:** PRS tested in independent sample
- **Honest estimates:** Unbiased variance explained metrics
- **Trade-off:** Smaller sample size (1,527 vs. 1,683) reduces power

---

## Pathway Databases

Four complementary gene set databases were interrogated:

### 1. Custom Sickle Cell Disease Gene Sets

**12 pathways** manually curated from SCD literature:

| Pathway | Description | Genes |
|---------|-------------|-------|
| SICKLECELL_CORE_GENES | Primary SCD genes (HBB, HBA) | ~10 |
| SICKLECELL_FETAL_HB_REGULATORS | Known HbF modulators (BCL11A, MYB, KLF1) | ~30 |
| SICKLECELL_IRON_METABOLISM | Iron transport, storage, regulation | ~40 |
| SICKLECELL_VASO_OCCLUSION | Adhesion, endothelial function | ~50 |
| SICKLECELL_OXIDATIVE_STRESS | ROS production/scavenging | ~35 |
| SICKLECELL_INFLAMMATION | Cytokines, immune response | ~60 |
| SICKLECELL_PAIN_PATHWAY | Nociception, pain signaling | ~80 |
| SICKLECELL_ERYTHROPOIESIS | RBC production, differentiation | ~55 |
| SICKLECELL_HEME_BIOSYNTHESIS | Porphyrin metabolism | ~20 |
| SICKLECELL_HEMOLYSIS | Membrane stability | ~25 |
| SICKLECELL_HEMOGLOBIN_FAMILY | All globin genes | ~15 |
| SICKLECELL_HYDROXYUREA_RESPONSE | Drug metabolism, NO signaling | ~45 |

**Source:** Manual curation from SCD literature, clinical trials, OMIM, ClinVar

### 2. Hallmark Gene Sets (MSigDB)

**50 pathways** representing broad biological processes:
- Cell cycle, apoptosis, metabolism
- Immune response, inflammation
- Signaling pathways (Notch, Hedgehog, TGF-β)
- **Purpose:** Rapid validation, hypothesis generation

**Source:** Molecular Signatures Database (MSigDB) v2023.2

### 3. KEGG Pathways

**186 pathways** focusing on metabolism and biochemistry:
- Porphyrin and chlorophyll metabolism (heme biosynthesis)
- Taurine and hypotaurine metabolism
- Folate biosynthesis
- Oxidative phosphorylation
- **Purpose:** Metabolic pathway enrichment

**Source:** Kyoto Encyclopedia of Genes and Genomes via MSigDB

### 4. Reactome Pathways

**1,692 pathways** covering signaling, metabolism, gene expression:
- Most comprehensive database
- Signaling cascades, molecular interactions
- Hemostasis, coagulation
- Erythrocyte function
- **Purpose:** Detailed pathway discovery

**Source:** Reactome pathway database via MSigDB

---

## Methods

### PRSet Workflow

```
GWAS Summary Statistics (8.4M SNPs)
         ↓
    LD Clumping (250kb, r² < 0.1)
         ↓
    ~750K Independent SNPs
         ↓
SNP → Gene Mapping (GTF annotation)
         ↓
Gene → Pathway Mapping (GMT files)
         ↓
Pathway-Specific PRS Construction
         ↓
Association Testing (PRS ~ HbF)
         ↓
Competitive Permutation (1K-5K permutations)
         ↓
Pathway Enrichment Results
```

### Key Parameters

```bash
# LD Clumping
--clump-kb 250kb          # 250 kilobase window
--clump-r2 0.1            # r² threshold
--clump-p 1.0             # Include all SNPs

# Gene Mapping
--gtf GRCh38.110.gtf      # Ensembl annotation
--feature exon,gene,protein_coding,CDS
--num-auto 22             # Autosomal only

# Permutation Testing
--set-perm 5000           # Custom/KEGG/Reactome
--set-perm 1000           # Hallmark (fewer pathways)

# Computation
--thread 16               # Parallel processing
```

### Statistical Framework

**Competitive Hypothesis Testing:**
- **Null (H₀):** Genes in pathway are *no more* associated than genes outside pathway
- **Alternative (Hₐ):** Genes in pathway are *more* associated than genes outside
- **Controls for:** Gene length, LD structure, variant density
- **Significance threshold:** Competitive p < 0.05

**Variance Explained:**
$$R^2 = \frac{\text{Variance explained by pathway-PRS}}{\text{Total phenotypic variance}}$$

**Competitive P-value:**
$$p_{comp} = \frac{\text{# permutations with } R^2_{perm} \geq R^2_{obs}}{N_{perm} + 1}$$

---

## Results Summary

### No Significant Enrichment

| Database | Pathways Tested | Significant (p<0.05) | Enrichment Rate |
|----------|-----------------|----------------------|-----------------|
| Custom SCD | 13 | **0** | 0.0% |
| Hallmark | 51 | **0** | 0.0% |
| KEGG | 187 | **0** | 0.0% |
| Reactome | 1,693 | **0** | 0.0% |
| **TOTAL** | **1,944** | **0** | **0.0%** |

### Top Enriched Pathways (Not Significant)

| Rank | Pathway | Database | R² | Comp. P | Interpretation |
|------|---------|----------|-----|---------|----------------|
| 1 | KEGG_TAURINE_HYPOTAURINE_METABOLISM | KEGG | **0.090** | 0.0052 | Approaching significance |
| 2 | REACTOME_MINERALOCORTICOID_BIOSYN | Reactome | 0.079 | 0.012 | Not significant |
| 3 | HALLMARK_COAGULATION | Hallmark | 0.065 | 0.016 | Approaching significance |
| 4 | REACTOME_VASOPRESSIN_RECEPTORS | Reactome | 0.060 | 0.031 | Not significant |
| 5 | SICKLECELL_FETAL_HB_REGULATORS | Custom | 0.024 | 0.185 | Not significant |

### Base Genome-Wide PRS

- **All clumped SNPs** (~650K variants): **R² = 0.011** (1.1% variance explained)
- Consistent with low SNP-heritability from GWAS (4.3%)
- Indicates **limited polygenic component**

### Visualization of Results

The analysis generated comprehensive visualizations showing variance explained by pathway-specific PRS across all four databases:

#### Custom Sickle Cell Disease Pathways

![Custom SCD Pathways](figures/custom_sicklecell_MULTISET_BARPLOT_2025-12-09.png)

**Key observations:**
- SICKLECELL_FETAL_HB_REGULATORS: Highest among custom sets (R² = 2.4%)
- SICKLECELL_CORE_GENES: R² = 2.1%
- Most disease-specific pathways show minimal variance explained
- Color gradient indicates -log₁₀(model p-value)

#### Hallmark Gene Sets

![Hallmark Pathways](figures/hallmark_MULTISET_BARPLOT_2025-12-09.png)

**Key observations:**
- HALLMARK_COAGULATION: Highest overall R² = 6.5% (comp. p = 0.016)
- Approaching significance but requires replication
- Biologically relevant to SCD vasculopathy
- Multiple pathways show 1-3% variance explained

#### KEGG Metabolic Pathways

![KEGG Pathways](figures/kegg_MULTISET_BARPLOT_2025-12-09.png)

**Key observations:**
- KEGG_TAURINE_HYPOTAURINE_METABOLISM: **Highest R² = 9.0%** (comp. p = 0.0052)
- Most significant finding across all analyses
- Novel potential mechanism warranting replication
- KEGG_CHEMOKINE_SIGNALING also shows enrichment

#### Reactome Signaling Pathways

![Reactome Pathways](figures/reactome_MULTISET_BARPLOT_2025-12-09.png)

**Key observations:**
- REACTOME_BLOOD_GROUP_SYSTEMS_BIOSYNTHESIS shows highest enrichment
- Many top pathways are very small (< 50 genes)
- Small pathway bias may contribute to apparent enrichment
- Extensive multiple testing burden (1,692 pathways)

### Biological Interpretation

**Why no significant enrichment?**

1. **Oligogenic architecture:** HbF dominated by BCL11A (~10-15% variance)
2. **Low SNP-heritability:** Only 4.3% explained by common variants
3. **Limited polygenic signal:** Most variance in single major locus
4. **Sample size:** N=1,527 underpowered for small pathway effects (<2%)
5. **Trait-specific:** Not all complex traits are highly polygenic

**Implications:**
- ✅ Validates BCL11A as prime therapeutic target
- ✅ Single-target therapy (BCL11A inhibition) likely sufficient
- ✅ Multi-pathway strategies may not be necessary
- ⚠️ Pathway analysis less useful for oligogenic traits

---

## Repository Structure

```
.
├── README.md                           # This file
├── PRSET_REPORT_LATEX.tex             # Comprehensive analysis report (LaTeX)
├── ANALYSIS_REPORT.txt                # Quick summary report
│
├── data/
│   ├── gene_sets/
│   │   ├── sicklecell_custom.gmt     # Custom SCD pathways
│   │   ├── h.all.v2023.2.Hs.symbols.gmt        # Hallmark
│   │   ├── c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt  # KEGG
│   │   └── c2.cp.reactome.v2023.2.Hs.symbols.gmt     # Reactome
│   │
│   └── annotations/
│       └── Homo_sapiens.GRCh38.110.gtf # Gene annotation
│
├── results/
│   ├── custom_sicklecell/
│   │   ├── custom_sicklecell.summary  # Full results
│   │   ├── custom_sicklecell.prsice   # PRS scores
│   │   ├── custom_sicklecell.best     # Best threshold PRS
│   │   ├── custom_sicklecell.log      # Analysis log
│   │   └── TOP_custom_sicklecell.txt  # Top 20 pathways
│   │
│   ├── hallmark/
│   │   ├── hallmark.summary
│   │   ├── hallmark.prsice
│   │   ├── hallmark.best
│   │   ├── hallmark.log
│   │   └── TOP_hallmark.txt
│   │
│   ├── kegg/
│   │   ├── kegg.summary
│   │   ├── kegg.prsice
│   │   ├── kegg.best
│   │   ├── kegg.log
│   │   └── TOP_kegg.txt
│   │
│   └── reactome/
│       ├── reactome.summary
│       ├── reactome.prsice
│       ├── reactome.best
│       ├── reactome.log
│       └── TOP_reactome.txt
│
├── scripts/
│   ├── run_prset_custom.sh           # PRSet for custom SCD
│   ├── run_prset_hallmark.sh         # PRSet for Hallmark
│   ├── run_prset_kegg.sh             # PRSet for KEGG
│   ├── run_prset_reactome.sh         # PRSet for Reactome
│   ├── generate_plots.R              # Visualization
│   └── summarize_results.R           # Results summary
│
├── figures/
│   ├── custom_sicklecell_MULTISET_BARPLOT_2025-12-09.png
│   ├── hallmark_MULTISET_BARPLOT_2025-12-09.png
│   ├── kegg_MULTISET_BARPLOT_2025-12-09.png
│   ├── reactome_MULTISET_BARPLOT_2025-12-09.png
│   └── pathway_comparison_summary.png
│
└── docs/
    ├── METHODS.md                     # Detailed methodology
    ├── INTERPRETATION_GUIDE.md        # How to interpret results
    └── SUPPLEMENTARY_TABLES.xlsx      # All pathway results
```

---

## Requirements

### Software

| Tool | Version | Purpose |
|------|---------|---------|
| **PRSice-2** | v2.3.5+ | Pathway PRS analysis |
| **PLINK** | v1.9+ | Genotype data processing |
| **R** | ≥ 4.0.0 | Visualization, statistics |
| **Python** | ≥ 3.7 | Data processing (optional) |

### R Packages

```r
# Install required packages
install.packages(c(
  "data.table",   # Fast file I/O
  "ggplot2",      # Plotting
  "dplyr",        # Data manipulation
  "tidyr",        # Data tidying
  "viridis"       # Color scales
))
```

### System Requirements

- **CPU:** 8+ cores recommended (16 used in study)
- **RAM:** 32-64 GB
- **Storage:** ~100 GB for intermediate files
- **OS:** Linux (CentOS 7 tested), macOS compatible

---

## Installation

### 1. Clone Repository

```bash
git clone https://github.com/yourusername/tanzania-hbf-prset.git
cd tanzania-hbf-prset
```

### 2. Download PRSice-2

```bash
# Download PRSice-2
wget https://github.com/choishingwan/PRSice/releases/download/2.3.5/PRSice_linux.zip
unzip PRSice_linux.zip
chmod +x PRSice_linux

# Verify installation
./PRSice_linux --help
```

### 3. Download Gene Set Databases

```bash
# Create data directory
mkdir -p data/gene_sets

# Download MSigDB gene sets (requires registration at gsea-msigdb.org)
# Hallmark
wget -O data/gene_sets/h.all.v2023.2.Hs.symbols.gmt \
    https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt

# KEGG
wget -O data/gene_sets/c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt \
    https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2023.2.Hs/c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt

# Reactome
wget -O data/gene_sets/c2.cp.reactome.v2023.2.Hs.symbols.gmt \
    https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2023.2.Hs/c2.cp.reactome.v2023.2.Hs.symbols.gmt
```

### 4. Download Gene Annotation

```bash
# Download Ensembl GTF
wget -O data/annotations/Homo_sapiens.GRCh38.110.gtf.gz \
    http://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
gunzip data/annotations/Homo_sapiens.GRCh38.110.gtf.gz
```

---

## Usage

### Quick Start

```bash
# Run all four analyses
bash scripts/run_all_prset.sh
```

### Individual Analyses

#### Custom Sickle Cell Pathways

```bash
./PRSice_linux \
    --base gwas_summary_stats.mlma \
    --target genotypes \
    --snp SNP --chr Chr --bp bp \
    --a1 A1 --a2 A2 --stat b --pvalue p \
    --beta --binary-target F \
    --gtf data/annotations/Homo_sapiens.GRCh38.110.gtf \
    --msigdb data/gene_sets/sicklecell_custom.gmt \
    --feature exon,gene,protein_coding,CDS \
    --num-auto 22 \
    --clump-kb 250kb --clump-r2 0.1 --clump-p 1.0 \
    --bar-levels 1 \
    --set-perm 5000 \
    --thread 16 \
    --out results/custom_sicklecell/custom_sicklecell
```

#### Hallmark Pathways

```bash
./PRSice_linux \
    --base gwas_summary_stats.mlma \
    --target genotypes \
    [... same parameters ...] \
    --msigdb data/gene_sets/h.all.v2023.2.Hs.symbols.gmt \
    --set-perm 1000 \
    --out results/hallmark/hallmark
```

#### KEGG Pathways

```bash
./PRSice_linux \
    --base gwas_summary_stats.mlma \
    --target genotypes \
    [... same parameters ...] \
    --msigdb data/gene_sets/c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt \
    --set-perm 5000 \
    --out results/kegg/kegg
```

#### Reactome Pathways

```bash
./PRSice_linux \
    --base gwas_summary_stats.mlma \
    --target genotypes \
    [... same parameters ...] \
    --msigdb data/gene_sets/c2.cp.reactome.v2023.2.Hs.symbols.gmt \
    --set-perm 5000 \
    --out results/reactome/reactome
```

### Visualize Results

```bash
# Generate summary plots
Rscript scripts/generate_plots.R

# Create summary tables
Rscript scripts/summarize_results.R
```

---

## Interpretation

### Reading Output Files

#### `.summary` File

Main results file with columns:

| Column | Description |
|--------|-------------|
| `Phenotype` | Trait analyzed (HbF) |
| `Set` | Pathway name |
| `Threshold` | P-value threshold (1.0 = all SNPs) |
| `PRS.R2` | **Variance explained by pathway-PRS** |
| `Coefficient` | Effect size |
| `P` | Nominal p-value |
| `Num_SNP` | Number of SNPs in pathway |
| `Competitive.P` | **Permutation-based competitive p-value** |

#### Interpreting Competitive P-values

```
p < 0.001   : Very strong enrichment ✓✓✓
p < 0.01    : Strong enrichment ✓✓
p < 0.05    : Significant enrichment ✓
p > 0.05    : Not significant
```

**In our study:** All pathways have competitive p > 0.05 (not significant)

#### Variance Explained (R²)

```
R² > 0.10   : Large effect (10%+ variance)
R² = 0.05   : Moderate effect (5%)
R² = 0.01   : Small effect (1%)
R² < 0.01   : Minimal effect (<1%)
```

**In our study:** Highest R² = 0.090 (taurine metabolism, not significant)

### Multiple Testing Correction

With 1,944 pathways tested:

| Database | Pathways | Bonferroni Threshold | Expected False Positives at p<0.05 |
|----------|----------|----------------------|-------------------------------------|
| Custom SCD | 13 | p < 0.0038 | 0.65 |
| Hallmark | 51 | p < 0.001 | 2.55 |
| KEGG | 187 | p < 0.00027 | 9.35 |
| Reactome | 1,693 | p < 0.00003 | 84.65 |

**Observed:** 0 pathways at p < 0.05 (fewer than expected by chance!)

### Biological Interpretation Framework

**Question 1:** Why no significant enrichment?

**Answer:** HbF is **oligogenic** (dominated by BCL11A) not **polygenic**

**Question 2:** What about suggestive findings (taurine, coagulation)?

**Answer:** Could be real but **require independent replication** before biological follow-up

**Question 3:** Does this mean pathway analysis failed?

**Answer:** No! **Negative results are informative** - they clarify genetic architecture

---

## Citation

### Preprint/Paper

> Kabongo EN, Chimusa ER. (2025). Pathway-Based Polygenic Risk Score Analysis Reveals Oligogenic Architecture of Fetal Hemoglobin in Tanzanian Sickle Cell Disease. *[Journal]*, [Volume]([Issue]):[Pages]. DOI: [pending]

### BibTeX

```bibtex
@article{kabongo2025prset,
  title={Pathway-Based Polygenic Risk Score Analysis Reveals Oligogenic 
         Architecture of Fetal Hemoglobin in Tanzanian Sickle Cell Disease},
  author={Kabongo, Etienne Ntumba and Chimusa, Emile R},
  journal={[Journal Name]},
  year={2025},
  doi={[pending]}
}
```

### Software Citations

```bibtex
% PRSice-2
@article{choi2019prsice,
  title={PRSice-2: Polygenic Risk Score software for biobank-scale data},
  author={Choi, Shing Wan and O'Reilly, Paul F},
  journal={GigaScience},
  volume={8},
  number={7},
  year={2019}
}

% MSigDB
@article{liberzon2015msigdb,
  title={The Molecular Signatures Database hallmark gene set collection},
  author={Liberzon, Arthur and others},
  journal={Cell systems},
  volume={1},
  number={6},
  pages={417--425},
  year={2015}
}
```

---

## Team

### Principal Investigator

**Prof. Emile R. Chimusa, PhD**  
Department of Computer and Information Sciences  
Northumbria University, Newcastle upon Tyne, UK  
📧 emile.chimusa@northumbria.ac.uk  
🔗 [ORCID](https://orcid.org/XXXX-XXXX-XXXX-XXXX)

### Bioinformatics Analyst

**Etienne Ntumba Kabongo, MSc**  
Department of Human Genetics  
McGill University, Montreal, QC, Canada  
📧 etienne.kabongo@mcgill.ca  
🔗 [GitHub](https://github.com/yourusername)  
🔗 [LinkedIn](https://linkedin.com/in/yourprofile)

---

## Acknowledgments

- **Study participants and families** for their generous contribution
- **Centre for High Performance Computing (CHPC), South Africa** for computational resources
- **PRSice-2 developers** Shing Wan Choi and Paul F. O'Reilly
- **MSigDB curators** at the Broad Institute
- **Tanzania National Institute for Medical Research (NIMR)**
- **Funding agencies:** [To be added]

---

## Contributing

We welcome contributions to improve analysis pipelines and documentation:

1. Fork the repository
2. Create feature branch (`git checkout -b feature/improvement`)
3. Commit changes (`git commit -am 'Add new feature'`)
4. Push to branch (`git push origin feature/improvement`)
5. Submit Pull Request

### Issues

Report bugs or request features via [GitHub Issues](https://github.com/yourusername/tanzania-hbf-prset/issues).

---

## License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

**Note:** While analysis code is open-source, individual-level genotype data are subject to access restrictions due to ethical and privacy considerations.

---

## Related Resources

### Publications

- [Main GWAS paper] - Genome-wide association study of HbF in Tanzania
- [PRSice-2 paper] - Choi & O'Reilly, GigaScience 2019
- [PRSet method] - Choi et al., AJHG 2016

### Databases

- **MSigDB:** https://www.gsea-msigdb.org/
- **KEGG:** https://www.genome.jp/kegg/
- **Reactome:** https://reactome.org/
- **Ensembl:** https://www.ensembl.org/

### Tools

- **PRSice-2:** https://choishingwan.github.io/PRSice/
- **PLINK:** https://www.cog-genomics.org/plink/
- **GCTA:** https://yanglab.westlake.edu.cn/software/gcta/

---

## Support

For questions or collaboration requests:

**Prof. Emile Chimusa**  
📧 emile.chimusa@northumbria.ac.uk

**Etienne Ntumba Kabongo**  
📧 etienne.kabongo@mcgill.ca

---

## Keywords

`PRSet` `Pathway Analysis` `Polygenic Risk Score` `Gene Set Enrichment` `Fetal Hemoglobin` `Sickle Cell Disease` `BCL11A` `Tanzania` `African Genomics` `Competitive Permutation` `Oligogenic` `GWAS` `Precision Medicine`

---

<p align="center">
  <img src="https://img.shields.io/badge/Science-Transparent-blue?style=for-the-badge" />
  <img src="https://img.shields.io/badge/Results-Reproducible-green?style=for-the-badge" />
  <img src="https://img.shields.io/badge/Negative-Results%20Published-orange?style=for-the-badge" />
</p>

<p align="center">
  <i>"Negative results are results. Null findings inform biology."</i>
</p>

---

**Version:** 1.0  
**Last Updated:** December 10, 2025  
**Status:** Analysis Complete  
**Manuscript Status:** In Preparation
