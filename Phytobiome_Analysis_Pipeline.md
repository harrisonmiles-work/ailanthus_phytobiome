# Phytobiome Analysis Pipeline Documentation

**Project:** Verticillium nonalfalfae effects on Ailanthus altissima endophyte communities  
**Date:** April 2026  
**Authors:** H.H. Miles, C.J. Fearer  

## Project Overview

This project examined bacterial (16S rRNA) and fungal (ITS) endophyte communities in Ailanthus altissima before and after inoculation with biocontrol agent Verticillium nonalfalfae VnAa140. The study used a paired before-after experimental design with mixed-effects models to account for tree-level random effects.

**Key Finding:** VnAa140 treatment significantly reduced bacterial phylogenetic diversity (Faith's PD) with a large effect size (η² = 0.718), suggesting the biocontrol agent affects phylogenetic breadth while preserving taxonomic diversity.

## Experimental Design

- **Treatments:** Control vs VnAa140 (n=3 trees each)
- **Timepoints:** Pre-treatment and Post-treatment (30 days)
- **Sample structure:** Paired samples from individual trees 
- **Sequencing targets:** 16S V5-V7 (bacteria) and ITS1F-ITS2 (fungi)
- **Statistical approach:** Linear mixed-effects models with Tree as random effect

## Data Processing Pipeline

### 1. QIIME2 Processing
- **Platform:** Illumina NovaSeq X Plus 2×150
- **Denoising:** DADA2 for ASV inference
- **Taxonomy:** SILVA 138.2 (16S), UNITE v10.0 (ITS)
- **Rarefaction:** 350K reads (16S), 670K reads (ITS)
- **Tree construction:** Rooted phylogenetic trees for Faith's PD

### 2. Key QIIME2 Exports
Essential files for R analysis:
```
qiime2_exports/
├── sample-metadata-updated.tsv
├── 16S/_export/
│   ├── feature-table.tsv
│   ├── taxonomy.tsv
│   └── tree.nwk
└── ITS/_export/
    ├── feature-table.tsv
    ├── taxonomy.tsv
    └── tree.nwk
```

### 3. R Analysis Workflow

#### Required Packages
```r
# CRAN packages
library(ape)         # phylogenetic trees
library(dplyr)       # data manipulation
library(ggplot2)     # plotting
library(tidyr)       # data reshaping
library(readr)       # file reading
library(vegan)       # ecological statistics
library(picante)     # phylogenetic diversity
library(lme4)        # mixed-effects models
library(phangorn)    # phylogenetic analysis
library(emmeans)     # post-hoc comparisons
library(lmerTest)    # p-values for mixed models
library(patchwork)   # plot composition

# Bioconductor
library(phyloseq)    # microbiome data handling
```

#### Key Analysis Steps
1. **Data Import:** Load QIIME2 exports into phyloseq objects
2. **Quality Control:** Filter samples (>1000 reads) and ASVs (>1 sample)
3. **Diversity Analysis:** Calculate Faith's PD, Shannon, Simpson
4. **Statistical Testing:** Mixed-effects models with Tree as random effect
5. **Community Analysis:** PERMANOVA with UniFrac distances
6. **Visualization:** NMDS ordination and diversity plots

## Statistical Framework

### Primary Model Structure
```r
# Mixed-effects model accounting for paired design
modelPD <- lmer(PD ~ Treatment * Timepoint + (1 | Tree), data = div2)
```

### Key Statistical Results
- **Primary outcome:** Faith's Phylogenetic Diversity significantly reduced by VnAa140
- **Effect size:** η² = 0.718 (large effect - 72% variance explained)
- **Statistical significance:** F = 6.20, p = 0.037
- **Pattern:** Treatment inhibited natural temporal increase in phylogenetic diversity

### Important Statistical Insights
1. **Main effect vs. Interaction:** Significant treatment main effect without significant interaction indicates consistent treatment effect across timepoints
2. **No baseline differences:** Pre-treatment groups were not significantly different (good experimental design)
3. **Effect size calculation:** Use Sum of Squares for treatment effect [1] not interaction [3]

## File Management

### Essential Files for Reproduction (8 files total)
1. **Analysis code:** `phytobiome_analysis_final.Rmd`
2. **Metadata:** `sample-metadata-updated.tsv`  
3. **16S data:** `feature-table.tsv`, `taxonomy.tsv`, `tree.nwk`
4. **ITS data:** `feature-table.tsv`, `taxonomy.tsv`, `tree.nwk`

### Files NOT Needed for Publication
- Raw FASTQ files (uploaded to NCBI SRA)
- 95+ intermediate QIIME2 files (.qza/.qzv)
- BIOM format files (TSV versions used)
- Database reference files

## Key Code Patterns

### Data Loading and Phyloseq Object Creation
```r
# Load data with proper handling of QIIME2 exports
data_16s <- read.table("../qiime2_exports/16S/_export/feature-table.tsv",
                      header = TRUE, sep = "\t", row.names = 1,
                      comment.char = "", check.names = FALSE,
                      quote = "", skip = 1)

# Process taxonomy strings and remove QIIME2 prefixes
tax_strings_16s <- taxa_16s[["Taxon"]]
parts_16s <- strsplit(tax_strings_16s, ";\\s*")
# Remove k__, p__, c__ prefixes
x <- sub("^[dkpcofgs]__", "", x)
```

### Tree Preprocessing
```r
# Ensure phylogenetic tree is properly formatted
tr_16s <- ape::reorder.phylo(tr_16s, "postorder")
tr_16s <- ape::multi2di(tr_16s)  # Convert to bifurcating
if (!ape::is.rooted(tr_16s)) {
  tr_16s <- phangorn::midpoint(tr_16s)
}
# Fix zero/negative branch lengths
tr_16s$edge.length[is.na(tr_16s$edge.length) | tr_16s$edge.length <= 0] <- 1e-8
```

### Statistical Analysis
```r
# Calculate Faith's PD with proper tree handling
comm.pd_16s <- pd(t(as.data.frame(otu_table(PO_16s))), phy_tree(PO_16s), include.root = TRUE)

# Mixed-effects model
modelPD_16s <- lmer(PD ~ Treatment * Timepoint + (1 | Tree), data = div2_16s)

# Effect size calculation (for significant effects)
if(anova_pd_16s$`Pr(>F)`[1] < 0.05) {
  ss_total <- sum(anova_pd_16s$`Sum Sq`)
  eta_sq <- anova_pd_16s$`Sum Sq`[1] / ss_total  # [1] = treatment main effect
  cat("Effect size for Treatment main effect (η²):", round(eta_sq, 3), "\n")
}
```

### PERMANOVA with Blocking
```r
# Account for paired design in community analysis
set.seed(100)
perm_16s <- how(nperm = 999)
setBlocks(perm_16s) <- sampledf_16s$Tree  # Block by tree
pmv_16s <- adonis2(dist.mat_16s ~ Treatment * Timepoint,
                   data = sampledf_16s, permutations = perm_16s)
```

## Common Issues and Solutions

### 1. Package Dependencies
**Issue:** Some packages loaded but not used  
**Solution:** Audit actual package usage with grep searches, remove unused libraries

### 2. Effect Size Calculation
**Issue:** Original code calculated η² for interaction (which wasn't significant)  
**Solution:** Calculate η² for significant main effects using correct index

### 3. Tree Processing
**Issue:** Phylogenetic trees need specific formatting for Faith's PD  
**Solution:** Apply reordering, rooting, bifurcation, and branch length corrections

### 4. Statistical Interpretation
**Issue:** Confusion between main effects and interactions  
**Solution:** Understand that significant main effect without interaction = consistent treatment effect

## Data Submission Checklist

### Public Repository
- ✅ Raw sequencing data: NCBI SRA PRJNA1449437
- ✅ Sample accessions: SRX32807467-SRX32807478

### Journal Submission Package
- ✅ Manuscript with proper statistical reporting
- ✅ High-resolution figures (PNG format, >150 DPI)
- ✅ Analysis code (R Markdown)
- ✅ Essential processed data files
- ✅ Phylogenetic trees for reproducibility

## Best Practices Learned

1. **Always include Tree ID as random effect** for paired sampling designs
2. **Report effect sizes (η²)** alongside statistical significance
3. **Use letter notation (a,b)** rather than asterisks for significance marking
4. **Verify baseline similarity** before claiming treatment effects
5. **Keep minimal file set** for publication (avoid overwhelming reviewers)
6. **Document rarefaction depths** clearly (unusually high depths need justification)
7. **Include accession numbers** in Methods section for data availability

## Future Pipeline Applications

This pipeline can be adapted for:
- Other biocontrol-microbiome interaction studies
- Paired before-after treatment designs in microbiology
- Any phyloseq-based microbiome analysis requiring phylogenetic diversity
- Studies needing mixed-effects models for non-independent sampling

## Software Versions Used
- QIIME2: 2025.7
- R: 4.5.2
- Key packages: phyloseq (1.54.2), lme4 (2.0-1), picante (1.8.2), vegan (2.7-3)

## Contact
Harrison H. Miles (harrisonmiles@vt.edu)  
Virginia Tech Department of Forest Resources and Environmental Conservation

---
*Pipeline documented April 2026 for future reference and reproducibility*