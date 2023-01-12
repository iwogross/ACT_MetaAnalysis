# ACT_MetaAnalysis

Version controlled and editable source for the data and code supporting the paper "The fitness consequences of wildlife conservation translocations: A meta-analysis:" by Gross, Wilson, & Wolak.

## Data

### Data citation
If you use the data, please cite the data package TBD

### Data metadata

Column header descriptions for the primary dataset analyzed in the paper.

This systematic review and meta-analysis examines the efficacy of animal conservation translocations using direct comparisons of animal performance between translocated and wild-resident organisms. Each data row characterizes a single comparison between a translocated and wild-conspecific cohort in the wild. Comparisons are standardized using various effect size estimators. 

Column headings for the datasets reflect variables defined and discussed in the main manuscript and supplementary materials. In brief, these are:

 - `estimateID`: Unique comparison (row) indentifier
 - `estimator`: Designated effect size estimator (standardized mean difference [SMD] or log odds ratio [LOR])
 - `studyID`: Unique study identifier
 - `groupID`: Within-study identifier of unique animal cohorts
 - `sex`: Designated sex of compared cohorts (female-only, mixed-sex, male-only)
 - `common.name`: Layman's species moniker
 - `phylum`: Taxonomic phylum designation
 - `class`: Taxonomic class designation
 - `order`: Taxonomic order designation         
 - `fam`: Taxonomic family designation
 - `genus`: Taxonomic genus designation
 - `species`: Taxonomic species designation
 - `enrich`: Ordinal scale of within-study enrichment strategies provided to released cohort (simplified to binomial `enrich.bin` term below)
 - `capt.gens`: Number of generations of captive ancestry of released cohort (simplified to trinomial `capt.bin` term below)
 - `scope`: Fitness correlate or component measured as a performance metric (e.g., constitution, survival, reproduction)
 - `strategy`: _in situ_ vs. _ex situ_ translocation strategy
 - `n1i`: Wild cohort sample size (SMD)
 - `m1i`: Wild cohort effect mean (SMD)
 - `sd1i`: Wild cohort effect standard deviation (SMD)
- `n1i`: Translocated cohort sample size (SMD)
 - `m1i`: Translocated cohort effect mean (SMD)
 - `sd1i`: Translocated cohort effect standard deviation (SMD)
 - `A`: Wild cohort positive outcomes (LOR)
 - `C`: Translocated cohort positive outcomes (LOR)
 - `B`: Wild cohort negative outcomes (LOR)
 - `D`: Translocated cohort negative outcomes (LOR)
 - `PrintYear`: Publication year
 - `enrich.present`: Study-specific presence/absence of enrichment
 - `capt.present`: Study-specific presence/absence of captive ancestry
 - `enrich.bin`: Presence vs. absence of ernichmen
 - `capt.bin`: Categorical descriptor of the number of captive generations (0 gens, < 1 gen, >1 gen)

## Changes
For ease of reference, an overview of significant changes to be noted below. Tag with commits or issues, where appropriate.