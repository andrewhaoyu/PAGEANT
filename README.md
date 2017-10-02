## PAGEANT

## PAGEANT: Power Analysis for Genetic Association Tests

## Introduction

The application allows rapid power and sample size analysis for a variety of genetic association tests by specification of a few key parameters [Derkach et al. 2017]. Power and sample size calculations can be done at the level of single variant for simple trend test, at the level of a genes/regions for various complex aggregated tests [Neale et al. 2011, Derkach et al. 2013, Wu et al. 2011, Madsen and Browning 2009] and at the level of the whole genome for the assessment of overall yield of a study. The calculations currently uses underlying distribution of gene size and minor allele frequencies of variants observed in the in the public data for 60,000 individuals from Exome Aggregation Consortium [Lek et al. 2016]

## Power and sample size for association test at the level of a single variant or a single gene/region

### Essential Input Parameters


1)	EV: For continuous trait, EV represents % of phenotypic variance explained the variants in a gene (or by a single variant for single-variant tests). For binary trait, EV represents sum of squares of log-odds-ratios (in standardized unit) for the variants in a gene (or for a single variant). For example, if power calculation is desired for a locus which may have 5 causal variants each with an OR=1.2, then EV should be set to 5× (ln(1.2)^2)=0.17;

2) alpha = level of the test;

3) Sample size: Total sample size for a continuous trait or the number of cases and number of controls for a case-control study; This parameter required in power calculations.

4) Power target: targeted average power of the test


### Optional Input Parameters

1)	Total number of variants (J): The total number of variants under study within a gene/region. This is a key parameter in power calculation for the gene-level tests and when it’s not specified the application evaluates distribution of power according to distribution observed for in the ExAC database.With J=1, gene based power calculation simplifies to  single variant one;

2)	Proportion of causal variants (J_c/J) : Assumed proportion of causal variants in a locus as a ratio to the total number of variants. This parameter is required for burden test and a more accurate second-order approximation of the variance component test. For burden tests, it’s assumed that all causal variants are either deleterious or protective and by default proportion of causal variants is set to 20%;

3)Range of EV: Instead of a single EV, the user can specify a range of EV over which power calculation is desired ;

### Output

The application conducts power analysis and sample size under three different models for genetic architecture assuming (S1): MAF is independent of EV ;(S2) MAF is independent of genetic effects measured in the unit of per copy of an allele (beta^2=EV/(2MAF(1-MAF))); and (S3) MAF is negatively correlated with genetic effect through the function β=−log10(MAF). When a single EV is specified, for each genetic architecture, it returns a distribution of power or sample size and key summary measures (mean, median, 25th and 75th percentiles). This distribution corresponds to uncertainty association with various additional parameters, such as number of variants within a gene and minor allele frequencies. Empirical distributions for these two parameters are displayed. If a range of EV is specified, plots and table for average power or sample size over the range of specified EV is returned. Fast option runs genome-wide calculations within 3 minutes and provides rough estimates and . Intermediate option runs genome-wide calculations within 6 minutes and provides more accurate estimates over 50 possible effect size distributions. Lastly, most accurate option runs genome-wide calculations within 15 minutes and provides very accurate estimates over 100 possible effect size distributions.

## Genome-level power calculation

### Essential Input parameter

1) M: Hypothesized number of underlying causal loci (or variants if analysis to be done based on single variant test)

2) GEV: Total EV explained by  loci in genome -wide study (see definition of EV above)

### Optional Input Parameter

1) m: The number of causal loci for which probability of discovery to be calculated (see output)
2) Level of complexity: The number of models and iterations used to estimate range of expected number of discoveries and probabilities (see output). There are three options: fast, intermediate and most accurate. Fast option runs genome-wide calculations within 3 minutes and provides rough estimates. Intermediate option runs genome-wide calculations within 6 minutes and provides more accurate estimates. Lastly, most accurate option runs genome-wide calculations within 15 minutes and provides very accurate estimates.



### Output

Expected number of discoveries: The application returns expected number of discoveries where the expectation is calculated across the M loci accounting for uncertainty associated with distribution of number of variants per locus (J), allele frequencies and the distributions of EVs the loci explains. Currently, the distribution of J and MAF in these calculations are obtained from those observed in the ExAC database. In addition, it is assumed the effect size distribution follows a L-shaped gamma distribution with mean specified as .  The application calculates a range of expected number of discoveries based on the range of the dispersion parameter of the underlying gamma distribution for the effect size distribution and the corresponding maximum and minimum values are returned.  
Probability of discoveries: This returns maximum and minimum probability of a certain number of discoveries (m) for values of specified by the user. 
           
 ## Additional notes
 
 Currently power calculations are based on distribution of number of variants per gene and minor allele frequencies observed in the public data for 60,000 individuals from Exome Aggregation Consortium [Lek et al. 2016]. It will underestimate the total number of variants per gene/region for whole-genome study and may overestimate power for gene-level tests. 

Contributor
===========
PAGEANT is designed by Andriy Derkach, Haoyu Zhang and Nilanjan Chatterjee

