PAGEANT
==========
PAGEANT: Power Analysis for Genetic Association Tests
==========
Introduction
============
The application allows rapid power analysis for a variety of genetic association tests by specification of a few key parameters [Derkach et al. 2017]. Power calculations can be done at the level of single variant for simple trend test, at the level of a genes/regions for various complex aggregated tests [Neale et al. 2011, Derkach et al. 2013, Wu et al. 2011, Madsen and Browning 2009] and at the level of the whole genome for the assessment of overall yield of a study. The calculations currently uses underlying distribution of gene size and minor allele frequencies of variants observed in the in the public data for 60,000 individuals from Exome Aggregation Consortium [Lek et al. 2016]

Power for association test at the level of a single variant or a single gene/region
-----------------------------------------------
# Essential Input Parameters


1)	EV= percent of variation in a trait explained by loci (gene); 2) alpha = level of the test; 3) Number of cases/Number of controls in case-control study; 4) Total Sample Size in continuous trait study. 

Optional Input Parameters
=========================
1)	 Number of Variants = number of causal variants in a locus. 2) Range of EV= range of percentages of variation in a trait explained by loci; 3) Proportion of Causal = proportion of causal variants in a locus; 

Description
=============
Average power for specific locus is calculated from specified number of variants in a locus and proportion of variation in a trait explained by a locus. Average power for genome-wide study are calculated for specified proportion of variation in a trait explained by a locus, while averaging over different gene/region sizes across the genome. Currently these power calculations are based on public data for 60,000 individuals from Exome Aggregation Consortium [Lek et al. 2016]. It may underestimate number of variants per gene in whole-genome study, as a result it may overestimate an average power. We recommend to specify several values of Proportion of Causal for sensitivity analysis. 
Power calculations estimate an average power under three relationships: 1) there is no relationship between MAF and % of variations explained by a variant; 2) there is no relationship between MAF and effect size (log-OR) and 3) effect size (log-OR) is proportional to log10(MAF). 

Contributor
===========
PAGEANT is designed by Andriy Derkach, Haoyu Zhang and Nilanjan Chatterjee

