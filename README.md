A public repository for methylation risk scores derived from UCLA's electronic health record system. Included are the weights for a collection of medications, lab panels, and Phecodes, for which there was a significant association between the observed value and the 10-fold cross-validated predicted value. Importantly, the weights shared here are trained on the entire methylation dataset of 1050 ethnically diverse individuals, whereas the results in the manuscript are reported on 827 individuals for whom we have a set of measured baseline features (e.g. age, sex, genetic PCs). This is because in the manuscript we wished to compare ththe utility of adding genomics to the set of baseline features alone (thus requiring nested baseline features). Since previous works have shown that fitting an MRS using penalized regression with fixed covariates often performs worse than fitting just the penalized CpGs alone, we share here MRS built on only the CpGs[1]. 

We include mapping files for the Phecodes and medications in order to facilitate proper comparisons. The Phecodes were generated on ICD10 codes, and the mapping file contains comma-separated information on the phecode (which is reported in the weights) as well as the trait it maps to. Medications were grouped into subclasses, and the comma-separated file contains the code (which is reported in the weights) and the pharmaceutical subclass to which the code maps.

If you have any requests of specific outcomes (such as training on a very specific medication or within a certain age group, sex, or ethnicity), please do not hesitate to send us an email at mjthompson at ucla dot edu, and we will happily send you summary statistics as well as the weight of interest.

[1] Trejo Banos, D., McCartney, D.L., Patxot, M. et al. Bayesian reassessment of the epigenetic architecture of complex traits. Nat Commun 11, 2865 (2020). https://doi.org/10.1038/s41467-020-16520-1
