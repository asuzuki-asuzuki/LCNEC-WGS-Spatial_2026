# Analysis of public data from the previous study

## Requirements
R (>=4.0; tested with 4.3.2)  

The data from the LCNEC study (Supplementary Data 1, 6 and 11, George J et al. 2018 <I>Nat Commun</I>) were used for the analysis.

## Analysis
### Mutation status of ADAMTS family genes  
- [mut_ADAMTS.pl](./mut_ADAMTS.pl): generating the table with mutation status of ADAMTS family genes and survival information.
- [survival_mut.R](./survival_mut.R): Kaplan-Meier analysis, comparing between mutation (+) and other cases

### Expression patterns of Module 17 (DP-specific module)  
- [exp_module17.pl](./exp_module17.pl): generating the table with expression levels of Module 17.
- [heatmap_Module17.R](./heatmap_Module17.R): conducting clsutering analysis and generating the heatmap.
- [exp_module17_surv.pl](./exp_module17_surv.pl): generating the table with clusters (expression patterns of Module 17) and survival information.
- [survival_exp.R](./survival_exp.R): Kaplan-Meier analysis, comparing among clusters of Module 17 expression. 
