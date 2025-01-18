# Preprocessing summary statistics

After downloading the summary statistics, we need the following steps to pre-process it to conduct the analysis in the paper

- Make sure each GWAS file is a ".csv" or ".txt" file containing a data frame

- Make sure the GWAS file has at least 6 columns with these column names: SNP, effect_allele, other_allele, beta, se, pval. The SNP column contains rsID for each SNP. Both the effect_allele and other_allele columns need to have capital letters. The beta column contains the estimated effect size for continuous traits and log odds ratio for binary trait, and the se column is the standard deviation of the corresponding beta. If names are different, then one should rename the columns
