In this work we use how to ldsc to compute genetic correlation.

First one should clone the ```ldsc``` repo and activate the environment

```
git clone https://github.com/bulik/ldsc.git
cd ldsc
conda env create --file environment.yml
source activate ldsc
```

Then one should downlaod the required allele files using the following command

```
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2

```

Then one can use ldsc to estimate genetic correlation from GWAS files. Make sure the GWAS txt files have the following columns 

1. SNP: A unique identifier (e.g., the rs number)
2. a1: Allele 1 (effect allele)
3. a2: Allele 2 (non-effect allele)
4. N: Sample size (which often varies from SNP to SNP)
5. P: A P-value
6. BETA: A signed summary statistic (beta, OR, log odds, Z-score, etc)


Then one can use ```munge``` to convert data into .sumstat format. For example, suppose we are interested in the genetic correlation between childhood body size and adult bmi, then one can run

```
python munge_sumstats.py \
--sumstats childhood_body_size.txt \
--out body_size \
--merge-alleles w_hm3.snplist

```


```
python munge_sumstats.py \
--sumstats adult_bmi.txt \
--out body_size \
--merge-alleles w_hm3.snplist

```

Finally we can use ldsc to compute genetic correlation.

```
python ldsc.py \
--rg childhood_body_size.sumstats.gz, adult_bmi.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out scz_bip
```
