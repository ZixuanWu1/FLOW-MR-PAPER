In this work we use ldsc to compute genetic correlation.


For instance, we use ```munge``` to convert data into .sumstat format

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

Then one can use ldsc to compute genetic correlation

```
python ldsc.py \
--rg childhood_body_size.sumstats.gz, adult_bmi.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out scz_bip
```
