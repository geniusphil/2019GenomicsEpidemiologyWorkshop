# COSMIC

### For annovar db use

## GRCh38/cosmic/v89/
### COSMIC Mutation Data
* [CosmicMutantExport.tsv.gz](https://cog.sanger.ac.uk/cosmic/GRCh38/cosmic/v89/CosmicMutantExport.tsv.gz)

### COSMIC Non coding variants
* [CosmicNCV.tsv](https://cog.sanger.ac.uk/cosmic/GRCh38/cosmic/v89/CosmicNCV.tsv)

### VCF Fils (coding and Non-conding mutations)
* [CosmicCodingMuts.vcf.gz](https://cog.sanger.ac.uk/cosmic/GRCh38/cosmic/v89/VCF/CosmicCodingMuts.vcf.gz) 
    
* [CosmicNonCodingVariants.vcf.gz](https://cog.sanger.ac.uk/cosmic/GRCh38/cosmic/v89/VCF/CosmicNonCodingVariants.vcf.gz)


```
$prepare_annovar_user_mod.pl -dbtype cosmic CosmicMutantExport.tsv -vcf CosmicCodingMuts.vcf > hg38_cosmic89_coding.txt
$prepare_annovar_user_mod.pl -dbtype cosmic CosmicNCV.tsv -vcf CosmicNonCodingVariants.vcf > hg38_cosmic89_noncoding.txt
```

---

## GRCh37/cosmic/v89/
### COSMIC mutation data
* [CosmicMutantExport.tsv.gz](https://cog.sanger.ac.uk/cosmic/GRCh37/cosmic/v89/CosmicMutantExport.tsv.gz)

### COSMIC Non coding variants
* [CosmicNCV.tsv](https://cog.sanger.ac.uk/cosmic/GRCh37/cosmic/v89/CosmicNCV.tsv)

### VCF Fils (coding and Non-conding mutations)
* [CosmicCodingMuts.vcf.gz](https://cog.sanger.ac.uk/cosmic/GRCh37/cosmic/v89/VCF/CosmicCodingMuts.vcf.gz) 

* [CosmicNonCodingVariants.vcf.gz](https://cog.sanger.ac.uk/cosmic/GRCh37/cosmic/v89/VCF/CosmicNonCodingVariants.vcf.gz)

```
prepare_annovar_user_mod.pl -dbtype cosmic CosmicMutantExport.tsv -vcf CosmicCodingMuts.vcf > hg19_cosmic89_coding.txt
prepare_annovar_user_mod.pl -dbtype cosmic CosmicNCV.tsv -vcf CosmicNonCodingVariants.vcf > hg19_cosmic89_noncoding.txt
```