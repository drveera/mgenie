 ## Important thing
 Make sure you've installed the mgenie first. (see home page readme of this repository)
 
 ## Check if you can call the script
when you type 

`mgenie predixcan qc-samples -h`

you should see a help message like this

```
usage:
 predixcan qc-samples [options] --bfile=FILE --ref=FILE

options:
 --bfile=FILE      plink file basename or list
 --ref=FILE        referrence plink file (only basename)
 --out=PREFIX      outname prefix [default: predixcan]
 --skip-sex        skip sex checks
 --cluster=NAME    cluster name  [default: open]
 --nojob           run in front end
 --dry-run         just show the codes
 --int             control job submission from frontend
 --njobs=NUMBER    number of parallel jobs; applicable only when running
                    in front end
 ```
 
 ## How to run the pipeline
To run this pipeline, you need two files mainly, 

1. plink file of your samples (.bed, .bim, .fam) 
2. plink file of your referrence genotypes, example, 1000g

you have to pass only the basename (without extension) to the pipeline. 

For example, i have the following files, lets say,

```
cmc.bed
cmc.bim
cmc.fam
1kg.bed
1kg.bim
1kg.fam
```
then, i'll run the pipeline using the following command,

`genie predixcan qc-samples --bfile cmc --ref 1kg --out cmc --cluster minerva`

## Referrence 1000 genome plink file in minerva cluster
if you are working in minerva, and if you have access to project `va-biobank`, then you are free to use the referrence file I generated for this pipeline. It's located in `/project/va-biobank/Veera/genie/resources/qc-dna/1000g/*bed`

## if you don't have chromosome X and sex information in the fam file 
First step in the pipeline is to check if the sex in the fam file and the sex determined using the genotypes are concordant. The discordant samples will be removed during the process. 
if you don't have chromosome X in the plink file or if you don't have sex information in the fam file, then use the argument `--skip-sex`
otherwise your jobs will fail. 

## what the pipeline does?

1. Checks if the sex info and genotype sex are same, if not, it removes the sample. 
2. it runs a heterozygosity test and removes the samples with F>0.2 
3. Then, it thins the plink file to 20,000 markers
4. it runs IBD analysis, identifies related pairs, i.e. pairs with piHAT > 0.2 and removes one individual in the pair randomly. 
5. merges the analysis samples with the referrence samples. it first identifies only the common markers between two bim files, subset both the analysis samples and ref samples to only the common markers and merges them to single plink file
6. the merged file is thinned to a random 30,000 markers
7. runs PCA, plots the PC1 and PC2, then using the ref samples, it determines the center point, draws an ellipsoid with 8SD and removes the samples that lie outside this ellipsoid. 
8. repeats PCA only with the selected samples, so you can use these PCS as covariates in your analysis
9. prepares the final plink file with only selected samples
10. split the plink file in individual chromosomes, 1-22
11. convert each chromosome file to vcf file
12. bgzip and tabix index the vcf file
13. these vcf.gz files can be uploaded straight to michigan server for imputation. 
