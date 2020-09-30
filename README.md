## MR-JTI

For questions:  
Dan Zhou <zdangm@gmail.com>  
Eric R. Gamazon <ericgamazon@gmail.com>  

## PREREQUISITES
To run JTI only, the R packages 'glmnet' and 'optparse' will be needed.  
To run MR-JTI, the R packages 'glmnet', 'optparse', and 'HDCI' will be needed.  

## JTI TUTORIAL

+ Training (You can also skip this step and directly download the pre-trained models.)  

Input file format:  
1. Genotype files  
JTI uses genotype files in plink format (.bed, .bim, .fam).  

2. Expression file  
The expression file contains five elements as listed below. Headers are required.  
tissue: tissue name  (required)  
sampleid: sampleid (SID) matched with the SID of genotype file.=, i.e. the second column of .fam file.  (required)  
exp: the measured expression level. The residual of normalized expression after adjusting for covariates is recommended here.  (required)  
exp_w: the expression similarity between the target tissue and the tissue in the second column.  (required)  
dhs_w: the DNase I hypersensitivity site (DHS) profile similarity between the target tissue and the tissue in the second column.  (required)  

3. Gene annotation  
JTI uses GENCODE v32.GRCh37 as the annotation for locating gene position  
A simplied annotation file is provied with the example files.  

An example of JTI model training  

```
script_dir='/replace_with_script_dir'  
example_file_dir='/replace_with_the_example_file_dir'  

Rscript ${script_dir}/JTI.r \
--tissue=Adipose_Subcutaneous \
--geneid=ENSG00000182957 \
--genotype_path=${example_file_dir}/jti_example_geno \
--expression_path=${example_file_dir}/jti_example_exp.txt \
--tmp_folder=${example_file_dir}/tmp \
--gencode_path ${example_file_dir}/gencode.v32.GRCh37.txt \
--out_path ${example_file_dir}
```

Options for JTI model training
```
--tissue  The target tissue name. (required)  
--geneid  ENSG geneid. Please provide the real ENSG geneid which will be used to find the chromosome and postion for the gene. (required)  
--plink_path  Path to plink software. If not specified, you will need to add the path of plink to $PATH in ~/.bash_profile  
--genotype_path Genotype file in plink bfile format '.bed, .bim, .fam' (required)  
--tmp_folder  The folder for temporary files. Will be cleaned up after model training. (required)  
--gencode_path  The gene annotation file. A GENCODE annotation file is provided with the example files. (required)  
--out_path  The output directory. The output file will be names as geneid_targettissue.txt (required)   
```

Results format  
The output file contains the following columns.  
gene: geneid  
rsid: snpid  
chr_bp: chromosome_position  
ref_allele: reference allele (uncounted allele when generating the dosage file.)  
counted_allele: counted_allele (counted allele when generating the dosage file.)  
weight: weight for each counted allele  
r2: cross-validation r2. The square of the correlation between the predicted and observed expression levels.  
p:  cross-validation p-value. The significance of the correlation test (correlation between the predicted and observed expression levels)  
lambda: The final hyperparameter.  

+ Pretrained prediction models (GTEx v8) are available. For summary-stat-based analysis, SNP-SNP covariance matrices are also provided.


+ GWAS summary statistic based TWAS  
For convenience, we redirect you to the S-PrediXcan (https://github.com/hakyimlab/MetaXcan/blob/master/software/SPrediXcan.py), which has been developed to perform PrediXcan using summary statistics (Barbeira, Alvaro N., et al. "Exploring the phenotypic consequences of tissue specific gene expression variation inferred from GWAS summary statistics." Nature communications 9.1 (2018): 1-20).  

# MR-JTI TUTORIAL

Input file format  
The input file contains six elements, as listed below. (The headers are required.) 

rsid: rsid. SNPs need to be pruned before runing MR-JTI. 
effect_allele:  The effect allele. Harmonization needs to be performed to make sure the effect alleles of eQTL and GWAS are correctly aligned.  
ldscore: The LD score of each SNP. GCTA was used to generate LD score. gcta64 --bfile test --ld-score --ld-wind 1000 --ld-rsq-cutoff 0.01 --out test  (required)  
eqtl_beta:  Effect size of eQTL  (required)  
eqtl_se: SE of eQTL effect size  (required)  
eqtl_p: eQTL p-value  (required)  
gwas_beta: GWAS effect size  (required)  
gwas_p: GWAS p-value  (required)  

A typical run
```
script_dir='/replace_with_script_dir'  
example_file_dir='/replace_with_the_example_file_dir'  

Rscript ${script_dir}/MR-JTI.r \
--df_path ${example_file_dir}/mrjti_example.txt \
--n_genes 100 \
--out_path ${example_file_dir}/mrjti_result.csv

```

Options for JTI model training
```
--df_path Path to dataframe of GWAS and eQTL summary statistics (required)  
--n_folds Number of cross-validation folds. Default=5.  
--n_snps  Minimum number of SNPs (obs) to run MR-JTI. Default=20.  
--n_bootstrap  Number of resampling times for bootstrap. Default=500. #we used B=100 for our specific dataset, the default (recommended by Efron) is more general.
--n_genes Total number of genes tested (Bonferroni correction will be applied). n_genes=1 denotes user requires only nominal significance level (i.e., p<0.05 will be considered as significant). 
--weighted Weighted by inverse variance. Default=F
--out_path Output path.  (required)  
```

Results format
MR-JTI generates the upper and lower estimates of the gene's effect on GWAS trait as well as the heterogeneity estimates.  
variable: Variables including the gene's effect and the heterogeneity effects  
beta: Point estimate of the effect size  
beta_CI_lower:  Bonferroni adjusted confidence interval (CI), lower  
beta_CI_upper:  Bonferroni adjusted CI, upper  
CI_significance:  Significant if the CI does not overlap the null hypothesis (i.e. 0).  


