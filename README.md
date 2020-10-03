## MR-JTI

For questions:  
Dan Zhou <zdangm@gmail.com>  
Eric R. Gamazon <ericgamazon@gmail.com>  

## PREREQUISITES
To run JTI only, the R packages 'glmnet' and 'optparse' will be needed.  
To run MR-JTI, the R packages 'glmnet', 'optparse', and 'HDCI' will be needed.  

## JTI TUTORIAL
### Training 
(You can also skip this step and directly download the pre-trained models.)  

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


## Pre-trained prediction models (in GTEx v8) are available. 
For summary-statistics-based analysis, SNP-SNP covariance matrices are also provided. Three models are provided here.

PrediXcan: 
Using elastic net, five-fold cross-validation was performed to determine the optimal hyper-parameter (lambda) with minimal cross-validation error. The prediction performance was evaluated using the Pearson’s correlation between the predicted expression (the same lambda but from five different models trained from five folds) and the observed expression. Here, we defined an ‘imputable’ gene as one with Pearson’s correlation r > 0.1 and P < 0.05 (see Online Methods and Supplementary Figure 3 of the paper). This choice of threshold considers not just the significance but a reasonable magnitude of the proportion of variance explained (PVE) by genetic variants. We filtered out genes with a negative correlation between the predicted and observed expression level (in contrast to the PrediXcan definition based on r2).  
Reference: Gamazon, Eric R., et al. "A gene-based association method for mapping traits using reference transcriptome data." Nature Genetics 47.9 (2015): 1091.(https://www.nature.com/articles/ng.3367)

UTMOST:
Using a sparse group-LASSO, Cross Tissue gene expression IMPutation (CTIMP, the model training part of the UTMOST framework) borrows information across tissues and significantly increases the prediction accuracy compared to the conventional PrediXcan. In the original version of UTMOST, fold-specific hyper-parameters were used. An imputable gene, which is used for downstream analysis, was defined by PFDR < 0.05 from the correlation test between the predicted and observed expression. Importantly, the original UTMOST used the re-trained model (in the entire dataset) to generate the predicted expression. This implies that the estimation of both the p-value and the correlation r may be inflated. The FDR correction does not fix the problem of inflated correlation.
In order to facilitate model evaluation and comparison of the different approaches, we proposed a modification of the model training of UTMOST  . Briefly, we used uniform hyper-parameters across different folds to make the hyper-parameters directly comparable. This modification provided an approximately unbiased estimate of prediction performance (as we confirmed in external test datasets), facilitating comparison with PrediXcan. Actually, UTMOST did provide a valid way for model evaluation (using the test sets). However, since the lambdas are not consistent across the different folds, the predicted expression in the test sets cannot be pooled to give an overall estimate of performance. As a result, the performance has to be estimated in each of the test sets which has only 20% (given the 5 folds) of the samples. Although the prediction performance is not inflated using this approach, the point estimate in the test set will not be very precise because of the substantially reduced sample size. Anyway, the original UTMOST did not use the evaluation in the test set to define a ‘heritable gene’ (i.e., an imputable gene) for downstream analysis. Instead, the original UTMOST used the re-trained model to define a heritable gene.  
Here we provide the pre-trained prediction models for imputable genes (same definition and consistent evaluation as the conventional PrediXcan) using our modification to UTMOST. In order to maximize the model accuracy, we also retained the model in the final step with all samples using optimized hyper-parameters. But the performance evaluation was performed in the cross-validation step. i.e., not in the retrained step.  
Reference: Hu, Yiming, et al. "A statistical framework for cross-tissue transcriptome-wide association analysis." Nature Genetics 51.3 (2019): 568-576.(https://www.nature.com/articles/s41588-019-0345-7)

JTI:
Joint-tissue imputation (JTI) borrows information across tissues to improve the prediction quality by leveraging expression and epigenetic (chromatin accessibility) similarity across different tissues / cell-types. The hyper-parameter tuning provides the flexibility to reduce to the conventional PrediXcan when the expression and regulatory profile in the target tissue is unique (with little to borrow from the other tissues). In contrast to UTMOST, which treats all tissues equivalently and does not leverage the tissue similarity, JTI seeks to exploit this similarity for improved prediction. The pre-trained JTI models are provided for imputable genes (same definition and consistent evaluation as the conventional PrediXcan).  
Reference: Zhou, Dan, et al. "A unified framework for joint-tissue transcriptome-wide association and Mendelian randomization analysis." Nature Genetics (2020). (https://www.nature.com/articles/s41588-020-0706-2)

### GWAS summary statistic based TWAS  
For convenience, we redirect you to the S-PrediXcan (https://github.com/hakyimlab/MetaXcan/blob/master/software/SPrediXcan.py), which has been developed to perform PrediXcan using summary statistics (Barbeira, Alvaro N., et al. "Exploring the phenotypic consequences of tissue specific gene expression variation inferred from GWAS summary statistics." Nature communications 9.1 (2018): 1-20).  

## MR-JTI TUTORIAL

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


## Cite MR-JTI
Zhou, Dan, et al. "A unified framework for joint-tissue transcriptome-wide association and Mendelian randomization analysis." Nature Genetics (2020). (https://www.nature.com/articles/s41588-020-0706-2)
