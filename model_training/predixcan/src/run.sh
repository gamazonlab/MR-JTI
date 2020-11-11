#!/bin/bash

#SLURM_ARRAY_TASK_ID=1  #1 to N

#main_dir=${your_dir}/predixcan

#step 1. model training 
#ml GCC OpenMPI R
#Rscript ${main_dir}/src/predixcan.r \
#        --model_training \
#        --main_dir ${main_dir} \
#        --plink_file_name ${main_dir}/chr22/chr22 \
#        --expression_file_name ${main_dir}/chr22/exp_chr22.txt \
#        --subjob_id ${SLURM_ARRAY_TASK_ID} \
#        --n_genes_for_each_subjob 100 \
#        --annotation_file_name ${main_dir}/chr22/gencode.v32.GRCh37.txt


#step 2. generate .db and .cov file
#ml GCC OpenMPI R

#model_name=chr22 #the predix of output files
#Rscript ${main_dir}/src/predixcan.r \
#         --generate_db_and_cov \
#         --main_dir ${main_dir} \
#         --plink_file_name ${main_dir}/chr22/chr22 \
#         --expression_file_name ${main_dir}/chr22/exp_chr22.txt \
#         --annotation_file_name ${main_dir}/chr22/gencode.v32.GRCh37.txt \
#         --output_file_name ${model_name}


#step 3. apply to gwas data
#ml GCC/5.4.0-2.26  OpenMPI/1.10.3 pandas/0.18.1-Python-2.7.12 numpy/1.11.1-Python-2.7.12 scipy/0.17.0-Python-2.7.12

#model_name=chr22

#python2.7 ${your_dir_for_metaxcan_or_Spredixcan}/MetaXcan.py \
#--model_db_path ${main_dir}/output/${model_name}.db \
#--covariance ${main_dir}/output/${model_name}.cov \
#--gwas_file ${main_dir}/chr22/LDL_with_Effect_chr22.txt \
#--snp_column MarkerName \
#--effect_allele_column Allele1 \
#--non_effect_allele_column Allele2 \
#--beta_column Effect \
#--pvalue_column GC.Pvalue \
#--output_file ${main_dir}/output/chr22.asso.csv



