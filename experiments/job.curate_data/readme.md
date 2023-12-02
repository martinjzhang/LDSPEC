Curate annotation data

Meta information
- curate_info.ipynb

Genotype and phenotype
- Get samples: get_sample_list.ipynb
- Get phenotype: get_sample_list.ipynb
- Get genotype: get_plink.sh

- Check files: baseline_check_results.ipynb

Single-SNP annotations
- From .bed files: baseline_s1_convert_bed.py; baseline_s1_convert_bed.sh; obtain most baseline-LD annotations from .bed files 
- Nucleotide diversity: baseline_s2_nucleotide_diversity.R; baseline_s2_nucleotide_diversity.sh; provide your own .bim file to baseline_s2_nucleotide_diversity.sh to get the corresponding nucleotide diversity information. baseline_s2_nucleotide_diversity.sh calls baseline_s2_nucleotide_diversity.R and computes the diversity based on UK10K data. It does not take a long time (<1hr for all 22 chromosomes).
- Recobination rate: baseline_s3_recomb_rate.sh; provide your own .bim file to baseline_s3_recomb_rate.sh. Please note that recombination_rate.sh contains 3 steps, as commented in the script. Each step should be finished in <1min for each chromosome. 
- CADD annotations: baseline_s4_CADD.download.sh; baseline_s4_CADD_annot.py; baseline_s4_CADD_annot.sh; provide your own .bim file to baseline_s4_CADD_annot.sh. Please note that baseline_s4_CADD_annot.sh contains 3 steps, as commented in the script. Step1 should be run in python and takes <1min. Step2 should be run in bash. It uses tabix to extract information and takes ~2hrs for each bim file. Step3 calls a python script from bash and takes 4hrs for all 22 chromosomes. 
- SnpEff annotations: baseline_s5_snpeff.sh; baseline_s5_snpeff.ipynb
- VEP annotations: baseline_s6_VEP.sh; baseline_s6_vep.ipynb
- LLD_AFR: baseline_s7_LLD_AFR.ipynb
- Combine baseline-LF single-SNP annotations: baseline_combine.ipynb

SNP-pair annotations
- Gene reference: pannot_reference.ipynb
- Creat SNP-pair annotations: get_pannot_basic_call.sh; get_pannot_call.sh; get_pannot_fct_call.sh; get_pannot_gene_call.sh
- Check annotations: pannot_check_basic.ipynb; pannot_check_gene.ipynb; pannot_check_fct.ipynb
- Distance distritution: pannot_stats.ipynb
