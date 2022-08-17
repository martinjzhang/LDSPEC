s1 : obtain most baseline-LD annotations from .bed files 

s2 : nucleotide_diversity
    - supplement your own .bim file to nucleotide_diversity.sh to get the corresponding nucleotide diversity information. nucleotide_diversity.sh calls nucleotide_diversity.R and computes the diversity based on UK10K data. It does not take a long time (<1hr for all 22 chromosomes).

s3 : recombination rate
    -supplement your own .bim file to recombination_rate.sh. Please note that recombination_rate.sh contains 3 steps, as commented in the script. Each step should be finished in <1min for each chromosome. 
    
s4 : 

s5 : 

- CADD annotations: supplement your own .bim file to get_CADD_annot.sh. Please note that get_CADD_annot.sh contains 3 steps, as commented in the script. Step1 should be run in python and takes <1min. Step2 should be run in bash. It uses tabix to extract information and takes ~2hrs for each bim file. Step3 calls a python script from bash and takes 4hrs for all 22 chromosomes. 