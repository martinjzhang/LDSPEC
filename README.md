# LDSPEC
LDSPEC (LD SNP-pair effect correlation regression) is a method for estimating the correlation of causal disease effect sizes for pairs of nearby SNPs, depending on their functional annotations. 

### Reference
[Zhang, et al. "Pervasive correlations between causal disease effects of proximal SNPs vary with functional annotations and implicate stabilizing selection"](XX), preprint, 2023.

### Versions
- 0.0.3: beta version accompanying initial submission. 

### Code and data for the paper
See [experiments](https://github.com/martinjzhang/LDSPEC/tree/main/experiments) for more details. Data are at [figshare](https://figshare.com/projects/LD_SNP-pair_effect_correlation_regression_LDSPEC_/188052). 
- Download [baseline-SP annotations](https://figshare.com/articles/dataset/LDSPEC_data_release_120223/24716877) (165 single-SNP annotations and 136 SNP-pair annotations).
- Download corresponding [LD scores and directional LD scores](https://figshare.com/articles/dataset/LDSPEC_data_release_120223_ldscore/24716901)
- Download [GWAS summary statistics](https://figshare.com/articles/dataset/LDSPEC_data_release_120223_sumstats/24716952) for 70 UK Biobank diseases/traits.

## Using LDSPEC
### Installation 

    git clone https://github.com/martinjzhang/LDSPEC.git 
    cd LDSPEC
    pip install -e .

### Usage
LDSPEC consists of many steps and may be tricky to run. We are working hard to make it more user-friendly. If you want to use it, please check out this [self-contained example](https://github.com/martinjzhang/LDSPEC/blob/main/experiments/job.example/readme.md) first for steps of LDSPEC and inputs and outputs in each step. Let us know if you have any questions or suggestions for the software.

## Contact
- Martin Zhang (martinzh@andrew.cmu.edu).  
