# Functions LDSPEC CLI
Example code in `ldspec_pipeline.sh`

## Step 1: get SNP blocks

    job=get_snp_block : generate SNP blocks (10,000 SNPs per block)
        - Input : --job | --pgen_file | --prefix_out
        - Output : list of snp ranges (in the format of "snp_range"), e.g., `c1_s0_e10000`.


## Step 2: compute LD matrices
Compute LD matrix for each SNP block in the snp_range.txt file; can parallelize.

    job=compute_ld : compute LD between target SNPs in `snp_range` and reference SNPs within `win_size` of target SNPs
        - Input : --job | --pgen_file | --prefix_out | --snp_range | [--win_size]
        - Output : `_ld.npz` file; LD matrix.

## Step 3: compute LD score
    compute_score : compute LD and DLD scores.
        - Input : --job | --pgen_file | --ld_file | --annot_file | --prefix_out | [--win_size] | [--flag_cross_term]
        - Output : `_score.tsv.gz` file; LD and DLD scores.

## Step 4: combine LD scores
    combine_score : concatenate scores from the same CHR 
        - Input : --score_file | --snp_range_file | --prefix_out
        - Output : concatenated score files by CHR

## Step 5: compute AVGR
    
    compute_avgr : compute average LD (avgr) for each pannot. "--ld_file" should contain all LD files
        - Input : --job | --pgen_file | --ld_file | --annot_file | --prefix_out
        - Output : Average LD for each pannot.

## Step 6: regress
    regress : estimate LDSPEC parameters.
        - Input : --job | --pgen_file | --annot_file  | --score_file | --sumstats_file| --avgr_file | --prefix_out |
        [--flag_cross_term] | [--flag_nofil_snp]
        - Output : LD-SPEC result.
