# Topy example of LDSPEC pipeline using LDSPEC CLI
Toy pipeline showcasing steps for running LDSPEC.
- Example code in `ldspec_pipeline.sh`, example [input files](https://github.com/martinjzhang/LDSPEC/tree/main/ldspec/data), and example [results](https://github.com/martinjzhang/LDSPEC/tree/main/ldspec/data/ldspec_res).
- See [CLI_ldspec.py](https://github.com/martinjzhang/LDSPEC/blob/main/CLI_ldspec.py) for more details.

### File formats
- `--job`: one of `get_snp_block`, `compute_ld`, `compute_score`, `compute_avgr`, `regress`.
- `--pgen_file`: prefix of PLINK2 pgen files. `@` refers to CHR number.
- `--ld_file`: LDSPEC LD files. 
- `--annot_file`: LDSPEC annotation files.
- `--score_file`: LDSPEC LD score files.
- `--snp_range_file`: LDSPEC SNP range file. 
- `--sumstats_file`: Summary statistics files, same as in LDSC.
- `--avgr_file`: LDSPEC `.avgr` file for average r in SNP-pair annotations. 
- `--prefix_out`: Output prefix.
- `--snp_range`: SNP range.
- `--win_size`: window size.
- `--flag_nofil_snp`: If to filter SNPs in regression.

## Step 1: get SNP blocks
Create a `.txt` file with SNP ranges. Output example: [ldspec_s1.snp_range.txt](https://github.com/martinjzhang/LDSPEC/blob/main/ldspec/data/ldspec_res/ldspec_s1.snp_range.txt).

    job=get_snp_block: generate SNP blocks (10,000 SNPs per block)
    - Input : --job | --pgen_file | --prefix_out
    - Output : list of snp ranges (in the format of "snp_range"), e.g., `c1_s0_e10000`.


## Step 2: compute LD matrices
Compute LD matrix for each SNP block in the `snp_range.txt` file; can parallelize.

    job=compute_ld : compute LD between target SNPs in `snp_range` and reference SNPs within `win_size` of target SNPs
    - Input : --job | --pgen_file | --prefix_out | --snp_range | [--win_size]
    - Output : `_ld.npz` file; LD matrix.

## Step 3: compute LD score
Compute LD score for each SNP block in the snp_range.txt file; can parallelize.

`--annot_file` takes comma seperated `.annot.gz` and `.pannot_mat.npz` files, or a `.txt` file of file paths.

    compute_score : compute LD and DLD scores.
    - Input : --job | --pgen_file | --ld_file | --annot_file | --prefix_out | [--win_size] | [--flag_cross_term]
    - Output : `_score.tsv.gz` file; LD and DLD scores.

## Step 4: combine LD scores
`--score_file`: `@` represents `snp_range`

    combine_score : concatenate scores from the same CHR 
    - Input : --score_file | --snp_range_file | --prefix_out
    - Output : concatenated score files by CHR

## Step 5: compute AVGR
Compute average LD (avgr) for each pannot.

    compute_avgr : compute average LD (avgr) for each pannot. "--ld_file" should contain all LD files
    - Input : --job | --pgen_file | --ld_file | --annot_file | --prefix_out
    - Output : Average LD for each pannot.

## Step 6: regress
Estimate LDSPEC results. `--avgr_file` supports comma-separated `.avgr` files.

    regress : estimate LDSPEC parameters.
    - Input : --job | --pgen_file | --annot_file  | --score_file | --sumstats_file| --avgr_file | --prefix_out | [--flag_cross_term] | [--flag_nofil_snp]
    - Output : LD-SPEC result.
