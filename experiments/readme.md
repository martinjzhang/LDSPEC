## Data curation (Table 1): `job.curate_data`
Curate data for genotype, diseases/traits, SNP annotations. 
See details in [job.curate_data#readme](./job.curate_data#readme).

## Simulation (Fig. 1): `job.simulation`
- Generate simulation parameters: `get_config.ipynb`.
- Generate simulation data: `s1_simulate_paper.sh`.
- Run LDSPEC on the simulation data: `s2_regress_paper.sh`, `s2_regress_paper.vary_winsize.sh`.
- Make figures: `paper_simulations.main.ipynb`, `paper_simulations.other_model.ipynb`, `paper_simulations.vary_winsize.ipynb`.

## UKB data analyses (Figs. 2-5): `job.analysis_imp_geno_chimp`
LDSPEC pipeline:
- Compute LD: `gdreg_s1_get_ld.all_1e7.sh`, `gdreg_s1_get_ld.all_1e7.sh`.
- Compute sumstats: `ldspec_s2_get_sumstats.sh`.
- Compute average R for SNP-pair anntoations: `ldspec_s3_get_avgr.sh`.
- Compute LD scores for single-SNP annotations: `ldspec_s3_get_score_pannot.sh`, `ldspec_s3_get_score_pannot.vary_winsize.sh`.
- Run regression: `ldspec_s4_regress.sh`.

Compute score correlations: `score_corr.ipynb`

Figures
- Figures 2,3,5: `paper_main.ipynb`
- Figure 4: `paper_main.h2enrich_vs_cor.ipynb`
- Supp. Figures: `paper_main.comparison.ipynb`, `paper_main.individual_trait.ipynb`.

## Evolutionary simulations (Fig. 6): `job.evosimu`

## Example LDSPEC analysis: `job.example`
