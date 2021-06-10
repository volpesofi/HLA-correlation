# HLA-correlation

From hla typing results (`input/fished_all_tops.tsv`) we produce a count matrix with patients on the row and HLA on the column.
The allele count can be: 
  * 2 (homozygosis for hla_i) 
  * 1 (heterozygosis for hla_i) 
  * 0 (hla_i not present)

From the metadata (`input/metadata_Vo.xlsx`), we subset the rows with information for WGS samples 
and we add to the count matrix as many columns as the metadata information with the phenotypes we wish to assess the correlation for.

The correlation between phenotype and genotype is then performed with two R functions:
  * `lm(var ~ hla)` - to correlate the genotype (hla) with continuos variables (var), _i.e._ "Abbot_semiquantitative", "Roche_Total_ICO", "Diasorin_IgG_semiquantitative" and the residuals.
  * `glm(var ~ hla , family = "binomial")` – to correlate the genotype (hla) with discrete variable (var), _i.e._ "Groundtruth_GTA", "Groundtruth_direct_contacts_GTB", "Groundtruth_indirect_contacts_GTC", "swabs"

The script outputs 
 * .tsv tables with the statistical information of the linear correlations (in a folder called `tables_fished_tops`). In particular, the files 
`results_hla_a_allmodels.tsv`, `results_hla_b_allmodels.tsv`, `results_hla_c_allmodels.tsv` contain the Estimate and canonical pvalue for all the corralations.
 * several plots (heapmaps of -Log10(canonical pvalue), boxplots and mosaic plots), in a folder `called plots_fished_tops`.
 * The .Rdata objects 
