# pgsc_calc_wdl

## pgsc_calc

WDL wrapper for the [pgsc_calc](https://pgsc-calc.readthedocs.io/en/latest/) workflow

The first step of the workflow [formats the input genomes](https://pgsc-calc.readthedocs.io/en/latest/how-to/prepare.html) in PLINK2 format.

The next step [runs the pgsc_calc nextflow workflow](https://pgsc-calc.readthedocs.io/en/latest/getting-started.html) inside a docker container. 

input | description
--- | ---
vcf | Array of VCF files. If provided, will be converted to pgen/pvar/psam. If not provided, use pgen/pvar/psam inputs instead.
pgen | Array of pgen files
pvar | Array of pvar files
psam | Array of psam files
chromosome | Array of chromosome strings (1-22, X, Y) corresponding to `vcf` or `pgen/pvar/psam`. If there is one file with multiple chromosomes, this input should be an empty string (`[""]`)
target_build | `"GRCh38"` (default) or `"GRCh37"`
pgs_id | PGS catalog IDs to calculate (e.g. `"PGS001229, PGS000802"`)
scorefile | Score file in the [required format](https://pgsc-calc.readthedocs.io/en/latest/how-to/calculate_custom.html).
ancestry_ref_panel | Google bucket path of a reference panel file (e.g. `"gs://fc-a8511200-791a-4375-bccf-fbe41ac3f9f6/pgsc_HGDP+1kGP_v1.tar.zst"`) to [perform ancestry adjustment](https://pgsc-calc.readthedocs.io/en/latest/explanation/geneticancestry.html). If not provided, no ancestry adjustment is performed.
sampleset_name | Name of the sampleset; used to construct output file names (default `"cohort"`)
arguments | [Additional arguments](https://pgsc-calc.readthedocs.io/en/latest/reference/params.html#param-ref) to pass to psgc_calc

Output files from pgsc_calc are described [here](https://pgsc-calc.readthedocs.io/en/latest/explanation/output.html#interpret).


## pgsc_calc_prepare_genomes

Standalone workflow to convert VCF to pgen/pvar/psam. The pvar generated from the workflow will have variant ids in the form chr:pos:ref:alt without the "chr" prefix.

input | description
--- | ---
vcf | Array of VCF files
merge_chroms | Boolean for whether to merge files to a single set of output files with all chromosomes
snps_only | Boolean for whether to keep only SNPs in output

output description
--- | ---
pgen | Array of pgen files
pvar | Array of pvar files
psam | Array of psam files


## calc_scores

Calculate scores without using Nextflow. Use pgsc_calc_prepare_genomes first to generate files.

The first two columns of the score file are expected to be "ID" (variant id) and "effect_allele". All subsequent columns are score weights, with the header of each column containing the score name.

If ancestry_adjust is true, a file with PCs must be supplied. The PC file is expected to have column "IID" for the sample ID and columns starting with "PC" for the PCs. For each score in the scorefile, the code fits a model regressing the score on the PCs. The coefficients from this model are used to adjust the score. Both unadjusted and adjusted scores are returned by the workflow.

If your pvar file contains variant ids starting with a "chr" prefix but your score file does not, set add_chr_prefix to "true" to add the "chr" prefix to the score file to match the pvar. pgsc_calc_prepare_genomes generates a pvar file without the "chr" prefix, so if you prepared your genomes with that workflow, "add_chr_prefix" should be set to "false" (default).

input | description
--- | ---
scorefile | score file
pgen | pgen file
pvar | pvar file
psam | psam file
harmonize_scorefile | Boolean for whether to harmonize scorefile to consistent effect allele (default true)
add_chr_prefix | Boolean for whether to add "chr" prefix to scorefile variant ids to match pvar (default false)
ancestry_adjust | Boolean for whether to adjust scores for ancestry using PCs
pcs | optional file with PCs to adjust for ancestry
subset_variants | optional list of variants to subset from the scorefile before running

output description
--- | ---
scores | sscore file output by plink
adjusted_scores | if ancestry_adjust is true, a file with adjusted scores is returned
variants | variants included in sscore
overlap | TSV file with fraction of overlapping variants and ratio of squared weights for each score
