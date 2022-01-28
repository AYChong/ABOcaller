# ABOcaller
ABOcaller: a python script to easily call ABO blood type from genotype data


## Requirements
- Python3
- pandas
- numpy


## Background
### ABO types
The ABO gene encodes the glycosyltransferase enzymes responsible for the presentation of A or B antigens on cell surfaces. The most common ABO alleles can be defined using four variable positions, and genetic inference of ABO blood type generally relies on a combination of these (A/B alleles: rs8176743, rs8176746, and rs8176747; O allele: rs8176719). ABOcaller uses rs8176747 (AB) and rs8176719 (O) to infer ABO blood type, and at minimum, only requires that these variants are included in the input genotype file. 

Imputation servers may return phased VCF files with SNP IDs in the format "chr:position:ref:alt". These can be used directly by including the "--rs8176719" and "--rs8176747" flags and the appropriate SNP IDs. 

If the data is not already phased, we recommend phasing the region spanning a minimum of 1Mb across the ABO gene region using IMPUTE2 or SHAPEIT and using the resulting haplotype file as input alongside the appropriate sample file. Unphased data can be used by including the flag “--unphased”. 

### FUT2 secretor status
The FUT2 gene encodes the Fucosyltransferase 2 enzyme responsible for the expression of ABO antigens on mucosal endothelial cells. For the determination of FUT2 secretor status, a VCF file or .gen file containing genotypes for rs601338 can be provided as input. As only a single variant is used to determine status, data does not need to be phased.

### Variant positions
| rsID  | GRCh37 | GRCh38 |
| ------------- | ------------- | ------------- |
| rs8176747  | 9:136131315  | 9:133255928  |
| rs8176719  | 9:136132909  | 9:133257521  |
| rs601338   | 19:49206674  | 19:48703417  |

## Usage

```
 usage: ABOcaller.py [-h] -i INFILE [-s SAMPLE] -o OUTFILE -m {ABO,Se}                                                                                            [--rs8176719 O_VARIANT] [--rs8176747 AB_VARIANT]
                     [--rs601338 SE_VARIANT] [--gp_cutoff GP_THRESHOLD]
                     [--unphased]

required arguments:
  -i INFILE, --genotype_file INFILE
                        File containing genotypes for either ABO or FUT2
                        variants
  -o OUTFILE, --output_file OUTFILE
                        File name for the resulting blood type calls
  -m {ABO,Se}, --mode {ABO,Se}
                        Type of blood group to call. Valid options are ABO and
                        Se

optional arguments:
  -s SAMPLE, --sample_file SAMPLE
                        File containing sample IDs. Required if providing .hap
                        or .gen input
  --rs8176719 O_VARIANT
                        Alternate variant ID for ABO variant rs8176719.
                        Default is rs8176719.
  --rs8176747 AB_VARIANT
                        Alternate variant ID for ABO variant rs8176747.
                        Default is rs8176747.
  --rs601338 SE_VARIANT
                        Alternate variant ID for FUT2 variant rs601338.
                        Default is rs601338.
  --gp_cutoff GP_THRESHOLD
                        Posterior probability threshold for calling genotypes
                        from imputed data. Used when parsing .vcf files and
                        .gen files. Default is 0 (no threshold).
  --unphased            Use to indicate that input data in vcf format is
                        unphased. Default behaviour is to assume vcf data is
                        phased

```
### Output files 

Using ABOcaller in either ABO mode or Se mode will return a tab separated output file containing either ABO  or Se genotypes and blood type, one line per individual, and columns as outlined below: 

- SID: Sample ID as provided in the VCF header, or concatenated as FID_IID from gen/haps input 

- ABO_haplotype: ABO blood type as a two letter haplotype, with inferred ABO status for each strand 

- ABO_bloodtype: ABO blood type in O/A/B/AB format 

- O_genotype: Genotype at rs8176719 or the given O allele variant 

- AB_genotype: Genotype at rs8176747 or the given A/B allele variant 

- Se_genotype: Genotype at rs601338 or the given secretor status variant 

- Se_status: Secretor status with secretors coded as 0 and non-secretors coded as 1 

### Examples
ABO status from vcf file
```
python ABOcaller.py \
  -i ./example_data/ABO.vcf \
  -o ABO.tsv \
  -m ABO
```

ABO status from vcf file with alternative SNP IDs
```
python ABOcaller.py \
  -i ./example_data/ABO.vcf \
  -o ABO.tsv \
  -m ABO
  --rs8176719 chr9:133257521:T:TC
  --rs8176747 chr9:133255928:C:G
```

ABO status from unphased data
```
python ABOcaller.py \
  -i ./example_data/ABO.haps \
  -s ./example_data/ABO.sample \
  -o ABO.tsv \
  -m ABO \
  --unphased
```

Se status from vcf file
```
python ABOcaller.py \
  -i ./example_data/Se.vcf \
  -o Se.tsv \
  -m Se
```

Se status from .gen file using cutoff of 0.8
```
python ABOcaller.py \
  -i ./example_data/Se.gen \
  -s ./example_data/Se.sample \
  -o Se.gen.tsv \
  -m Se \
  --gp_cutoff 0.8
```

