# ABOcaller
ABOcaller: a python script to easily call ABO blood type from genotype data


## Requirements
- Python3
- pandas
- numpy

## Usage
```
usage: ABO_Se.haps_rework.py [-h] -i INFILE [-s SAMPLE] -o OUTFILE -m {ABO,Se}
                             [--o_variant RS8176719] [--ab_variant RS8176747]
                             [--se_variant RS601338] [--gp_cutoff GP_THRESHOLD]
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
  --rs8176719 RS8176719
                        Alternate snp ID for ABO variant rs8176719. Default is
                        rs8176719.
  --rs8176747 RS8176747
                        Alternate snp ID for ABO variant rs8176747. Default is
                        rs8176747.
  --rs601338 RS601338   Alternate snp ID for FUT2 variant rs601338. Default is
                        rs601338.
  --gp_cutoff GP_THRESHOLD
                        Posterior probability threshold for calling genotypes
                        from imputed data. Used when parsing .vcf files and
                        .gen files. Default is 0.8.
  --unphased            Use to indicate that input data in vcf format is
                        unphased. Default behaviour is to assume vcf data is
                        phased
```


### Examples
ABO status from vcf file
```
python ABO_Se.py \
  -i ABO.vcf \
  -o ABO.tsv \
  -m ABO
```

Se status from vcf file
```
python ABO_Se.py \
  -i Se.vcf \
  -o Se.tsv \
  -m Se
```

ABO status from unphased data
```
python ABO_Se.py \
  -i ABO.haps \
  -s ABO.sample \
  -o ABO.tsv \
  -m ABO \
--unphased
```

Se status from .gen file using default cutoff of 0.8
```
python ABO_Se.py \
  -i Se.gen \
  -s Se.sample \
  -o Se.gen.tsv \
  -m Se \
  --gp_cutoff
```

