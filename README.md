# VCF Parser #

[![Build Status](https://travis-ci.org/moonso/vcf_parser.svg)](https://travis-ci.org/moonso/vcf_parser)

Small library for parsing vcf files. Based on [PyVCF](https://github.com/jamescasbon/PyVCF)

usage examples:
    
    - Split variants lines with multiple alleles
    - Create vcf files from a python environment and print them
    - Parse through variants and retrieve relevant information

## Installation ##


    pip install vcf_parser

or 

```bash
>git clone https://github.com/moonso/vcf_parser.git
>cd vcf_parser
>python setup.py install
```

## Usage ##


If used within a python environment:

```python3
>from vcf_parser import VCFParser
>my_parser = VCFParser(infile='infile.vcf', split_variants=True)
>for variant in my_parser:
    print(variant)
```

or used as a command line tool

    vcf_parser examples/test_vcf.vcf --split

Prints a new vcf with splitted variants to screen.

Vcf parser is really a lightweight version of [PyVCF](https://github.com/jamescasbon/PyVCF) with most of it's code borrowed and modified from there.

The idea was to make a faster and more flexible tool that mostly work with python dictionaries.

It is easy to access information for each variant, edit the information and edit the headers.

## Basic function ##


Returns dictionary with the vcf info for each variant.
To split the multiallelic calls(and accurate splitting of INFO field including the VEP CSQ fiels) use:
    
    my_parser = VCFParser(infile='infile.vcf', split_variants=True)

The ordinary vcf entrys is stored by there header names, like
    
    variant['CHROM']
    variant['ALT']

etc.

The genotype information is converted to a genotype object and stored in a dictionary

    variant['genotypes']

and looks like:

```python
'genotypes': {'father': <genotype_object>,
              'mother': <genotype_object>,
              'proband': <genotype_object>}
```
The genotype class have the following attributes for phrasing common questions:

    - genotype STRING (Same as in VCF-standard)
    - allele_1 STRING (Base on allele 1)
    - allele_2 STRING (Base on allele 2)
    - nocall BOOL
    - heterozygote BOOL 
    - homo_alt BOOL (If individual is homozygote alternative)
    - homo_ref BOOL (If individual is homozygote reference)
    - has_variant BOOL (If individual is called and not homozygote reference)
    - ref_depth INT
    - alt_depth INT
    - phred_likelihoods LIST with FLOAT
    - depth_of_coverage INT
    - genotype_quality FLOAT
    - phased BOOL

Vep information, if present, is parsed into

    variant['vep_dict']

and looks like (depending on how vep was run):

    'vep_info': {<alternative_allele>: {
                            'Allele': 'G',
                            'Amino_acids': '',
                            'CDS_position': '',
                            'Codons': '',
                            'Consequence': 'downstream_gene_variant',
                            'DISTANCE': '3084',
                            'EXON': '',
                            'Existing_variation': '',
                            'Feature': 'NM_015658.3',
                            'Feature_type': 'Transcript',
                            'Gene': '26155',
                            'HGVSc': '',
                            'HGVSp': '',
                            'INTRON': '',
                            'PolyPhen': '',
                            'Protein_position': '',
                            'SIFT': '',
                            'STRAND': '-1',
                            'SYMBOL': 'NOC2L',
                            'SYMBOL_SOURCE': '',
                            'cDNA_position': ''},
                  <alternative_allele>: {
                              'Allele': 'G',
                             'Amino_acids': '',
                             'CDS_position': '',
                             'Codons': '',
                             'Consequence': 'intron_variant',
                             'DISTANCE': '',
                             'EXON': '',
                             'Existing_variation': '',
                             'Feature': 'NM_152486.2',
                             'Feature_type': 'Transcript',
                             'Gene': '148398',
                             'HGVSc': 'NM_152486.2:c.707-25A>G',
                             'HGVSp': '',
                             'INTRON': '7/13',
                             'PolyPhen': '',
                             'Protein_position': '',
                             'SIFT': '',
                             'STRAND': '1',
                             'SYMBOL': 'SAMD11',
                             'SYMBOL_SOURCE': '',
                             'cDNA_position': ''
                         }
                    'gene_ids':set([SAMD1, NOC2L])
                     }

INFO field is parsed into a dictionary
The keys are the names of the info field and values are lists separated on ','.

    variant['info_dict]

and looks like

    'info_dict': {'AC': ['1'],
                   'AF': ['0.167'],
                   'AN': ['6'],
                   'BaseQRankSum': ['2.286'],
                   'DB': [],
                   'DP': ['1306'],
                   'FS': ['1.539'],
                   'InbreedingCoeff': ['0.1379'],
                   'MQ': ['39.83'],
                   'MQ0': ['0'],
                   'MQRankSum': ['-2.146'],
                   'POSITIVE_TRAIN_SITE': [],
                   'QD': ['29.57'],
                   'ReadPosRankSum': ['0.897'],
                   'VQSLOD': ['4.52'],
                   'culprit': ['FS'],
                   'set': ['variant']}


### Print a vcf in it´s original format: ###

    my_parser = parser.VCFParser(infile='infile.vcf')
    for line in my_parser.metadata.print_header():
        print(line)
    for variant in my_parser:
	    print('\t'.join([[variant[head] for head in my_parser.header]))

###Add metadata information:###

Adding INFO field:

        my_parser.metadata.add_info('my_new_id', <number>, <type>, <description>)

Where 'number', 'type' and 'description' follows the [VCF](http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41) specification.  

