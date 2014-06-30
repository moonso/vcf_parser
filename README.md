#VCF Parser#
Small library for parsing vcf files. Based on [PyVCF](https://github.com/jamescasbon/PyVCF)

    pip install vcf_parser

```python3
    from vcf_parser import parser
    my_parser = parser.VCFParser('infile.vcf')
    for variant in my_parser:
        print(variant)
```

Returns dictionary with the vcf info for each variant.
The genotype information is converted to a genotype object and stored in a dictionary

    variant['ind_dict']

and looks like:

    'ind_dict': {'father': genotype_object,
                  'mother': genotype_object,
                  'proband': genotype_object}

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

Vep information, if present is parsed into

    variant['vep_dict']

and looks like:

    'vep_info': {'NOC2L': {'Allele': 'G',
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
                  'SAMD11': {'Allele': 'G',
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
                             'cDNA_position': ''}}

INFO field is parsed into 

    variant['info_dict]

and looks like

    'info_dict': {'AC': '1',
                   'AF': '0.167',
                   'AN': '6',
                   'BaseQRankSum': '2.286',
                   'DB': True,
                   'DP': '1306',
                   'FS': '1.539',
                   'InbreedingCoeff': '0.1379',
                   'MQ': '39.83',
                   'MQ0': '0',
                   'MQRankSum': '-2.146',
                   'POSITIVE_TRAIN_SITE': True,
                   'QD': '29.57',
                   'ReadPosRankSum': '0.897',
                   'VQSLOD': '4.52',
                   'culprit': 'FS',
                   'set': 'variant'}

