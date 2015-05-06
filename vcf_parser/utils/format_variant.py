from __future__ import unicode_literals
from logging import getLogger

from vcf_parser import Genotype
from vcf_parser.utils import (build_info_dict, build_vep_annotation, 
build_compounds_dict, build_rank_score_dict, build_models_dict)

def is_number(s):
    """
    Take a string and determin if it is a number
    
    Arguments:
        s (str): A string
    
    Returns:
        bool: True if it is a number, False otherwise
    """
    try:
        float(s)
        return True
    except ValueError:
        return False
    

def format_variant(line, header_parser, skip_info_check=False):
    """
    Yield the variant in the right format. 
    
    If the variants should be splitted on alternative alles one variant 
    for each alternative will be yielded.
    
    Arguments:
        line (str): A string that represents a variant line in the vcf format
    
    Yields:
        variant (dict): A dictionary with the variant information. The number
                        of variants yielded depends on if split variant is used
                        and how many alternatives there are
    """
    logger = getLogger(__name__)

    individuals = []

    vcf_header = header_parser.header

    individuals = header_parser.individuals

    variant_line = line.rstrip().split('\t')

    logger.debug("Checking if variant line is malformed")
    if len(vcf_header) != len(variant_line):
        raise SyntaxError("One of the variant lines is malformed: {0}".format(
            line
        ))

    variant = dict(zip(vcf_header, variant_line))

    # A dictionary with the vep information
    variant['vep_info'] = {}
    # A dictionary with the genetic models (family ids as keys)
    variant['genetic_models'] = {}
    # A dictionary with genotype objects (individual ids as keys)
    variant['genotypes'] = {}
    # A dictionary with the compounds (family ids as keys)
    variant['compound_variants'] = {}
    # A dictionary with the rank scores (family ids as keys)
    variant['rank_scores'] = {}
    
    variant['individual_scores'] = {}
    
    alternatives = variant['ALT'].split(',')
    
    info_dict = build_info_dict(variant.get('INFO', ''))
    
    # Check that the entry is on the proper format_
    if not skip_info_check:
        for info in info_dict:
            annotation = info_dict[info]
            extra_info = header_parser.extra_info.get(info, None)
            if not extra_info:
                logger.warning("The INFO field {0} is not specified in vcf"\
                " header. {1}".format(info, line))
            else:
                number = extra_info['Number']
                if is_number(number):
                    number_of_entrys = float(number)
                    if number_of_entrys != 0:
                        if len(annotation) != number_of_entrys:
                            raise SyntaxError("Info field {0} has the wrong "\
                            "number of entrys according to the vcf header."\
                            "Vcf header line: {1}".format(
                                '='.join([info, ','.join(annotation)]), 
                                header_parser.extra_info.get(info, None)
                            ))
                elif number == 'A':
                    if len(annotation) != len(alternatives):
                        raise SyntaxError("Info field {0} has the wrong "\
                        "number of entrys according to the vcf header"\
                        "Vcf header line: {1}".format(
                            '='.join([info, ','.join(annotation)]),
                            header_parser.extra_info.get(info, None)
                        ))
                elif number == 'R':
                    if len(annotation) != (len(alternatives) + 1):
                        raise SyntaxError("Info field {0} has the wrong "\
                        "number of entrys according to the vcf header".format(
                            '='.join([info, ','.join(annotation)])
                        ))
                elif number == 'G':
                    if len(annotation) != len(individuals):
                        raise SyntaxError("Info field {0} has the wrong "\
                        "number of entrys according to the vcf header".format(
                            '='.join([info, ','.join(annotation)])
                        ))
                        
                        
                
        
    
    variant['info_dict'] = info_dict
    #################### Some fields require special parsing ###########################
    
    ##### VEP ANNOTATIONS #####
    if 'CSQ' in info_dict:
        vep_columns = header_parser.vep_columns
        variant['vep_info'] = build_vep_annotation(
                    info_dict['CSQ'], 
                    variant['REF'], 
                    alternatives,
                    vep_columns
                )
    
    ##### GENMOD ANNOTATIONS #####
    
    if 'GeneticModels' in info_dict:
        variant['genetic_models'] = build_models_dict(
                        info_dict['GeneticModels'])

    if 'Compounds' in info_dict:
        variant['compound_variants'] = build_compounds_dict(
                         info_dict['Compounds'])

    if 'RankScore' in info_dict:
        variant['rank_scores'] = build_rank_score_dict(
                            info_dict['RankScore'])
    
    if 'IndividualRankScore' in info_dict:
        variant['individual_scores'] = build_rank_score_dict(
                                    info_dict['IndividualRankScore'])
    
    ##### GENOTYPE ANNOTATIONS #####
    
    gt_format = variant.get('FORMAT', '').split(':')
    
    genotype_dict = {}
    for individual in individuals:
        gt_info = variant[individual].split(':')
        gt_call = dict(zip(gt_format, gt_info))
        
        #Create a genotype object for this individual
        genotype_dict[individual] = Genotype(**gt_call)
    
    variant['genotypes'] = genotype_dict
    
    variant['variant_id'] = '_'.join(
                                [
                                    variant['CHROM'],
                                    variant['POS'],
                                    variant['REF'],
                                    alternatives[0]
                                ]
                            )
    
    return variant
