from __future__ import unicode_literals
from logging import getLogger

from vcf_parser import Genotype
from vcf_parser.utils import (build_info_dict, build_vep_annotation, 
build_compounds_dict, build_rank_score_dict)

def format_variant(line, vcf_header, vep_header = []):
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
    if len(vcf_header) > 9:
        individuals = vcf_header[9:]
    
    variant_line = line.rstrip().split('\t')
    variant = dict(zip(vcf_header, line.rstrip().split('\t')))
    
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
    
    variant['info_dict'] = info_dict
    #################### Some fields require special parsing ###########################
    
    ##### VEP ANNOTATIONS #####
    if 'CSQ' in info_dict:
        variant['vep_info'] = build_vep_annotation(
                    info_dict['CSQ'], 
                    variant['REF'], 
                    alternatives
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
