from logging import getLogger

def split_genotype(genotype, gt_format, alternative_number, allele_symbol = '0'):
    """
    Take a genotype call and make a new one that is working for the new
    splitted variant
    
    Arguments:
        genotype (str): The original genotype call
        gt_format (str): The format of the gt call
        alternative_number (int): What genotype call should we return
        allele_symbol (str): How should the unobserved allele be represented
                             when genotype is splitted
    
    Returns:
        new_genotype (str): A string that represents the new genotype
    """
    logger = getLogger(__name__)
    
    logger.info("Allele symbol {0}".format(allele_symbol))

    splitted_genotype = genotype.split(':')
    logger.debug("Parsing genotype {0}".format(splitted_genotype))
    splitted_gt_format = gt_format.split(':')
    logger.debug("Parsing gt format {0}".format(splitted_gt_format))
    new_genotype = []
    phased = False
    for number, genotype_info in enumerate(splitted_genotype):
        gt_info = splitted_gt_format[number]
        if gt_info == 'GT':
            if '/' in genotype_info:
                gt = genotype_info.split('/')
            else:
                gt = genotype_info.split('|')
                phased = True
            ref_allele = '.'
            alt_allele = '.'
            try:
                # Check the ref Allele
                if len(gt) == 2 and gt[0] != '.' and gt[1] != '.':
                    ref_allele = allele_symbol
                    alt_allele = allele_symbol
                    if gt[0] == gt[1]:
                        # In this case we have a homozygous call:
                        if int(gt[0]) == alternative_number + 1:
                            ref_allele = '1'
                            alt_allele = '1'
                    else:
                        if (int(gt[0]) == alternative_number + 1 or 
                            int(gt[1]) == alternative_number + 1):
                            alt_allele = '1'                        
                else:
                # We now know that at least one of the alleles are uncalled
                    if gt[0] != '.':
                        if int(gt[0]) == alternative_number + 1:
                            ref_allele = '1'
                        else:
                            ref_allele = '0'
                    elif len(gt) == 2 and gt[1] != '.':
                        if int(gt[1]) == alternative_number + 1:
                            alt_allele = '1'
                        else:
                            alt_allele = '0'
            except (ValueError, KeyError):
                pass
            
            if len(gt) == 2:
                if phased:
                    new_genotype.append('|'.join([ref_allele,alt_allele]))
                else:
                    new_genotype.append('/'.join([ref_allele,alt_allele]))
            else:
                new_genotype.append(ref_allele)
        
        elif gt_info == 'AD':
            ad = []
            # The reference depth will allways be the original depth now
            ad.append(genotype_info.split(',')[0])
            try:
                ad.append(genotype_info.split(',')[alternative_number+1])
            except IndexError:
                ad.append('0')
            new_genotype.append(','.join(ad))
        elif gt_info == 'DP':
            new_genotype.append(genotype_info)
        elif gt_info == 'PL':
            new_genotype.append(genotype_info)
        else:
            # There are several cases that we do not know how to handle yet so we just add the information
            new_genotype.append(genotype_info)
            
    return ':'.join(new_genotype)
