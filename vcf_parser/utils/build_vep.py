from logging import getLogger

def build_vep_string(vep_info, vep_columns):
    """
    Build a vep string formatted string.
    
    Take a list with vep annotations and build a new vep string
    
    Args:
        vep_info (list): A list with vep annotation dictionaries
        vep_columns (list): A list with the vep column names found in the
        header of the vcf
    
    Returns:
        string: A string with the proper vep annotations
    
    """
    logger = getLogger(__name__)
    logger.debug("Building vep string from {0}".format(vep_info))
    logger.debug("Found vep headers {0}".format(vep_columns))
    vep_strings = []
    for vep_annotation in vep_info:
        try:
            vep_info_list = [
                vep_annotation[vep_key] for vep_key in vep_columns
            ]
        except KeyError:
            raise SyntaxError("Vep entry does not correspond to vep headers")
        
        vep_strings.append('|'.join(vep_info_list))
    return ','.join(vep_strings)

def build_vep_annotation(csq_info, reference, alternatives, vep_columns):
    """
    Build a dictionary with the vep information from the vep annotation.
    
    Indels are handled different by vep depending on the number of 
    alternative alleles there is for a variant.
    
    If only one alternative:
        
        Insertion: vep represents the alternative by removing the first 
        base from the vcf alternative.
        
        Deletion: vep represents the alternative with '-'
    
    If there are several alternatives:
        
        Insertion: 
        vep represents the alternative by removing the first 
        base from the vcf alternative(Like above).
        
        Deletion: 
        If there are multiple alternative deletions vep represents them by 
        removing the first base from the vcf alternative.
        If the vcf line looks like:
            1   970549  .   TGGG    TG,TGG
        vep annotation for alternatives will be: G,GG
    
    Args:
        csq_info (list): A list with the raw vep annotations from the vcf line.
        reference (str): A string that represents the vcf reference
        alternatives (list): A list of strings that represents the vcf formated
                             alternatives
        vep_columns (list): A list of strings that represents the vep comluns
                            defined in the vcf header.
    
    Returns:
        vep_dict (dict): A dictionary with the alternative alleles (in vcf form)
                         as keys and a list of annotations for each alternative 
                         alleles. 
                         One key named 'gene_ids', 
                         value is a set with the genes found. 
    """
    logger = getLogger(__name__)

    # The keys in the vep dict are the vcf formatted alternatives, values are the
    # dictionaries with vep annotations
    vep_dict = {}

    # If we have several alternatives we need to check what types of 
    # alternatives we have
    vep_to_vcf = {}
    number_of_deletions = 0
    for alternative in alternatives:
        if len(alternative) < len(reference):
            number_of_deletions += 1

    logger.debug("Number of deletions found: {0}".format(number_of_deletions))
    for alternative in alternatives:
        # We store the annotations with keys from the vcf alternatives
        vep_dict[alternative] = []

        # If substitutuion reference and alternative have the same length
        if len(alternative) == len(reference):
             vep_to_vcf[alternative] = alternative
        # If deletion alternative is shorter that the reference
        else:
            # If there is a deletion then the alternative will be '-' in vep entry
            if len(alternative) == 1:
                vep_to_vcf['-'] = alternative
            else:
                vep_to_vcf[alternative[1:]] = alternative

    for vep_annotation in csq_info:
        logger.debug("Parsing vep annotation: {0}".format(vep_annotation))
        splitted_vep = vep_annotation.split('|')
        
        if len(splitted_vep) != len(vep_columns):
            raise SyntaxError("Csq info for variant does not match csq info in "\
                            "header. {0}, {1}".format(
            '|'.join(splitted_vep), '|'.join(vep_columns)))
        
        # Build the vep dict:
        vep_info = dict(zip(vep_columns, splitted_vep))
        
        # If no allele is found we can not determine what allele
        if vep_info.get('Allele', None):
            vep_allele = vep_info['Allele']
            try:
                vcf_allele = vep_to_vcf[vep_allele]
            except KeyError as e:
                vcf_allele = vep_allele
    
            if vcf_allele in vep_dict:
                vep_dict[vcf_allele].append(vep_info)
            else:
                vep_dict[vcf_allele] = [vep_info]
        else:
            logger.warning("No allele found in vep annotation! Skipping annotation")
            


    return vep_dict
