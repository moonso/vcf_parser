from logging import getLogger

def build_compounds_dict(compounds):
    """
    Take a list with annotated compound variants for each family and 
    returns a dictionary with family_id as key and a list of dictionarys
    that holds the information about the compounds.
    
    Args:
        compounds    : A list that can be either on the form 
                        [
                            '1:1_23_A_C|1_24_T_A',
                            '2:1_24_T_A'
                        ]
                        or if the compounds are scored:
                        [
                            '1:1_23_A_C>24|1_24_T_A>19',
                            '2:1_24_T_A>17'
                        ]
    
    Returns:
        parsed_compounds : A dictionary on the form
                                {
                                    1:[
                                        {
                                            'variant_id':'1_23_A_C',
                                            'compound_score':24
                                        },
                                        {
                                            'variant_id':'1_24_T_A',
                                            'compound_score:'19
                                        },
                                    ],
                                    2:[
                                        {'variant_id':'1_24_T_A',
                                         'compound_score':17
                                        }
                                    ]
                                }
    
    """
    logger = getLogger(__name__)
    logger.debug("Parsing compounds: {0}".format(compounds))
    
    parsed_compounds = {}
    for family_info in compounds:
        logger.debug("Parsing entry {0}".format(family_info))
        splitted_family_info = family_info.split(':')
        family_id = splitted_family_info[0]
        logger.debug("Found family {0}".format(family_id))
        parsed_compounds[family_id] = []
        compound_list = splitted_family_info[1].split('|')
        for compound in compound_list:
            compound_id = compound.split('>')[0]
            try:
                compound_score = compound.split('>')[1]
            except IndexError:
                compound_score = None
            parsed_compounds[family_id].append(
                {
                    'variant_id': compound_id,
                    'compound_score': compound_score
                }
            )
    
    return parsed_compounds
