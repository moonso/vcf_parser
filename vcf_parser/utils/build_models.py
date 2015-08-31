from logging import getLogger

def build_models_dict(annotated_models):
    """
    Take a list with annotated genetic inheritance patterns for each
    family and returns a dictionary with family_id as key and a list of
    genetic models as value.
    
    Args:
        annotated_models    : A list on the form ['1:AD','2:AR_comp|AD_dn']
    
    Returns:
        parsed_models       : A dictionary on the form
                                {
                                    1:['AD'],
                                    2:['AD_dn','AR_comp']
                                }
    
    """
    logger = getLogger(__name__)
    logger.debug("Parsing models {0}".format(annotated_models)
    )
    parsed_models = {}
    for family_annotation in annotated_models:
        family_id = family_annotation.split(':')[0]
        logger.debug("Parsing family {0}".format(family_id))
        models = family_annotation.split(':')[1].split('|')
        parsed_models[family_id] = models
        logger.debug("Adding models {0}".format(models))
    
    return parsed_models
