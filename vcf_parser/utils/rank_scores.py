from logging import getLogger

def build_rank_score_dict(rank_scores):
    """
    Take a list with annotated rank scores for each family and returns a 
    dictionary with family_id as key and a list of genetic models as value.
    
    Args:
        rank_scores    : A list on the form ['1:12','2:20']
    
    Returns:
        scores       : A dictionary with family id:s as key and scores as value
                                {
                                    '1':'12',
                                    '2':'20'
                                }
    
    """
    logger = getLogger(__name__)
    logger.debug("Checking rank scores: {0}".format(rank_scores))
    scores = {}
    for family in rank_scores:
        entry = family.split(':')
        try:
            family_id = entry[0]
            logger.debug("Extracting rank score for family:{0}".format(family_id))
            score = entry[1]
            logger.debug("Score:{0}".format(score))
        except Exception:
            raise SyntaxError("Malformed rank score input")
            
        scores[family_id] = score
    
    return scores
