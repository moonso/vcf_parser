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


def check_info_annotation(annotation, info, extra_info, alternatives, individuals=[]):
    """
    Check if the info annotation corresponds to the metadata specification
    
    Arguments:
        annotation (list): The annotation from the vcf file
        info (str): Name of the info field
        extra_info (dict): The metadata specification
        alternatives (list): A list with the alternative variants
        individuals (list): a list with the individuals
    
    Returns:
        bool: If the annotation is correct or not
    """
    
    number = extra_info['Number']
    if is_number(number):
        number_of_entrys = float(number)
        if number_of_entrys != 0:
            if len(annotation) != number_of_entrys:
                raise SyntaxError("Info field {0} has the wrong "\
                "number of entrys according to the vcf header."\
                " Vcf header specifies {1} should have {2} entry(s)".format(
                    '='.join([info, ','.join(annotation)]), 
                    info,
                    number
                ))
    elif number == 'A':
        if len(annotation) != len(alternatives):
            raise SyntaxError("Info field {0} has the wrong "\
            "number of entrys according to the vcf header."\
            "Vcf header specifies {1} should have {2} entry(s)".format(
                    '='.join([info, ','.join(annotation)]), 
                    info,
                    number
            ))
    elif number == 'R':
        if len(annotation) != (len(alternatives) + 1):
            raise SyntaxError("Info field {0} has the wrong "\
            "number of entrys according to the vcf header."\
            "Vcf header specifies {1} should have {2} entry(s)".format(
                    '='.join([info, ','.join(annotation)]), 
                    info,
                    number
            ))
    elif number == 'G':
        if len(annotation) != len(individuals):
            raise SyntaxError("Info field {0} has the wrong "\
            "number of entrys according to the vcf header."\
            "Vcf header specifies {1} should have {2} entry(s)".format(
                    '='.join([info, ','.join(annotation)]), 
                    info,
                    number
            ))
    return True
