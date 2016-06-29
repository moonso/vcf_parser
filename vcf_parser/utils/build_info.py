import sys
import logging

if sys.version_info < (2, 7):
    from ordereddict import OrderedDict
else:
    from collections import OrderedDict


def build_info_string(info):
    """
    Build a new vcf INFO string based on the information in the info_dict.
    
    The info is a dictionary with vcf info keys as keys and lists of vcf values
    as values. If there is no value False is value in info
    
    Args:
        info (dict): A dictionary with information from the vcf file
    
    Returns:
        String: A string that is on the proper vcf format for the INFO column
    
    """
    info_list = []
    
    for annotation in info:
        
        if info[annotation]:
            info_list.append('='.join([annotation, ','.join(info[annotation])]))
        else:
            info_list.append(annotation)
    
    return ';'.join(info_list)

def build_info_dict(vcf_info):
    """
    Build a dictionary from the info of a vcf line
    
    The dictionary will have the info keys as keys and info values as values.
    Values will allways be lists that are splitted on ','
    
    Arguments:
        vcf_info (str): A string with vcf info
    
    Returns:
        info_dict (OrderedDict): A ordered dictionary with the vcf info keys as 
                                 keys and lists of values as values
    """
    logger = logging.getLogger(__name__)
    logger.debug("Building info dict")
    info_dict = OrderedDict()
    
    for info in vcf_info.split(';'):
        info = info.split('=')
        if len(info) > 1:
            # If the INFO entry is like key=value, we store the value as a list
            info[1] = '='.join(info[1:])
            info_dict[info[0]] = info[1].split(',')
        else:
            info_dict[info[0]] = []
    
    return info_dict
    