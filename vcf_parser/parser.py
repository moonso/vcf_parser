#!/usr/bin/env python
# encoding: utf-8
"""
vcf_parser.py


Parse a vcf file.

Includes a header class for storing information about the headers.
Create variant objects and a dictionary with individuals that have a dictionary with genotypes for each variant.

Thanks to PyVCF for heaader parser and more...:

Copyright (c) 2011-2012, Population Genetics Technologies Ltd, All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this
list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the Population Genetics Technologies Ltd nor the names of
its contributors may be used to endorse or promote products derived from this
software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


Copyright (c) 2011 John Dougherty

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.



Created by Måns Magnusson on 2013-01-17.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

from __future__ import print_function, unicode_literals

import sys
import os
import gzip
import re
import pkg_resources
import click
import locale

from codecs import open, getreader, getwriter


if sys.version_info < (2, 7):
    from ordereddict import OrderedDict
else:
    from collections import OrderedDict

from pprint import pprint as pp

from vcf_parser import genotype
import vcf_parser

# Version = pkg_resources.require("vcf_parser")[0].version
Version = '0.7.3'

class HeaderParser(object):
    """Parses a file with family info and creates a family object with individuals."""
    def __init__(self):
        super(HeaderParser, self).__init__()
        self.info_lines=[]
        self.info_dict=OrderedDict()
        #This is a dictionary cantaining specific information about the info fields
        #It will have info name as key and then another dictionary with ID, Number, Type and Description
        self.extra_info = {}
        
        self.filter_lines=[]
        self.filter_dict=OrderedDict()
        
        self.contig_lines=[]
        self.contig_dict=OrderedDict()
        
        self.format_lines=[]
        self.format_dict=OrderedDict()
        
        self.alt_lines=[]
        self.alt_dict=OrderedDict()
        
        self.other_lines=[]
        self.other_dict=OrderedDict()
        
        self.header=[]
        self.header_keys={'info' : ['ID', 'Number', 'Type', 'Description'], 
                            'form' : ['ID', 'Number', 'Type', 'Description'], 
                            'filt' : ['ID', 'Description'],
                            'alt' : ['ID', 'Description'],
                            'contig' : ['ID', 'length']}
        self.fileformat = ''
        self.line_counter = 0
        self.individuals = []
        self.vep_columns = []
        self.info_pattern = re.compile(r'''\#\#INFO=<
            ID=(?P<id>[^,]+),
            Number=(?P<number>-?\d+|\.|[AG]),
            Type=(?P<type>Integer|Float|Flag|Character|String),
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE)
        self.filter_pattern = re.compile(r'''\#\#FILTER=<
            ID=(?P<id>[^,]+),
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE)
        self.contig_pattern = re.compile(r'''\#\#contig=<
            ID=(?P<id>[^,]+),
            .*
            length=(?P<length>-?\d+)
            .*
            >''', re.VERBOSE)
        self.format_pattern = re.compile(r'''\#\#FORMAT=<
            ID=(?P<id>.+),
            Number=(?P<number>-?\d+|\.|[AG]),
            Type=(?P<type>.+),
            Description="(?P<desc>.*)"
            >''', re.VERBOSE)
        self.alt_pattern = re.compile(r'''\#\#ALT=<
            ID=(?P<id>[^,]+),
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE)
        self.meta_pattern = re.compile(r'''##(?P<key>.+?)=(?P<val>.+)''')
    
    def parse_meta_data(self, line):
        """Parse a vcf metadataline"""
        line = line.rstrip()
        line_info = line[2:].split('=')
        match = False
        
        if line_info[0] == 'fileformat':
            self.fileformat = line_info[1]
        elif line_info[0] == 'INFO':
            match = self.info_pattern.match(line)
            if not match:
                print("\nOne of the INFO lines is malformed:\n%s\n" % line)
                raise SyntaxError()
            matches = [match.group('id'), match.group('number'), match.group('type'), match.group('desc')]
            # extra_info is a dictionary to check the metadata about the INFO values:
            self.extra_info[matches[0]] = dict(zip(self.header_keys['info'][1:], matches[1:]))
            info_line = dict(list(zip(self.header_keys['info'],matches)))
            if len(info_line['Description'].split('Format:')) > 1:
                info_line['Format'] = [info.strip() for info in info_line['Description'].split('Format:')][-1]
            self.info_lines.append(info_line)
            # Store the vep columns:
            if info_line['ID'] == 'CSQ':
                self.vep_columns = info_line.get('Format', '').split('|')
            self.info_dict[match.group('id')] = line
        elif line_info[0] == 'FILTER':
            match = self.filter_pattern.match(line)
            if not match:
                print("\nOne of the FILTER lines is malformed:\n%s\n" % line)
                raise SyntaxError()
            matches = [match.group('id'), match.group('desc')]
            self.filter_lines.append(dict(list(zip(self.header_keys['filt'],matches))))
            self.filter_dict[match.group('id')] = line
        elif line_info[0] == 'contig':
            match = self.contig_pattern.match(line)
            if not match:
                print("One of the contig lines is malformed:\n %s" % line)
                raise SyntaxError()
            matches = [match.group('id'), match.group('length')]
            self.contig_lines.append(dict(list(zip(self.header_keys['contig'],matches))))
            self.contig_dict[match.group('id')] = line
        elif line_info[0] == 'FORMAT':
            match = self.format_pattern.match(line)
            if not match:
                print("\nOne of the FORMAT lines is malformed:\n%s\n" % line)
                raise SyntaxError()
            matches = [match.group('id'), match.group('number'), match.group('type'), match.group('desc')]
            self.format_lines.append(dict(list(zip(self.header_keys['form'],matches))))
            self.format_dict[match.group('id')] = line
        elif line_info[0] == 'ALT':
            match = self.alt_pattern.match(line)
            if not match:
                print("\nOne of the ALT lines is malformed:\n%s\n" % line)
                raise SyntaxError()
            matches = [match.group('id'), match.group('desc')]
            self.alt_lines.append(dict(list(zip(self.header_keys['alt'],matches))))
            self.alt_dict[match.group('id')] = line
        else:
            match = self.meta_pattern.match(line)
            if not match:
                print("\nOne of the meta data lines is malformed:\n%s\n" % line)
                raise SyntaxError()
            self.other_lines.append({match.group('key'): match.group('val')})
            self.other_dict[match.group('key')] = line
    
    def parse_header_line(self, line):
        """docstring for parse_header_line"""
        self.header = line[1:].rstrip().split('\t')
        if len(self.header) < 9:
            self.header = line[1:].rstrip().split()
        self.individuals = self.header[9:]
    
    def print_header(self):
        """Returns a list with the header lines if proper format"""
        lines_to_print = []
        lines_to_print.append('##fileformat='+self.fileformat)
        for filt in self.filter_dict:
            lines_to_print.append(self.filter_dict[filt])
        for form in self.format_dict:
            lines_to_print.append(self.format_dict[form])
        for info in self.info_dict:
            lines_to_print.append(self.info_dict[info])
        for contig in self.contig_dict:
            lines_to_print.append(self.contig_dict[contig])
        for alt in self.alt_dict:
            lines_to_print.append(self.alt_dict[alt])
        for other in self.other_dict:
            lines_to_print.append(self.other_dict[other])
        lines_to_print.append('#'+ '\t'.join(self.header))
        return lines_to_print
    
    def add_info(self, info_id, number, entry_type, description):
        """Add an info line to the header."""
        info_line = '##INFO=<ID='+info_id+',Number='+str(number)+',Type='+entry_type+',Description="'+description+'">'
        self.info_dict[info_id] = info_line
        return

    def add_version_tracking(self, info_id, version, date, command_line=''):
        """Add a line with information about which software that was run and when to the header."""
        other_line = ['##Software=<ID=%s' % info_id ,'Version=%s' % version, 'Date="%s"' % date, 'CommandLineOptions="%s">' % command_line]
        self.other_dict[info_id] = ','.join(other_line)
        return


####            Parser:         ####


class VCFParser(object):
    """docstring for VCFParser"""
    def __init__(self, infile=None, fsock=None, split_variants=False):
        super(VCFParser, self).__init__()
        
        if not (fsock or infile):
            raise Exception('You must provide at least fsock or filename')
        
        # print('Hello' ,infile, type(infile), fsock, type(fsock))
        if fsock:
            if not infile and hasattr(fsock, 'name'):
                if sys.version_info < (3, 0):
                    sys.stdin = getreader('utf-8')(fsock)            
                self.vcf = sys.stdin
        
        else:
            file_name, file_extension = os.path.splitext(infile)
            if file_extension == '.gz':
                self.vcf = getreader('utf-8')(gzip.open(infile), errors='replace')
            elif file_extension == '.vcf':
                self.vcf = open(infile, mode='r', encoding='utf-8', errors='replace')
            else:
                print("""File is not in a supported format!\n Or use correct ending(.vcf or .vcf.gz)""")
                raise SyntaxError()
        
        self.split_variants = split_variants
        self.metadata = HeaderParser()
        self.individuals = []
        self.header = []
        
        self.next_line = self.vcf.readline().rstrip()
        self.current_line = self.next_line
        self.metadata.parse_meta_data(self.next_line)
        self.beginning = True
        
        while self.next_line.startswith('#'):
            if self.next_line.startswith('##'):
                self.metadata.parse_meta_data(self.next_line)
            elif self.next_line.startswith('#'):
                self.metadata.parse_header_line(self.next_line)
            self.next_line = self.vcf.readline().rstrip()
        self.individuals = self.metadata.individuals
        self.header = self.metadata.header
            
    def __iter__(self):
        for line in self.vcf:
            line = line.rstrip()
            variants = []
            if self.beginning:
                first_variant = self.format_variant(self.next_line)
                if not (self.split_variants and len(first_variant['ALT'].split(',')) > 1):
                    variants.append(first_variant)
                else:
                    for variant in self.make_splitted_variants(first_variant):
                        variants.append(variant)
                
                self.beginning = False
            
            if len(line.split('\t')) >= 8:
                variant = self.format_variant(line)
                if not (self.split_variants and len(variant['ALT'].split(',')) > 1):
                    variants.append(variant)
                
                else:
                    for splitted_variant in self.make_splitted_variants(variant):
                        variants.append(splitted_variant)
            
            for variant in variants:
                yield variant
    
    def build_new_info_string(self, info_dict):
        """
        Build a new INFO string based on the information in the info_dict.
        
        Args:
            info_dict: A dictionary with information from the vcf file
        
        Returns:
            String: A string that is on the proper vcf format for the INFO column
        """
        info_list = []
        for info in info_dict:
            if info_dict[info]:
                info_list.append('='.join([info, ','.join(info_dict[info])]))
            else:
                info_list.append(info)
        return ';'.join(info_list)
            
    def build_vep_annotation(self, csq_info, reference, alternatives):
        """
        Build a dictionary with the vep information from the vep annotation.
        Indels are handled different by vep depending on the number of 
        alternative alleles there is for a variant.
        If only one alternative:
            
            Insertion: vep represents the alternative by removing the first 
            base from the vcf alternative.
            
            Deletion: vep represents the alternative with '-'
        
        If there are several alternatives:
            
            Insertion: vep represents the alternative by removing the first 
            base from the vcf alternative(Like above).
            
            Deletion: If there is only one deletion among the alternatives it 
            will allways be represented with '-'.
            If there are multiple alternative deletions vep represents them by 
            removing the first base from the vcf alternative.
            If the vcf line looks like:
                1   970549  .   TGGG    TG,TGG
            vep annotation for alternatives will be: G,GG
        
        Args:
            csq_info: A list with the vep annotations from the vcf line.
        
        Returns:
            vep_dict (dict): A dictionary with the alternative alleles as keys 
            and a list of annotations for each alternative alleles. 
            One key named gene_ids, value is a set with the genes found. 
        """
        # These are the vep terms for the insertions and deletions
        VEP_INSERTION = 'feature_elongation'
        VEP_DELETION = 'feature_truncation'
        vep_dict = {'gene_ids' : set([])}
        deletions = []
        vep_deletions = []
        insertions = []
        vep_insertions = []
        substitutions = []
        # If we have several alternatives we need to check what types of alternatives we have
        for alternative in alternatives:
            vep_dict[alternative] = []
            if len(alternative) == len(reference):
                 substitutions.append(alternative)
            elif len(alternative) < len(reference):
                deletions.append(alternative)
                if len(alternative) == 1:
                    vep_deletions.append('-')
                else:
                    vep_deletions.append(alternative[1:])
            elif len(alternative) > len(reference):
                insertions.append(alternative)
                vep_insertions.append(alternative[1:])
        for vep_annotation in csq_info:
            vep_info = dict(zip(self.metadata.vep_columns, vep_annotation.split('|')))
            # If we only have one alternative then all vep annotations will represent that alternative:
            vep_allele = vep_info.get('Allele', '-')
            alternative_allele = vep_allele
            if len(alternatives) == 1:
                vep_dict[alternatives[0]].append(vep_info)
            else:
                consequences = set(vep_info.get('Consequence', '').split('&'))
                # If there is only one deletion it can be represented by '-'
                if vep_allele == '-':
                    alternative_allele = reference[0]
                
                elif VEP_INSERTION in consequences or VEP_DELETION in consequences:
                    alternative_allele = reference[0] + vep_allele
                
                else:
                    if vep_allele in vep_deletions or vep_allele in vep_insertions:
                        alternative_allele = reference[0] + vep_allele
                    
                    else:
                        if vep_allele in substitutions:
                            alternative_allele = vep_allele
            
                if alternative_allele in vep_dict:
                    vep_dict[alternative_allele].append(vep_info)
                else:
                    vep_dict[alternative_allele] = [vep_info]
            
            # Save the gene annotations for this variant:
            vep_dict['gene_ids'].add(vep_info.get('SYMBOL','-'))
        return vep_dict
    
    def build_new_vep_string(self, vep_info):
        """
        Take a list with vep annotations and build a new vep string
        
        Args:
            vep_info (list): A list with vep annotations
        
        Returns:
            string: A string with the proper vep annotations
        
        """
        vep_strings = []
        for vep_annotation in vep_info:
            vep_info_list = [vep_annotation[vep_key] for vep_key in self.metadata.vep_columns]
            vep_strings.append('|'.join(vep_info_list))
        return ','.join(vep_strings)
    
    def split_genotype(self, genotype, gt_format, alternative_number):
        """
        Take a genotype call and make a new one that is working for the new
        splitted variant
        """
        # print(genotype, gt_format, alternative_number)
        splitted_genotype = genotype.split(':')
        splitted_gt_format = gt_format.split(':')
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
                    if gt[0] != '.' and gt[1] != '.':
                        ref_allele = '0'
                        alt_allele = '0'
                        if gt[0] == gt[1]:
                            # In this case we have a homozygous call:
                            if int(gt[0]) == alternative_number + 1:
                                ref_allele = '1'
                                alt_allele = '1'
                        else:
                            if int(gt[0]) == alternative_number + 1 or int(gt[1]) == alternative_number + 1:
                                alt_allele = '1'                        
                    else:
                    # We now know that at least one of the alleles are uncalled
                        if gt[0] != '.':
                            if int(gt[0]) == alternative_number + 1:
                                ref_allele = '1'
                            else:
                                ref_allele = '0'
                        elif gt[1] != '.':
                            if int(gt[1]) == alternative_number + 1:
                                alt_allele = '1'
                            else:
                                alt_allele = '0'
                except (ValueError, KeyError):
                    pass
                
                if phased:
                    new_genotype.append('|'.join([ref_allele,alt_allele]))
                else:
                    new_genotype.append('/'.join([ref_allele,alt_allele]))            
            elif gt_info == 'AD':
                ad = []
                # The reference depth will allways be the original depth now
                ad.append(genotype_info.split(',')[0])
                ad.append(genotype_info.split(',')[alternative_number+1])
                new_genotype.append(','.join(ad))
            elif gt_info == 'DP':
                new_genotype.append(genotype_info)
            elif gt_info == 'PL':
                new_genotype.append(genotype_info)
            else:
                # There are several cases that we do not know how to handle yet so we just add the information
                new_genotype.append(genotype_info)
                
        return ':'.join(new_genotype)
            
    
    def make_splitted_variants(self, variant_dict):
        """
        Checks if there are multiple alternative alleles and splitts the 
        variant.
        If there are multiple alternatives the info fields, vep annotations 
        and genotype calls will be splitted in the correct way
        
        Args:
            variant_dict: a dictionary with the varianinformation
        
        Returns:
            List: A list of variant dictionaries with the splitted information
        """
        variants = []
        alternatives = variant_dict['ALT'].split(',')
        reference = variant_dict['REF']
        # Go through each of the alternative alleles:
        for alternative_number, alternative in enumerate(alternatives):
            variant = {}
            info_dict = OrderedDict()
            # This is a dict on the form {ALT:[{vep_info_dict}]}
            vep_dict = {}
            genotype_dict = {}
            variant['CHROM'] = variant_dict['CHROM']
            variant['POS'] = variant_dict['POS']
            try:
                # There will not allways be one rsID for each alternative
                variant['ID'] = variant_dict['ID'].split(';')[alternative_number]
            # If only one id is present for multiple alleles they all get the same ID
            except IndexError:
                variant['ID'] = variant_dict['ID']
            variant['REF'] = variant_dict['REF']
            variant['ALT'] = alternative
            variant['QUAL'] = variant_dict['QUAL']
            variant['FILTER'] = variant_dict['FILTER']
            gt_format = variant_dict['FORMAT']
            variant['FORMAT'] = gt_format
            
            
            for info in variant_dict['info_dict']:
                if info:
                    # Check if the info field have one entry per allele:
                    try:
                        number_of_values = self.metadata.extra_info[info]['Number']
                    except KeyError:
                        print(""""\nOne of the FILTER lines is missing in vcf 
                                header: %s \n""" % info, file=sys.stderr)
                        raise 
                    # If there if one value per allele we need to split it in
                    # the proper way
                    if number_of_values == 'A':
                        try:
                            # When we split the alleles we only want to annotate with the correct number
                            info_dict[info] = [variant_dict['info_dict'][info][alternative_number]]
                        except IndexError:
                            # If there is only one annotation we choose that one
                            info_dict[info] = [variant_dict['info_dict'][info][0]]
                    # Choose the right vep info from the old variant
                    elif info == 'CSQ':
                        try:
                            vep_dict[alternative] = variant_dict['vep_info'][alternative]
                            info_dict['CSQ'] = [self.build_new_vep_string(variant_dict['vep_info'][alternative])]
                        except KeyError:
                            pass
                    else:
                        info_dict[info] = variant_dict['info_dict'][info]
                    
                else:
                    info_dict[info] = False
            
            variant['INFO'] = self.build_new_info_string(info_dict)
            
            for individual in variant_dict['genotypes']:
                new_genotype = self.split_genotype(variant_dict[individual], gt_format, alternative_number)
                variant[individual] = new_genotype
                genotype_dict[individual] = genotype.Genotype(**dict(zip(gt_format.split(':'), variant[individual].split(':'))))
                
            variant['info_dict'] = info_dict
            variant['vep_info'] = vep_dict
            variant['genotypes'] = genotype_dict
            variant['variant_id'] = '_'.join([variant['CHROM'],
                                        variant['POS'],
                                        variant['REF'],
                                        alternative])
            
            yield variant
            
    def build_models_dict(self, annotated_models):
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
        parsed_models = {}
        for family_annotation in annotated_models:
            family_id = family_annotation.split(':')[0]
            models = family_annotation.split(':')[1].split('|')
            parsed_models[family_id] = models
        
        return parsed_models
    
    def build_rank_score_dict(self, rank_scores):
        """
        Take a list with annotated rank scores for each family and returns a 
        dictionary with family_id as key and a list of genetic models as value.
        
        Args:
            rank_scores    : A list on the form ['1:12','2:20']
        
        Returns:
            scores       : A dictionary on the form
                                    {
                                        1:12,
                                        2:20
                                    }
        
        """
        scores = {}
        for family in rank_scores:
            entry = family.split(':')
            family_id = entry[0]
            score = entry[1]
            scores[family_id] = score
        
        return scores
    
    def build_compounds_dict(self, compounds):
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
                                                'compound_score:24
                                            },
                                            {
                                                'variant_id':'1_24_T_A',
                                                'compound_score:19
                                            },
                                        ],
                                        2:[
                                            {'variant_id':'1_24_T_A',
                                             'compound_score':17
                                            }
                                        ]
                                    }
        
        """
        parsed_compounds = {}
        for family in compounds:
            family_id = family.split(':')[0]
            parsed_compounds[family_id] = []
            compound_list = family.split(':')[1].split('|')
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
    
    def format_variant(self, line):
        """
        Yield the variant in the right format. If the variants should be splitted on alternative alles
        one variant for each alternative will be yielded.
        """
        variant_line = line.rstrip().split('\t')
        variant = dict(zip(self.header, line.rstrip().split('\t')))
        
        info_dict = OrderedDict()
        vep_dict = {}
        models_dict = {}
        genotype_dict = {}
        compunds_dict = {}
        rank_score_dict = {}
        individual_score_dict = {}
        alternatives = variant['ALT'].split(',')
        for info in variant.get('INFO', '').split(';'):
            info = info.split('=')
            if len(info) > 1:
                #If the INFO entry is like key=value, we store the value as a list
                info_dict[info[0]] = info[1].split(',')
            else:
                info_dict[info[0]] = False
        
        #################### Some fields require special parsing ###########################
        
        ##### VEP ANNOTATIONS #####
        if 'CSQ' in info_dict:
            vep_dict = self.build_vep_annotation(
                                            info_dict['CSQ'], 
                                            variant['REF'], 
                                            alternatives
                                            )
        
        ##### GENMOD ANNOTATIONS #####
        
        if 'GeneticModels' in info_dict:
            models_dict = self.build_models_dict(
                                            info_dict['GeneticModels']
                                            )
        if 'Compounds' in info_dict:
            compunds_dict = self.build_compounds_dict(
                                                info_dict['Compounds']
                                                )
        if 'RankScore' in info_dict:
            rank_score_dict = self.build_rank_score_dict(
                                                info_dict['RankScore']
                                                )
        
        if 'IndividualRankScore' in info_dict:
            individual_score_dict = self.build_rank_score_dict(
                                                info_dict['IndividualRankScore']
                                                )
        
        ##### GENOTYPE ANNOTATIONS #####
        
        gt_format = variant.get('FORMAT', '').split(':')
        
        for individual in self.individuals:
            genotype_dict[individual] = genotype.Genotype(
                                                    **dict(
                                                        zip(
                                                        gt_format,
                                                        variant[individual].split(':')
                                                        )
                                                    )
                                                )
        
        variant['genotypes'] = genotype_dict
        variant['info_dict'] = info_dict
        variant['variant_id'] = '_'.join(
                                    [
                                        variant['CHROM'],
                                        variant['POS'],
                                        variant['REF'],
                                        alternatives[0]
                                    ]
                                )
        
        variant['vep_info'] = vep_dict
        variant['genetic_models'] = models_dict
        variant['compound_variants'] = compunds_dict
        variant['rank_scores'] = rank_score_dict
        variant['individual_scores'] = individual_score_dict
        
        return variant
    
    
    def __str__(self):
        """return the headers header lines to screen."""
        return '\n'.join(self.metadata.print_header())
        

@click.command()
@click.argument('variant_file',
        type=click.Path(),
        metavar='<vcf_file> or -'
)
@click.option('--vep', 
                    is_flag=True,
                    help='If variants are annotated with the Variant Effect Predictor.'
)
@click.option('-s' ,'--split', 
                    is_flag=True,
                    help='Split the variants with multiallelic calls.'
)
def cli(variant_file, vep, split):
    """Parses a vcf file.\n
        \n
        Usage:\n
            parser infile.vcf\n
        If pipe:\n
            parser - 
    """
    from datetime import datetime
    if variant_file == '-':
        my_parser = VCFParser(fsock=sys.stdin, split_variants=split)
    else:
        my_parser = VCFParser(infile = variant_file, split_variants=split)
    start = datetime.now()
    nr_of_variants = 0
    # for line in my_parser.metadata.print_header():
    #     print(line)
    for variant in my_parser:
        pp(variant)
        nr_of_variants += 1
    # print('Number of variants: %s' % nr_of_variants)
    # print('Time to parse: %s' % str(datetime.now()-start))
    # pp(my_parser.metadata.extra_info)
    

if __name__ == '__main__':
    cli()
