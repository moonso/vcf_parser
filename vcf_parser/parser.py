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

from __future__ import print_function
from __future__ import unicode_literals

import sys
import os
import gzip
import re
import argparse
from codecs import open, getreader
if sys.version_info < (2, 7):
    from ordereddict import OrderedDict
else:
    from collections import OrderedDict

from pprint import pprint as pp

from vcf_parser import genotype

class HeaderParser(object):
    """Parses a file with family info and creates a family object with individuals."""
    def __init__(self):
        super(HeaderParser, self).__init__()
        self.info_lines=[]
        self.info_dict=OrderedDict()
        
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
                raise SyntaxError("One of the INFO lines is malformed: %s" % line)
            matches = [match.group('id'), match.group('number'), match.group('type'), match.group('desc')]
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
                raise SyntaxError("One of the FILTER lines is malformed: %s" % line)
            matches = [match.group('id'), match.group('desc')]
            self.filter_lines.append(dict(list(zip(self.header_keys['filt'],matches))))
            self.filter_dict[match.group('id')] = line
        elif line_info[0] == 'contig':
            match = self.contig_pattern.match(line)
            if not match:
                raise SyntaxError("One of the contig lines is malformed: %s" % line)
            matches = [match.group('id'), match.group('length')]
            self.contig_lines.append(dict(list(zip(self.header_keys['contig'],matches))))
            self.contig_dict[match.group('id')] = line
        elif line_info[0] == 'FORMAT':
            match = self.format_pattern.match(line)
            if not match:
                raise SyntaxError("One of the FORMAT lines is malformed: %s" % line)
            matches = [match.group('id'), match.group('number'), match.group('type'), match.group('desc')]
            self.format_lines.append(dict(list(zip(self.header_keys['form'],matches))))
            self.format_dict[match.group('id')] = line
        elif line_info[0] == 'ALT':
            match = self.alt_pattern.match(line)
            if not match:
                raise SyntaxError("One of the ALT lines is malformed: %s" % line)
            matches = [match.group('id'), match.group('desc')]
            self.alt_lines.append(dict(list(zip(self.header_keys['alt'],matches))))
            self.alt_dict[match.group('id')] = line
        else:
            match = self.meta_pattern.match(line)
            if not match:
                raise SyntaxError("One of the meta data lines is malformed: %s" % line)
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


####            Parser:         ####


class VCFParser(object):
    """docstring for VCFParser"""
    def __init__(self, infile=None):
        super(VCFParser, self).__init__()
        self.infile = infile
        self.metadata = HeaderParser()
        self.individuals = []
        self.header = []
        file_name, file_extension = os.path.splitext(self.infile)
        if file_extension == '.gz':
            self.vcf = getreader('utf-8')(gzip.open(infile), errors='replace')
        elif file_extension == '.vcf':
            self.vcf = open(infile, mode='r', encoding='utf-8', errors='replace')
        else:
            raise SyntaxError('File is not in a supported format! Or use correct ending (.vcf or .vcf.gz)')
        
        self.next_line = self.vcf.readline().rstrip()
        self.current_line = self.next_line
        self.metadata.parse_meta_data(self.next_line)
        
        while self.next_line.startswith('#'):
            self.next_line = self.vcf.readline().rstrip()
            if self.next_line.startswith('##'):
                self.metadata.parse_meta_data(self.next_line)
            elif self.next_line.startswith('#'):
                self.metadata.parse_header_line(self.next_line)
        self.individuals = self.metadata.individuals
        self.header = self.metadata.header
            
    def __iter__(self):
        return self
    
    def next(self):
        """Works for python 2.x"""
        return self.__next__()
    
    def __next__(self):
        
        variant = dict(zip(self.header, self.next_line.split('\t')))
        self.next_line = self.vcf.readline().rstrip()
        
        if len(variant) < 8:
            raise StopIteration
        else:
            info_dict = {}
            vep_dict = {}
            genotype_dict = {}
            for info in variant.get('INFO', '').split(';'):
                info = info.split('=')
                if info[0] == 'CSQ':
                    for annotation in info[1].split(','):
                        vep_info = dict(zip(self.metadata.vep_columns, annotation.split('|')))
                        vep_dict[vep_info.get('SYMBOL','-')] = vep_info
                elif len(info) > 1:
                    info_dict[info[0]] = info[1]
                else:
                    info_dict[info[0]] = True
                    
            gt_format = variant.get('FORMAT', '').split(':')
            for individual in self.individuals:
                genotype_dict[individual] = genotype.Genotype(**dict(zip(gt_format, variant[individual].split(':'))))
                # ind_dict[individual] = dict(zip(gt_format, variant[individual].split(':')))
            
            variant['genotypes'] = genotype_dict            
            variant['info_dict'] = info_dict
            variant['variant_id'] = '_'.join([variant['CHROM'], variant['POS'], variant['REF'], variant['ALT']])
            variant['vep_info'] = vep_dict
            return variant
        
    
    def __str__(self):
        """return the headers header lines to screen."""
        return '\n'.join(self.metadata.print_header())
        

def main():
    from datetime import datetime
    parser = argparse.ArgumentParser(description="Parse vcf headers.")
    parser.add_argument('variant_file', 
                            type=str, 
                            nargs=1, 
                            help='A file with variant information.'
    )
    parser.add_argument('--vep', '-vep', 
                            action="store_true", 
                            help='If VEP info is provided.')
    args = parser.parse_args()
    infile = args.variant_file[0]
    my_parser = VCFParser(infile)
    start = datetime.now()
    # my_parser.parse()
    # print('parsing')
    # my_parser.metadataparser.add_info('ANN', '.', 'String', 'Annotates what feature(s) this variant belongs to.')
    # my_parser.metadataparser.add_info('Comp', '.', 'String', "':'-separated list of compound pairs for this variant.")
    # my_parser.metadataparser.add_info('GM', '.', 'String', "':'-separated list of genetic models for this variant.")
    # print(my_parser)
    nr_of_variants = 0
    for variant in my_parser:
        print('\t'.join([variant[head] for head in my_parser.header]).encode('utf-8'))
        if args.vep:
            pp(variant['vep_info'])
        nr_of_variants += 1
    print('Number of variants: %s' % nr_of_variants)
    print('Time to parse: %s' % str(datetime.now()-start))
    # print my_parser.__dict__
    # for line in my_parser.metadata:
    #     print line, my_parser.metadata[line]
    # print '\t'.join(my_parser.header)
    # print my_parser.line_counter
    # print my_parser.individuals
    # for individual in my_parser.individuals:
    #     for genotype in my_parser.individuals[individual]:
    #         print individual, genotype, my_parser.individuals[individual][genotype]
    


if __name__ == '__main__':
    main()
