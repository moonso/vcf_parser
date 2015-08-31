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

import sys
import os
import gzip
import re
import pkg_resources
import click
import locale
import logging


from codecs import open, getreader


from vcf_parser import (Genotype, HeaderParser)
from vcf_parser.utils import (format_variant, split_variants)

####            Parser:         ####

class VCFParser(object):
    """docstring for VCFParser"""
    def __init__(self, infile=None, fsock=None, split_variants=False, 
                check_info=False, allele_symbol='0', fileformat = None):
        super(VCFParser, self).__init__()
        self.logger = logging.getLogger(__name__)
        
        self.vcf = None
        self.logger.debug("Set self.vcf to:{0}".format(self.vcf))
        self.beginning = True
        self.infile = infile
        self.fsock = fsock
        self.split_variants = split_variants
        self.logger.info("Split variants = {0}".format(self.split_variants))
        self.fileformat = fileformat
        
        self.check_info = check_info
        self.logger.info("check info = {0}".format(self.check_info))

        self.allele_symbol = allele_symbol
        self.logger.info("Allele symbol = {0}".format(self.allele_symbol))
        
        self.logger.info("Initializing HeaderParser")
        self.metadata = HeaderParser()
        # These are the individuals described in the header
        self.individuals = []
        # This is the header line of the vcf
        self.header = []
        
        # If there are no file or stream the user can add variants manually.
        # These will be added to self.variants
        self.variants = []
        
        if (fsock or infile):
        
            if fsock:
                if not infile and hasattr(fsock, 'name'):
                    self.logger.info("Reading vcf form stdin")
                    if sys.version_info < (3, 0):
                        self.logger.info("Using codecs to read stdin")
                        sys.stdin = getreader('utf-8')(fsock)
                    
                    self.vcf = sys.stdin
            
            else:
                self.logger.info("Reading vcf form file {0}".format(infile))
                file_name, file_extension = os.path.splitext(infile)
                if file_extension == '.gz':
                    self.logger.debug("Vcf is zipped")
                    self.vcf = getreader('utf-8')(gzip.open(infile), errors='replace')
                elif file_extension == '.vcf':
                    self.vcf = open(infile, mode='r', encoding='utf-8', errors='replace')
                else:
                    raise IOError("File is not in a supported format!\n"
                                        " Or use correct ending(.vcf or .vcf.gz)")
            
            self.logger.debug("Reading first line.")
            self.next_line = self.vcf.readline().rstrip()
            self.current_line = self.next_line
           
            # First line is allways a metadata line
            if not self.next_line.startswith('#'):
                raise IOError("VCF files allways have to start with a metadata line.")
            self.metadata.parse_meta_data(self.next_line)
            
            # Parse the metadata lines
            while self.next_line.startswith('#'):
                if self.next_line.startswith('##'):
                    self.metadata.parse_meta_data(self.next_line)
                elif self.next_line.startswith('#'):
                    self.metadata.parse_header_line(self.next_line)
                self.next_line = self.vcf.readline().rstrip()
            
            self.individuals = self.metadata.individuals
            self.logger.info("Setting self.individuals to {0}".format(
                self.individuals
            ))
            self.header = self.metadata.header
            self.vep_header = self.metadata.vep_columns
        
        else:
            if not self.fileformat:
                raise IOError("Please initialize with a fileformat.")
            else:
                self.metadata.fileformat = self.fileformat
    
    def add_variant(self, chrom, pos, rs_id, ref, alt, qual, filt, info, form=None, genotypes=[]):
        """
        Add a variant to the parser.
        
        This function is for building a vcf. It takes the relevant parameters 
        and make a vcf variant in the proper format.
        """
        variant_info = [chrom, pos, rs_id, ref, alt, qual, filt, info]
        if form:
            variant_info.append(form)
        for individual in genotypes:
            variant_info.append(individual)
        
        variant_line = '\t'.join(variant_info)
        variant = format_variant(
            line = variant_line, 
            header_parser = self.metadata, 
            check_info = self.check_info
        )
        
        if not (self.split_variants and len(variant['ALT'].split(',')) > 1):
            self.variants.append(variant)
            
        # If multiple alternative and split_variants we must split the variant                 
        else:
            for splitted_variant in split_variants(
                                                    variant_dict=variant, 
                                                    header_parser=self.metadata, 
                                                    allele_symbol=self.allele_symbol):
                self.variants.append(splitted_variant)
    
    def __iter__(self):
        
        if not self.metadata.fileformat:
            raise SyntaxError("Vcf must have fileformat defined")
        
        if self.vcf:
            
            # We need to treat the first case as an exception
            if self.beginning:
                variants = []
                first_variant = format_variant(
                    line = self.next_line, 
                    header_parser = self.metadata, 
                    check_info = self.check_info
                )
                
                if not (self.split_variants and len(first_variant['ALT'].split(',')) > 1):
                    variants.append(first_variant)
                else:
                    for splitted_variant in split_variants(
                                                            variant_dict=first_variant, 
                                                            header_parser=self.metadata, 
                                                            allele_symbol=self.allele_symbol):
                        variants.append(splitted_variant)

                
                for variant in variants:
                    yield variant
                
                self.beginning = False
                
            
            for line in self.vcf:
                line = line.rstrip()
                # These are the variant(s) found in one line of the vcf
                # If there are multiple alternatives and self.split_variants
                # There can be more than one variant in one line
                variants = []
                
                if not line.startswith('#') and len(line.split('\t')) >= 8:
                    variant = format_variant(
                        line = line, 
                        header_parser = self.metadata, 
                        check_info = self.check_info
                        )
                    
                    if not (self.split_variants and len(variant['ALT'].split(',')) > 1):
                        variants.append(variant)
                    
                    else:
                        for splitted_variant in split_variants(
                                    variant_dict=variant, 
                                    header_parser=self.metadata, 
                                    allele_symbol=self.allele_symbol):
                            variants.append(splitted_variant)
                
                for variant in variants:
                    yield variant
        
        else:
            for variant in self.variants:
                yield variant

    def __repr__(self):
        return "Parser(infile={0},fsock={1},split_variants={2})".format(
            self.infile, self.fsock, self.split_variants
        )

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
    from pprint import pprint as pp
    if variant_file == '-':
        my_parser = VCFParser(fsock=sys.stdin, split_variants=split)
    else:
        my_parser = VCFParser(infile = variant_file, split_variants=split)
    start = datetime.now()
    nr_of_variants = 0
    for line in my_parser.metadata.print_header():
        print(line)
    for variant in my_parser:
        pp(variant)
        nr_of_variants += 1
    print('Number of variants: %s' % nr_of_variants)
    # print('Time to parse: %s' % str(datetime.now()-start))
    # pp(my_parser.metadata.extra_info)
    

if __name__ == '__main__':
    cli()
