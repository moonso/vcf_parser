#!/usr/bin/env python
# encoding: utf-8
"""
vcf_parser

Command Line Interface for vcf_parser

Created by Måns Magnusson on 2014-12-12.
Copyright (c) 2014 __MoonsoInc__. All rights reserved.
"""

from __future__ import print_function

import sys
import os
import click

from pprint import pprint as pp
from datetime import datetime
from codecs import open

from vcf_parser import __version__, VCFParser


def print_version(ctx, param, value):
    """Callback function for printing version and exiting
    Args:
        ctx (object) : Current context
        param (object) : Click parameter(s)
        value (boolean) : Click parameter was supplied or not
    Returns:
        None:
    """
    if not value or ctx.resilient_parsing:
        return
    click.echo('vcf_parser version: ' + __version__)
    ctx.exit()


###         This is the main script         ###

@click.command()
@click.argument('variant_file', 
                    nargs=1, 
                    type=click.Path(),
                    metavar='<vcf_file> or -'
)
@click.option('-v', '--vep', 
                    is_flag=True,
                    help='If variants are annotated with the Variant Effect Predictor.'
)
@click.option('-s', '--split', 
                    is_flag=True,
                    help='Split the variants with multiallelic calls.'
)
@click.option('-o', '--outfile', 
                    type=click.Path(exists=False),
                    help='Path to a outfile.'
)
@click.option('-a',"--allele_symbol",
                default='0',
                help="The symbol that should be used when representing "\
                "unobserved alleles. Default is '0'"
)
@click.option('-v', '--verbose',
                is_flag=True,
                help='Increase output verbosity.'
)
@click.option('-i', '--check_info', 
                is_flag=True,
                help='Skip to check if info fields are correct for all variants.'
)
@click.option('--silent', 
                is_flag=True,
                help='Do not print vcf data.'
)
@click.option('--version',
                is_flag=True,
                callback=print_version,
                expose_value=False,
                is_eager=True
)
@click.option('-l', '--logfile',
                    type=click.Path(exists=False),
                    help="Path to log file. If none logging is "\
                          "printed to stderr."
)
@click.option('--loglevel',
                    type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR', 
                                        'CRITICAL']),
                    help="Set the level of log output."
)
def cli(variant_file, vep, split, outfile, verbose, silent, check_info,
        allele_symbol, logfile, loglevel):
    """
    Tool for parsing vcf files.
    
    Prints the vcf file to output. 
    If --split/-s is used all multiallelic calls will be splitted and printed 
    as single variant calls.
    For more information, please see github.com/moonso/vcf_parser.
    """
    from vcf_parser import logger, init_log

    if not loglevel:
        if verbose:
            loglevel = 'INFO'

    init_log(logger, logfile, loglevel)
    nr_of_variants = 0
    start = datetime.now()
    
    # with open(variant_file, 'r', encoding="utf-8") as f:
    #     for line in f:
    #         if not line.startswith('#'):
    #             nr_of_variants += 1
    
    if variant_file == '-':
        logger.info("Start parsing variants from stdin")
        my_parser = VCFParser(
            fsock=sys.stdin, 
            split_variants=split,
            check_info=check_info, 
            allele_symbol=allele_symbol
        )
    else:
        logger.info("Start parsing variants from file {0}".format(variant_file))
        my_parser = VCFParser(
            infile = variant_file,
            split_variants=split, 
            check_info=check_info, 
            allele_symbol=allele_symbol
        )

    if outfile:
        f = open(outfile, 'w', encoding='utf-8')
        logger.info("Printing vcf to file {0}".format(outfile))
    
    if not silent:
        logger.info("Printing vcf to stdout")
    else:
        logger.info("Skip printing since silent is active")
    
    for line in my_parser.metadata.print_header():
        if outfile:
            f.write(line+'\n')
        else:
            if not silent:
                print(line)
    try:
        for variant in my_parser:
            variant_line = '\t'.join([variant[head] for head in my_parser.header])
            if outfile:
                f.write(variant_line + '\n')
            else:
                if not silent:
                    print(variant_line)
            nr_of_variants += 1
    except SyntaxError as e:
        print(e)

    logger.info('Number of variants: {0}'.format(nr_of_variants))
    logger.info('Time to parse file: {0}'.format(str(datetime.now() - start)))
    # print('Number of variants: {0}'.format(nr_of_variants))
    # print('Time to parse file: {0}'.format(str(datetime.now() - start)))
    

if __name__ == '__main__':
    cli()
