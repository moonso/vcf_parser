from vcf_parser import HeaderParser
import pytest

def test_parse_vcf_lines():
    """
    Test how the header parser behaves with simple vcf lines
    """
    
    header_parser = HeaderParser()
    
    header_lines = [
        '##fileformat=VCFv4.2',
        '##FILTER=<ID=LowQual,Description="Low quality">',
        '##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">',
        '##INFO=<ID=CNT,Number=A,Type=Integer,Description="Number of times '\
        'this allele was found in external db">',
        '##contig=<ID=1,length=249250621,assembly=b37>',
        '##INFO=<ID=DP_HIST,Number=R,Type=String,Description="Histogram for '\
        'DP; Mids: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|'\
        '62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5">',
        '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for'\
        ' the ref and alt alleles in the order listed">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=GQ,Number=1,Type=String,Description="GenotypeQuality">'
        '##reference=file:///human_g1k_v37.fasta',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfather\tmother\tproband'
    ]
    for line in header_lines:
        if line.startswith('##'):
            header_parser.parse_meta_data(line)
        elif line.startswith('#'):
            header_parser.parse_header_line(line)
    
    assert header_parser.fileformat == "VCFv4.2"
    assert header_parser.individuals == ['father','mother','proband']
    
    assert header_parser.vep_columns == []
    
    assert "MQ" in header_parser.extra_info
    assert header_parser.extra_info["MQ"]['Description'] == "RMS Mapping Quality"
    assert header_parser.extra_info["CNT"]['Number'] == "A"
    assert header_parser.extra_info["CNT"]['Type'] == "Integer"
    assert "CNT" in header_parser.extra_info
    assert "DP_HIST" in header_parser.extra_info
    
    assert "LowQual" in header_parser.filter_dict
    assert "1" in header_parser.contig_dict
    
    assert header_parser.header == [
        'CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',
        'father','mother','proband'
    ]

def test_malformed_lines():
    """
    Test how the header parser behaves with simple vcf lines
    """
    
    header_parser = HeaderParser()
    
    malformed_fileformat = '##fileformat'
    malformed_info_line = '##INFO=<ID=MQ,Number=1,Description="RMS Mapping Quality">'
    malformed_contig_line = '##contig=<ID=1,assembly=b37>'
    
    with pytest.raises(SyntaxError):
        header_parser.parse_meta_data(malformed_fileformat)

    with pytest.raises(SyntaxError):
        header_parser.parse_meta_data(malformed_info_line)

    with pytest.raises(SyntaxError):
        header_parser.parse_meta_data(malformed_contig_line)

def test_vep_columns():
    """
    Test how the vep columns are parsed
    """
    header_parser = HeaderParser()
    
    vep_info_line = '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence'\
    ' type as predicted by VEP. Format: Allele|Gene|Feature|Feature_type|Consequence">'
    
    header_parser.parse_meta_data(vep_info_line)
    
    assert header_parser.vep_columns == ['Allele','Gene','Feature','Feature_type','Consequence']
    
    


    