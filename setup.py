try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup
    
# For making things look nice on pypi:
# try:
#     import pypandoc
#     long_description = pypandoc.convert('README.md', 'rst')
# except (IOError, ImportError):

long_description = 'Tool for parsing Variant Call Format (VCF) files. Works like a lightweight version of PyVCF.'

setup(name='vcf_parser',
    version='1.5.8',
    description='Parsing vcf files',
    author = 'Mans Magnusson',
    author_email = 'mans.magnusson@scilifelab.se',
    url = 'http://github.com/moonso/vcf_parser',
    license = 'MIT License',
    install_requires=[
        'pytest', 
        'click'
    ],
    packages = [
        'vcf_parser',
        'vcf_parser.utils',
        'vcf_parser.cli',
    ],
    keywords = [
        'parser', 
        'vcf', 
        'variants'
    ],
    entry_points = {
        'console_scripts': [
            'vcf_parser = vcf_parser.cli.command_line:cli'
        ]
    },
    classifiers = [
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Development Status :: 4 - Beta",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Unix",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    long_description = long_description,
)