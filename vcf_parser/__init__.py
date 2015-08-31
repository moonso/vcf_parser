from __future__ import absolute_import

from pkg_resources import require

__version__ = require("vcf_parser")[0].version

from logging import getLogger
logger = getLogger(__name__)

from .header_parser import HeaderParser
from .log import init_log
from .genotype import Genotype
from .parser import VCFParser
