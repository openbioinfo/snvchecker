#!/usr/bin/env python
#coding=utf-8

import sys
import os

curdir = os.path.dirname(os.path.realpath(__file__))
modir = os.path.join(curdir,"../")
sys.path.insert(0,modir)
from snvchecker.hotcheck import hotcheck

def main(hotspot_bed, bamfile):
    GT = hotcheck.hotcheck(bamfile, hotspot_bed)
    for key, value in GT.iteritems():
        print ":".join(map(str,key)) + "\t" + "\t".join(map(str,value))
    return 0
    formatOut()

if __name__ == "__main__":
    from docopt import docopt
    usage = """
    Usage:
      snvchecker.py -v <variants> [-r <reference>] [-o <output>]  <bam>

    Options:
        -v,--variants=snp_bed_file            bed format in chr,start, end, ref-allele, alt-allele format
        -r,--reference=reference              reference genome 
        -o,--output=<variants_out_file        variants output in tab-seperated file format: chr,start, end, ref-allele, alt-allele, dp, alt-dp, gt, freq
    """
    args = docopt(usage)

    hotspot_bed = args['--variants']
    bamfile     = args['<bam>']
    main(hotspot_bed, bamfile)
