#!/usr/bin/env python
#coding=utf-8

import sys
import os

sys.path.insert(0,os.path.dirname(__file__) + "../")
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
        -v,--variants=variants_list           variants list in chr,locs,ref-alle,alt-alle format
        -r,--reference=reference              reference genome 
        -o,--output=<variants_out>            variants output in tab-seperated file format
    """
    args = docopt(usage)

    hotspot_bed = args['--variants']
    bamfile     = args['<bam>']
    main(hotspot_bed, bamfile)
