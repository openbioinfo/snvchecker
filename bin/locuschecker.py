#!/usr/bin/env python
#coding=utf-8

import sys
import os

curdir = os.path.dirname(os.path.realpath(__file__))
modir = os.path.join(curdir,"../")
sys.path.insert(0,modir)
from snvchecker.hotcheck import hotcheck


def main(hotspot_bed, bamfile, method, prefix, outdir, reference):
    outfile = os.path.join(
                os.path.realpath(outdir),
                "%s.gt.xls" %prefix)
    GT = hotcheck.hotcheck(bamfile, hotspot_bed, method, prefix, outdir, reference)
    fh = open(outfile, "w+")
    for key, value in GT.iteritems():
        # print ":".join(map(str,key)) + "\t" + "\t".join(map(str,value))
        fh.write(":".join(map(str,key)) + "\t" + "\t".join(map(str,value))+"\n"  )
    fh.close()
    return 0

if __name__ == "__main__":
    from docopt import docopt
    usage = """
    Usage:
      snvchecker.py (freebayes -r <reference> | sambamba [-r <reference>]) -v <variants>  [-o <output>]  [-p prefix] <bam>

    Options:
        -v,--variants=snp_bed_file        bed format in chr,start, end, ref-allele, alt-allele format
        -r,--reference=<reference>        reference genome [default: /4_disk/genomes/hg19/hg19.fa]
        -p, --prefix=<prefix>             prefix for all outfile [default: 1]
        -o,--outdir=<outdir>              outdir for variants output: gt file in tab-seperated file format: chr,start, end, ref-allele, alt-allele, dp, alt-dp, gt, freq [default: 1]
    """
    args = docopt(usage)
    hotspot_bed  =   args['--variants']
    bamfile      =   args['<bam>']
    method       =   'freebayes' if args['freebayes'] else 'sambamba'
    prefix       =   args['--prefix']
    outdir       =   args['--outdir']
    reference    =   args['--reference']
    main(hotspot_bed, bamfile, method, prefix, )
    print args
