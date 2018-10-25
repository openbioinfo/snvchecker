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
    fh.write("chrom\tstart\tend\tref\talt\tdp\tao\tfreq\tgt\n")
    for key, value in GT.iteritems():
        # print ":".join(map(str,key)) + "\t" + "\t".join(map(str,value))
        fh.write("\t".join(map(str,key)) + "\t" + "\t".join(map(str,value))+"\n"  )
    fh.close()
    return 0

if __name__ == "__main__":
    from docopt import docopt
    usage = """
    Usage:
      snvchecker.py (freebayes -r <reference> | sambamba [-r <reference>]) -v <variants>  [-o <output>]  [-p prefix] <bam>

    Options:
        -v,--variants=RiskInfoFile        risk infomation table (whith header)， this script will use following infomation: chrom, start, end, ref, alt. [default: {RiskInfoPath}]
        -r,--reference=<reference>        reference genome [default: /1_disk/ref/grch37.fa]
        -p, --prefix=<prefix>             prefix for all outfile [default: 1]
        -o,--outdir=<outdir>              outdir for variants output: gt file in tab-seperated file format: chr,start, end, ref-allele, alt-allele, dp, alt-dp, gt, freq [default: 1]
    """.format(RiskInfoPath = os.path.realpath(os.path.join(curdir,"../data", "RiskInfo_59site.tsv")) )
    args = docopt(usage)
    risk_info_table  =   args['--variants']
    bamfile          =   args['<bam>']
    method           =   'freebayes' if args['freebayes'] else 'sambamba'
    prefix           =   args['--prefix']
    outdir           =   args['--outdir']
    reference        =   args['--reference']
    main(risk_info_table, bamfile, method, prefix, outdir, reference)

