

def main():
    hotcheck()
    formatOut()

if __name__ == "__main__":
    from docopt import docopt
    usage = """
    Usage:
      snvchecker.py -v <variants> -r <reference> -o <output>  <bam>

    Options:
        -v,--variants=variants_list           variants list in chr,locs,ref-alle,alt-alle format
        -r,--reference=reference              reference genome 
        -o,--output=<variants_out>            variants output in tab-seperated file format
    """
    args = docopt(usage)
    print args

