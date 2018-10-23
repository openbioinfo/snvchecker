from __future__ import division
import os

JudgeGTFreqCutoff = {
   'WildType' : 0.1,
   'HomozygousMutant' : 0.8
}


# parse info from sambamba output, and convert to need data(dp, freq, gt) with uniform format
# input:
#      AlleleDP:    dict     : { tuple[chr,pos] => list[A_dp, C_dp, G_dp, T_dp]  }
#      hotspot_bed: bed file : chrom start end ref_allele alt_allele
# output:
#     hotspot_GT:  dict      : { tuple[chr,pos] => list[dp, ao, ref, freq, gt] }
def parseDPfile(AlleleDP, hotspot_bed):
    hotspot_GT = dict() #key: tuple-[chr,pos], value: list-[dp, ao, ref, freq, gt]

    # read hotspot bedfile to get ref_allele
    for line in open(hotspot_bed):
        lineArr = line.strip().split("\t")
        pos = tuple([lineArr[0], int(lineArr[2])])
        ref = lineArr[3]
        if not AlleleDP.has_key(pos):
            hotspot_GT[pos] = [0, 0, ref, "NA", "NA"]
            continue
        DP = sum(AlleleDP[pos])
        # make a dict to store dp for each base
        dp_base = dict(zip(['A','C','G','T'], AlleleDP[pos]))
        # make a dict to store freq for each base
        freq = dict(
               zip(['A','C','G','T'],  map(lambda x: x/DP, AlleleDP[pos]))
               )
        # get ref_allele freq
        freq_ref = freq[ref]

        # get alt_base and alt_freq
        alt_freq = dict(filter(lambda x: x[0]!=ref, freq.iteritems()))
        alt_max = max(alt_freq, key=alt_freq.get)
        alt_freq_max = alt_freq[alt_max]

        # judge whether Wildtype or Heterozygote or Homozygote
        key = tuple([lineArr[0], int(lineArr[2]), lineArr[3], lineArr[4]])
        if alt_freq_max   <= JudgeGTFreqCutoff['WildType']:
           hotspot_GT[key] = [DP, 0, ref, 1-freq_ref, ref+ref]
        elif alt_freq_max >= JudgeGTFreqCutoff['HomozygousMutant']:
           hotspot_GT[key] = [DP, dp_base[alt_max], ref, alt_freq_max, alt_max+alt_max]
        else:
           hotspot_GT[key] = [DP, dp_base[alt_max], ref, alt_freq_max, ref+alt_max]
    return hotspot_GT


# get genoType by Sambamba
# input:
#     bam: indexed bamfile
#     hotspot_bed: bed file : chrom start end ref_allele [alt_allele]
#     prefix: prefix for sambamba outfile
#     outdir: outdir for sambamba outfile
# output: hotspot_GT:   dict      : { tuple[chr,pos] => list[dp, ao, ref, freq, gt] }
def getGenoTypeInfobySambamba(bam, hotspot_bed, prefix, outdir):
    sambamba_out = "%s/%s.sambamba.dp" %(outdir, prefix)
    os.system("sambamba depth base --nthreads 4 --min-coverage 1 --regions {region} -o {outfile} --min-base-quality 20  {bam}".format(
                                              region = hotspot_bed, 
                                              outfile = sambamba_out,
                                              bam = bam))

    AlleleDP = dict() # key: tuple-[chr,pos], value: list-reads count for[A, C, G, T]
    for line in open(sambamba_out):
        if line[0:3]=="REF": continue
        lineArr = line.strip().split("\t")
        pos = tuple([lineArr[0],int(lineArr[1])+1 ])
        AlleleDP[pos] = map(int, lineArr[3:7])

    hotspot_GT = parseDPfile(AlleleDP, hotspot_bed)
    #if outdir == "/tmp": os.remove(sambamba_out)
    return hotspot_GT


# get genoType by Freebayes
# input:
#     bam: indexed bamfile
#     hotspot_bed: bed file : chrom start end ref_allele [alt_allele]
#     prefix: prefix for Freebayes vcf outfile
#     outdir: outdir for Freebayes vcf outfile
#     reference: reference fasta file
# output: hotspot_GT:   dict      : { tuple[chr,pos] => list[dp, ao, ref, freq, gt] }
def getGenoTypeInfobyFreebayes(bam, hotspot_bed, prefix, outdir, reference):
    Freebayes_out = "%s/%s.freebayes.vcf" %(outdir, prefix)
    os.system("freebayes --fasta-reference {reference} --vcf {outfile} --targets {region} --min-base-quality 20  --min-alternate-count 1 -F 0.005 --use-duplicate-reads --pooled-continuous --report-monomorphic --no-mnps --no-complex {bam}".format(
                                            reference = reference,
                                            outfile = Freebayes_out,
                                            region = hotspot_bed,
                                            bam = bam))

    hotspot_GT = dict() #key: tuple-[chr,pos], value: list-[dp, ao, ref, freq, gt]
    for line in open(Freebayes_out):
        if line[0] == "#": continue
        lineArr = line.strip().split("\t")
        pos = tuple([lineArr[0], int(lineArr[1])])
        ref = lineArr[3]
        alt = lineArr[4]
        info_arr = lineArr[-1].split(":") #GT:DP:RO:QR:AO:QA:GL -> 1/1:28988:4:144:28982:1010518:-90893,-8712.61,0
        gt_type = info_arr[0]
        dp = int(info_arr[1])
        ro = int(info_arr[2])
        ao = (dp - ro) if info_arr[4] == "." else int(info_arr[4])
        freq = ao/dp
        if gt_type == "0/1" or gt_type == "1/0":
            gt = ref+alt
        elif gt_type == "1/1":
            gt = alt+alt
        elif gt_type == "1/2" or gt_type == "2/1":
            gt = "".join(alt.split(","))
        else:
            gt = ref + ref
        key = tuple([lineArr[0], int(lineArr[1]), lineArr[2], lineArr[3]])
        hotspot_GT[pos] = [dp, ao, ref, freq, gt]
    #if outdir == "/tmp": os.remove(sambamba_out)
    return hotspot_GT


# get genoType by Freebayes/Sambamba
# input:
#     bam: indexed bamfile
#     hotspot_bed: bed file : chrom start end ref_allele [alt_allele]
#     method: call genotype method: string: "sambamba" or "freebayes" as so far
#     prefix: prefix for Freebayes vcf outfile: default = 1(auto-detect, first part cut by "."),otherwise use input
#     outdir: outdir for Freebayes vcf outfile: default = 1(auto-detect, same dir of bamfile), if NA or empty, then /tmp, otherwise use input
#     reference: reference fasta file: default = /4_disk/genomes/hg19/hg19.fa
# output: hotspot_GT:   dict      : { tuple[chr,pos] => list[dp, ao, ref, freq, gt] }
def hotcheck(bam,hotspot_bed, method = "sambamba", prefix = 1, outdir = 1, reference = "/4_disk/genomes/hg19/hg19.fa"):
    if prefix == 1:
        prefix = os.path.basename(bam).split(".")[0]
    if outdir == 1:
        outdir = os.path.dirname(os.path.realpath(bam))
    elif outdir == "" or outdir == "NA":
        outdir = "/tmp"

    if method == "sambamba":
        hotspot_GT = getGenoTypeInfobySambamba(bam, hotspot_bed, prefix, outdir)
    elif method == "freebayes":
        hotspot_GT = getGenoTypeInfobyFreebayes(bam, hotspot_bed, prefix, outdir, reference)
    else:
         hotspot_GT = getGenoTypeInfobySambamba(bam, hotspot_bed, prefix, outdir) # default sambamba method
    return hotspot_GT

