from __future__ import division
import os

def parseDPfile(AlleleDP, hotspot_bed):
    hotspot_GT = dict() #key: tuple-[chr,pos], value: list-[dp, ao, ref, freq, gt]
    for line in open(hotspot_bed):
        lineArr = line.strip().split("\t")
        pos = tuple([lineArr[0], int(lineArr[2])])
        ref = lineArr[3]
        if not AlleleDP.has_key(pos):
            hotspot_GT[pos] = [0, 0, ref, "NA", "NA"]
            continue
        DP = sum(AlleleDP[pos])
        dp_base = dict(zip(['A','C','G','T'], AlleleDP[pos]))
        freq = dict(
               zip(['A','C','G','T'],  map(lambda x: x/DP, AlleleDP[pos]))
               )
        freq_ref = freq[ref]
        alt_freq = dict(filter(lambda x: x[0]!=ref, freq.iteritems()))

        alt_max = max(alt_freq, key=alt_freq.get)
        alt_freq_max = alt_freq[alt_max]
        #if pos == tuple([1,155157715]):
        #    print AlleleDP[pos]
        #    print freq
        #    print alt_freq
        if alt_freq_max < 0.1:
           hotspot_GT[pos] = [DP, 0, ref, 1-freq_ref, ref+ref]
        elif alt_freq_max < 0.8:
           hotspot_GT[pos] = [DP, dp_base[alt_max], ref, alt_freq_max, ref+alt_max]
        else:
           hotspot_GT[pos] = [DP, dp_base[alt_max], ref, alt_freq_max, alt_max+alt_max]
    return hotspot_GT

def getGenoTypeInfobySambamba(bam,hotspot_bed):
    os.system("sambamba depth base --nthreads 4 --min-coverage 1 --regions {} -o /tmp/out2.dp --min-base-quality 20  {}".format(hotspot_bed, bam))
    AlleleDP = dict() # key: tuple-[chr,pos], value: list-reads count for[A, C, G, T]
    for line in open("/tmp/out2.dp"):
        if line[0:3]=="REF": continue
        lineArr = line.strip().split("\t")
        pos = tuple([lineArr[0],int(lineArr[1])+1 ])
        AlleleDP[pos] = map(int, lineArr[3:7])
    hotspot_GT = parseDPfile(AlleleDP, hotspot_bed)
    #os.remove("/tmp/out2.dp")
    return hotspot_GT

def getGenoTypeInfobyFreebayes(bam,hotspot_bed):
    #os.system("freebayes --fasta-reference /1_disk/ref/grch37.fa --vcf /tmp/out.vcf --targets {} --min-base-quality 20  --min-alternate-count 1 -F 0.005 --use-duplicate-reads --pooled-continuous --report-monomorphic --no-mnps --no-complex {}".format(hotspot_bed, bam))
    os.system("freebayes --fasta-reference /4_disk/genomes/hg19/hg19.fa --vcf /tmp/out.vcf --targets {} --min-base-quality 20  --min-alternate-count 1 -F 0.005 --use-duplicate-reads --pooled-continuous --report-monomorphic --no-mnps --no-complex {}".format(hotspot_bed, bam))
    hotspot_GT = dict() #key: tuple-[chr,pos], value: list-[dp, ao, ref, freq, gt]
    for line in open("/tmp/out.vcf"):
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
        hotspot_GT[pos] = [dp, ao, ref, freq, gt]
    #os.remove("/tmp/out.vcf")
    return hotspot_GT

def hotcheck(bam,hotspot_bed):
    #hotspot_GT = getGenoTypeInfobySambamba(bam,hotspot_bed)
    hotspot_GT = getGenoTypeInfobyFreebayes(bam,hotspot_bed)
    return hotspot_GT

