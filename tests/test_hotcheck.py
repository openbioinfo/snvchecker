import sys
sys.path.append("../")
from snvchecker.hotcheck.hotcheck import hotcheck

bam = "data/test.sorted.bam"
riskinfo = "data/test_RiskInfo.tsv"

def test_hotcheck():
    hotcheck(bam,riskinfo,method="sambamba")


if __name__ == "__main__":
    test_hotcheck()



