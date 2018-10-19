import sys
sys.path.append("../")
from snvchecker.hotcheck.hotcheck import hotcheck

bam = "data/T3.sorted.bam"
hot = "data/yigan.bed"

def test_hotcheck():
    hotcheck(bam,hot)


if __name__ == "__main__":
    test_hotcheck()



