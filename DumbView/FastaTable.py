from pbcore.io.base import *
from collections import namedtuple, OrderedDict
import mmap, numpy as np

FaiRecord = namedtuple("FaiRecord", ("name", "length", "offset", "lineWidth", "stride"))

def faiFilename(fastaFilename):
    return fastaFilename + ".fai"

def loadFastaIndex(filename):
    tbl = OrderedDict()
    for line in open(filename):
        name, length_, offset_, lineWidth_, blen_ = line.split()
        record = FaiRecord(name, int(length_), int(offset_),
                           int(lineWidth_), int(blen_))
        tbl[name] = record
    return tbl

def fileOffset(fai, name, pos):
    """
    Find the in-file position (in bytes) corresponding to the position
    in the named contig, using the FASTA index.
    """
    faiRecord = fai[name]
    assert 0 <= pos < faiRecord.length
    q, r = divmod(pos, faiRecord.lineWidth)
    offset = faiRecord.offset + q*faiRecord.stride + r
    return offset

class FastaTable(ReaderBase):
    def __init__(self, f):
        self.file = open(f, "r")
        self.faiFilename = faiFilename(f)
        self.fai = loadFastaIndex(self.faiFilename)

    def sequence(self, name, start, end):
        assert name in self.fai
        assert start < end
        startOffset = fileOffset(self.fai, name, start)
        endOffset   = fileOffset(self.fai, name, end)
        self.file.seek(startOffset)
        snip = self.file.read(endOffset-startOffset).replace("\n", "")
        return snip

    def length(self, name):
        assert name in self.fai
        return self.fai[name].length
