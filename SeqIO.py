## for python3

import re

class Seq:
    def __init__(self, sId, sDesc, sSeq):
        self.id = sId
        self.desc = ''
        if len(re.split('\s+', sDesc, 1)) > 1:
            self.desc = re.split('\s+', sDesc, 1)[1]
        self.seq = sSeq
    def getString(self):
        return '>'+self.id+' '+self.desc+'\n'+self.seq+'\n'

class SeqIO:
    def __init__(self, file):
        self.file = file
        self.end = 1
        sLast = ''
        for sLine in self.file:
            sLine = sLine.strip('\n')
            if re.search('^>', sLast) and len(sLine) > 0 and sLine[0] != '\s':
                self.header = sLast.lstrip('>').strip('\n')
                self.seq = sLine.strip()
                self.end = 0
                break
            sLast = sLine
    def __iter__(self):
        return self
    def __next__(self):
        return self.next()
    def end(self):
        return self.end
    def next(self):
        if self.end:
            raise StopIteration()
        sLine = ''
        while True:
            try:
                sLine = next(self.file).strip('\n')
                if re.search('^>',sLine):
                    seq = Seq(self.header.split()[0], self.header, self.seq)
                    self.header = sLine.lstrip('>').strip('\n')
                    self.seq = ''
                    return seq
                self.seq = self.seq+sLine.strip()
            except StopIteration:
                seq = Seq(self.header.split()[0], self.header, self.seq)
                self.header = None
                self.seq = None
                self.end = 1
                return seq
                break
        return 0
    def getDict(self):
        dSeqs = {}
        while not self.end:
            seq = self.next()
            dSeqs[seq.id] = seq
        return dSeqs
    def writeIndex(self, sDb=None):
        handle = open(self.file.name)
        fileOut = open(handle.name+'.index','w')
        iIndex = handle.tell()
        sLine = handle.readline()
        if sDb == 'uniprot':
            while sLine:
                if len(sLine) > 1 and sLine[0] == '>':
                    fileOut.write(sLine.split()[0].split('|')[2]+'\t'+str(iIndex)+'\n')
                iIndex = handle.tell()
                sLine = handle.readline()
        else:
            while sLine:
                if len(sLine) > 1 and sLine[0] == '>':
                    fileOut.write(sLine.split()[0][1:]+'\t'+str(iIndex)+'\n')
                iIndex = handle.tell()
                sLine = handle.readline()
        return 0
    def readIndex(self):
        self.dId2Index = {}
        for sLine in open(self.file.name+'.index'):
            (sId, sIndex) = sLine.strip().split('\t',1)
            self.dId2Index[sId] = int(sIndex)
        return 0
    def getSeq(self, sId, iIndex=None):
        if iIndex == None:
            iIndex = self.dId2Index[sId]
        self.file.seek( iIndex )
        sHeader = self.file.readline()[1:].strip()
        sSeq = ''
        for sLine in self.file:
            sLine = sLine.strip()
            if sLine == '' or sLine[0] == '>':
                break
            sSeq += sLine
        return Seq(sHeader.split()[0], sHeader, sSeq)

"""
seqio = SeqIO('test.fasta')
seqio.writeIndex()
seqio.readIndex()
seq = seqio.getSeq('XP_012898703.1')
print seq.id, seq.seq
"""

