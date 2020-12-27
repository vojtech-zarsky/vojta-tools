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
        if hasattr(file, 'read'):
            self.name = file.name
            self.file = file
            self.byte = False
        elif file[-3:] == '.gz':
            self.name = file
            import indexed_gzip
            self.file = indexed_gzip.IndexedGzipFile(self.name)
            self.byte = True
        else:
            self.name = file
            self.file = open(self.name)
            self.byte = False
        self.end = 1
        sLast = ''
        for sLine in self.file:
            if self.byte:
                sLine = sLine.decode()
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
                sLine = next(self.file)
                if self.byte:
                    sLine = sLine.decode()
                sLine = sLine.strip('\n')
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
        if self.byte:
            import indexed_gzip
            handle = indexed_gzip.IndexedGzipFile(self.name)
        else:
            handle = open(self.name)
        fileOut = open(self.name+'.index','w')
        iIndex = handle.tell()
        sLine = handle.readline()
        if self.byte:
            sLine = sLine.decode()
        if sDb == 'uniprot':
            while sLine:
                if len(sLine) > 1 and sLine[0] == '>':
                    fileOut.write(sLine.split()[0].split('|')[2]+'\t'+str(iIndex)+'\n')
                iIndex = handle.tell()
                sLine = handle.readline()
                if self.byte:
                    sLine = sLine.decode()
        else:
            while sLine:
                if len(sLine) > 1 and sLine[0] == '>':
                    fileOut.write(sLine.split()[0][1:]+'\t'+str(iIndex)+'\n')
                iIndex = handle.tell()
                sLine = handle.readline()
                if self.byte:
                    sLine = sLine.decode()
        return 0
    def readIndex(self):
        self.dId2Index = {}
        for sLine in open(self.name+'.index'):
            (sId, sIndex) = sLine.strip().split('\t',1)
            self.dId2Index[sId] = int(sIndex)
        return 0
    def getSeq(self, sId, iIndex=None):
        if iIndex == None:
            iIndex = self.dId2Index[sId]
        self.file.seek( iIndex )
        sLine = self.file.readline()
        if self.byte:
            sLine = sLine.decode()
        sHeader = sLine[1:].strip()
        sSeq = ''
        for sLine in self.file:
            if self.byte:
                sLine = sLine.decode()
            sLine = sLine.strip()
            if sLine == '' or sLine[0] == '>':
                break
            sSeq += sLine
        return Seq(sHeader.split()[0], sHeader, sSeq)

"""
import indexed_gzip
#f = indexed_gzip.IndexedGzipFile('../databases/uniprot/uniprot_sprot.201002.fasta.gz')
#print(dir(f))
#print(f.__str__)
seqio = SeqIO('temp.fasta')
#seqio = SeqIO(open('temp.fasta'))
for seq in seqio:
    print(seq.getString())
    break
seqio.writeIndex()
seqio.readIndex()
seq = seqio.getSeq('sp|Q6GZX3|002L_FRG3G')
print(seq.getString())
#"""

