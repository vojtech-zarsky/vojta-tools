import sys
import os
import re
import numpy as np
import matplotlib.pyplot as plt
import argparse

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
        print('seeking within file...')
        self.file.seek( iIndex )
        print('done seeking.')
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

def partitionAlnByGaps(sAlnFile, iPercentStep=1, iMaxPercentGaps=50):
    iPercentStep = int(iPercentStep)
    iMaxPercentGaps = int(iMaxPercentGaps)
    lSeqIds = []
    lColumn2Count = None
    for seq in SeqIO(sAlnFile):
        lSeqIds.append(seq.id)
        if lColumn2Count == None:
            lColumn2Count = [0]*len(seq.seq)
        for iPos in range(len(seq.seq)):
            if seq.seq[iPos].isalpha():
                lColumn2Count[iPos] += 1
    aColumn2Count = np.array(lColumn2Count)
    print('seqs:', len(lSeqIds), 'columns:', aColumn2Count.size)
    aColumn2Count = aColumn2Count/len(lSeqIds)
    tableOut = open('{}.partitionAlnByGaps.tsv'.format(sAlnFile.split('/')[-1]), 'w')
    tableOut.write('percent gaps cutoff\tcolumns added\tsum columns\tseqs kept\n')
    print('percent gaps cutoff\tcolumns added\tsum columns\tseqs kept')
    iLastColumnsToKeepLen = -1
    iSeqsKept = 0
    lColumnsToKeepPrint = []
    lSeqsToKeepPrint = []
    lPercentCutoffs = range(0, 100+iPercentStep, iPercentStep)[::-1] ##
    lColors = ['#917c6f', '#3771c8', '#ff8080'] ##
    for iPercentCutoff in lPercentCutoffs:
        lColumnsToKeep = []
        iCountColumnsInTheStep = 0
        for iPos in range(aColumn2Count.size):
            if aColumn2Count[iPos]*100 >= iPercentCutoff:
                lColumnsToKeep.append(iPos)
                if aColumn2Count[iPos]*100 < iPercentCutoff+iPercentStep:
                    iCountColumnsInTheStep += 1
        if len(lColumnsToKeep) > 0 and len(lColumnsToKeep) != iLastColumnsToKeepLen:
            iLastColumnsToKeepLen = len(lColumnsToKeep)
            lSeqsOut = []
            for seq in SeqIO(sAlnFile):
                sSeq = ''
                for iPosToKeep in lColumnsToKeep:
                    sSeq += seq.seq[iPosToKeep]
                if len(re.sub('[^a-zA-Z]','',sSeq)) >= len(sSeq)*iMaxPercentGaps/100:
                    sHeader = seq.id
                    if len(seq.desc.strip()) > 0:
                        sHeader += ' '+seq.desc.strip()
                    lSeqsOut.append('>{}\n{}\n'.format(sHeader, sSeq))
            iSeqsKept = len(lSeqsOut)
            if len(lSeqsOut) > 0:
                fastaOut = open('{}.partitionAlnByGaps.{}.fasta'.format(sAlnFile.split('/')[-1], iPercentCutoff), 'w')
                for sSeqOut in lSeqsOut:
                    fastaOut.write(sSeqOut)
        tableOut.write('{}\t{}\t{}\t{}\n'.format(iPercentCutoff, iCountColumnsInTheStep, len(lColumnsToKeep), iSeqsKept))
        print('{}\t{}\t{}\t{}'.format(iPercentCutoff, iCountColumnsInTheStep, len(lColumnsToKeep), iSeqsKept))
        plt.bar(iPercentCutoff, iCountColumnsInTheStep*100/len(aColumn2Count), color=lColors[0])
        lColumnsToKeepPrint.append( len(lColumnsToKeep)*100/len(aColumn2Count) )
        lSeqsToKeepPrint.append( iSeqsKept*100/len(lSeqIds) )
    plt.plot(lPercentCutoffs, lColumnsToKeepPrint, linestyle='-', linewidth=1,  marker='o', color=lColors[1])
    plt.plot(lPercentCutoffs, lSeqsToKeepPrint, linestyle='-', linewidth=1, marker='+', color=lColors[2])
    plt.text(80, 90, 'Columns added', color=lColors[0])
    plt.text(80, 80, 'o Columns', color=lColors[1])
    plt.text(80, 70, '+ Sequences', color=lColors[2])
    plt.xlabel('percent gaps in column cutoff')
    plt.ylabel('percent')
    plt.savefig('{}.partitionAlnByGaps.pdf'.format(sAlnFile.split('/')[-1]))
    return 0

argumentParser = argparse.ArgumentParser()
argumentParser.add_argument("--percent_step", default=5, help = "percent step (default: %(default)s)")
argumentParser.add_argument("--max_percent_gaps", help="max percent gaps allowed in a sequence (default: %(default)s)", default=50)
argumentParser.add_argument("alignment", help="alignment file")
args = argumentParser.parse_args()
partitionAlnByGaps(args.alignment, iPercentStep=args.percent_step, iMaxPercentGaps=args.max_percent_gaps)

