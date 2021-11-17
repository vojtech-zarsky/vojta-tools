
# python ../myclasses/pepaln2dna.py --nt GTA_capsid_nucleotides_Rickettsiales.fasta --pep_aln GTA_capsid_Rickettsiales_amino_acids_aligned.fasta --pep_aln_trimmed GTA_capsid_Rickettsiales_amino_acids_aligned.bmge.fasta

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

import os
import argparse

parser = argparse.ArgumentParser(prog='Build nucleotide sequence alignment based on a protein alignment. The nucleotide sequence id or description must contain the protein id, or the corresponding sequences must have the same order (specify the --order argument). Creates a .aln.fna file. (also .aln.trimmed.fna if the trimmed protein alignment is specified)')

parser.add_argument('--nt', required=True, help='Coding DNA matching the proteins in the .pep file.')
parser.add_argument('--pep_aln', required=True, help='Full-length protein alignment.')
parser.add_argument('--pep_aln_trimmed', required=False, help='Trimmed protein alignment.')
parser.add_argument('--order', action='store_true', required=False, help='Sequences have the same order in the nt and pep files.')
args = parser.parse_args()

def aln2columns(seqs):
    lAlnIds = []
    lAlnColumns = None
    if isinstance(seqs, str):
        seqs = SeqIO(seqs)
    for seq in seqs:
        seq.seq = seq.seq.upper() ## to upper case
        lAlnIds.append(seq.id)
        if lAlnColumns == None:
            lAlnColumns = []
            for iPos in range(len(seq.seq)):
                lAlnColumns.append([])
        for iPos in range(len(seq.seq)):
            lAlnColumns[iPos].append(seq.seq[iPos])
    return lAlnIds, lAlnColumns

def findTrimmedAlnPositions(sFullAlnFile, sTrimmedAlnFile):
    lTrimmedAlnIds, lTrimmedAlnColumns = aln2columns(sTrimmedAlnFile)
    lFullAlnSequences = []
    dFullAlnSequences = SeqIO(sFullAlnFile).getDict()
    for sTrimmedAlnId in lTrimmedAlnIds:
        bMatch = False
        ## look for exact match
        if sTrimmedAlnId in dFullAlnSequences:
            lFullAlnSequences.append(dFullAlnSequences[sTrimmedAlnId])
            bMatch = True
        if not bMatch:
            for sFullAlnId in dFullAlnSequences.keys():
                if sFullAlnId in sTrimmedAlnId:
                    lFullAlnSequences.append(dFullAlnSequences[sFullAlnId])
                    bMatch = True
        if not bMatch:
            print('Full length sequence matching {} not found!'.format(sTrimmedAlnId))
            raise
    lFullAlnIds, lFullAlnColumns = aln2columns(lFullAlnSequences)
    lTrimmedPositions = []
    iAlnPos = 0
    for lTrimmedAlnColumn in lTrimmedAlnColumns:
        while lTrimmedAlnColumn != lFullAlnColumns[iAlnPos]:
            iAlnPos += 1
        if lTrimmedAlnColumn == lFullAlnColumns[iAlnPos]:
            lTrimmedPositions.append(iAlnPos)
        iAlnPos += 1
    if len(lTrimmedPositions) != len(lTrimmedAlnColumns):
        print('Not all matching columns found in the original alignment!')
        raise
    return lTrimmedPositions

if __name__ == "__main__":
    ## pair seqs
    dNt = SeqIO(args.nt).getDict()
    lSeqs = []
    if not args.order:
        for seqPepAln in SeqIO(args.pep_aln):
            bMatch = False
            for seqNt in dNt.values():
                if seqPepAln.id in seqNt.id or seqPepAln.id in seqNt.desc:
                    lSeqs.append( [seqNt, seqPepAln] )
                    bMatch = True
                    break
            if not bMatch:
                print('No matching nucleotide sequence foud for {}!'.format(seqPepAln.id))
    else:
        for seqNt in SeqIO(args.nt):
            lSeqs.append( [seqNt] )
        iSeq = 0
        for seqPepAln in SeqIO(args.pep_aln):
            lSeqs[iSeq].append(seqPepAln)
            iSeq += 1
        for lSeqPair in lSeqs:
            if len(lSeqPair) < 2:
                print('No matching nucleotide sequence foud for {}!'.format(lSeqPair[0]))
    print('Found {} pairs.'.format(len(lSeqs)))
    fileOut = open('{}.aln.fna'.format(args.nt.rsplit('.',1)[0]), 'w')
    for (seqNt, seqPepAln) in lSeqs:
        seqNt.seq = seqNt.seq.replace('-','') ## remove dashes
        fileOut.write('>{}'.format(seqNt.id))
        if seqNt.desc != '':
            fileOut.write(' {}'.format(seqNt.desc))
        fileOut.write('\n')
        iNtPos = 0
        for iAlnPos in range(len(seqPepAln.seq)):
            if seqPepAln.seq[iAlnPos] == '-':
                fileOut.write('---')
            else:
                sFragment = seqNt.seq[iNtPos:iNtPos+3]
                #print(seqPepAln.id, iAlnPos, seqPepAln.seq[iAlnPos], iNtPos, sFragment)
                if len(sFragment) != 3:
                    raise
                fileOut.write(sFragment)
                iNtPos += 3
        fileOut.write('\n')
    fileOut.close()
    ## trimmed aln
    if args.pep_aln_trimmed != None:
        lTrimmedPositions = findTrimmedAlnPositions(args.pep_aln, args.pep_aln_trimmed)
        fileOut = open('{}.aln.trimmed.fna'.format(args.nt.rsplit('.',1)[0]), 'w')
        for (seqNt, seqPepAln) in lSeqs:
            seqNt.seq = seqNt.seq.replace('-','') ## remove dashes
            fileOut.write('>{}'.format(seqNt.id))
            if seqNt.desc != '':
                fileOut.write(' {}'.format(seqNt.desc))
            fileOut.write('\n')
            iNtPos = 0
            for iAlnPos in range(len(seqPepAln.seq)):
                if seqPepAln.seq[iAlnPos] == '-':
                    if iAlnPos in lTrimmedPositions:
                        fileOut.write('---')
                else:
                    sFragment = seqNt.seq[iNtPos:iNtPos+3]
                    #print(seqPepAln.id, iAlnPos, seqPepAln.seq[iAlnPos], iNtPos, sFragment)
                    if len(sFragment) != 3:
                        raise
                    if iAlnPos in lTrimmedPositions:
                        fileOut.write(sFragment)
                    iNtPos += 3
            fileOut.write('\n')
        fileOut.close()









