import re
import math

class Hit:
    def __init__(self, sQueryId, sHitId, sHitAcc, sHitDesc, iQueryStart, iQueryEnd, iEnvStart, iEnvEnd, iHitStart, iHitEnd, fE, fBits):#, sQuerySeq=None, sHitSeq=None):
        #self.db = sDb
        self.queryId = sQueryId
        self.hitId = sHitId
        self.hitAcc = sHitAcc
        self.hitDesc = sHitDesc
        #self.identity = fIdentity
        #self.alnLen = iAlnLen
        #self.mismatches = iMismatches
        #self.gapOpenings = iGapOpenings
        self.queryStart = int(iQueryStart)
        self.queryEnd = int(iQueryEnd)
        self.envStart = int(iEnvStart)
        self.envEnd = int(iEnvEnd)
        self.hitStart = int(iHitStart)
        self.hitEnd = int(iHitEnd)
        self.eval = float(fE)
        self.score = float(fBits)
        #if sAlnSeq != None:
        # self.alnSeq = sAlnSeq
        #if sQuerySeq != None:
        # self.querySeq = sQuerySeq
        #if sHitSeq != None:
        # self.hitSeq = sHitSeq
    def show(self):
        return self.queryId, self.hitId, self.eval, self.score, self.queryStart, self.queryEnd

class ParseHmmer:
    def __init__(self, file, nooverlap=False, filter=''):
        self.file = file
        self.sFormat = format
        self.end = 0
        self.iNooverlap = nooverlap
        self.sFilter = filter
        self.init()
    def __iter__(self):
        return self
    def __next__(self):
        if self.end:
            raise StopIteration
        return self.next()
    def init(self):
        self.lFirstLine = []
        for sLine in self.file:
            if sLine[0] == '#' or sLine.strip() == '':
                continue
            self.lFirstLine = re.split('\s+', sLine.strip(), 22)
            break
        if self.lFirstLine == []:
            self.end = 1
    def next(self):
        lHits = [self.lFirstLine]
        try:
            while True:
                sLine = next(self.file).strip()
                if sLine == '' or sLine[0] == '#':
                    continue
                lLine = re.split('\s+', sLine, 22)
                if lLine[3] != lHits[0][3]:
                    self.lFirstLine = lLine
                    break
                lHits.append(lLine)
        except:
            self.end = 1
        return self.parseDomtblout(lHits)
    
    def parseDomtblout(self, lHits):
        lHitsOut = []
        for lHit in lHits:
            if len(lHit) < 20:
                continue
            hit = Hit(lHit[3], lHit[0], lHit[1], lHit[22], lHit[17], lHit[18], lHit[19], lHit[20], lHit[15], lHit[16], lHit[12], lHit[13])
            lHitsOut.append(hit)
        
        if self.iNooverlap:
            lHitsTemp = []
            for hit in sorted(lHitsOut, key=lambda x:(x.eval, -x.score)):
                iStart1 = hit.queryStart
                iEnd1 = hit.queryEnd
                iCheck = 1
                for hit2 in lHitsTemp:
                    iStart2 = hit2.queryStart
                    iEnd2 = hit2.queryEnd
                    if (iStart1 < iStart2 and iEnd1 > iStart2) or (iStart1 >= iStart2 and iStart1 < iEnd2):
                        iCheck = 0
                        break
                if iCheck:
                    lHitsTemp.append(hit)
            lHitsOut = lHitsTemp
        return lHitsOut



