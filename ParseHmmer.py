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
    def __init__(self, file, sFormat='domtblout', nooverlap=False, filter='', iFirstNHits=None):
        self.file = file
        if isinstance(self.file, str):
            self.file = open(self.file)
        self.sFormat = sFormat
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
        if self.sFormat=='domtblout':
            self.lFirstLine = []
            for sLine in self.file:
                if sLine[0] == '#' or sLine.strip() == '':
                    continue
                self.lFirstLine = re.split('\s+', sLine.strip(), 22)
                break
            if self.lFirstLine == []:
                self.end = 1
        elif self.sFormat=='o':
            self.lHits = []
            #iCount = 0
            for sLine in self.file:
                #iCount += 1
                #if iCount%1000 == 0:
                #    print(iCount)
                self.lHits.append(sLine.strip())
                if sLine.strip() == '//':
                    break
    def next(self):
        if self.sFormat=='domtblout':
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
        elif self.sFormat=='o':
            lHits = self.lHits
            self.lHits = []
            bCheck = False
            for sLine in self.file:
                self.lHits.append(sLine.strip())
                if sLine.strip() == '//':
                    bCheck = True
                    break
            if not bCheck:
                self.end = 1
            return self.parseO(lHits)
    
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
    
    def parseO(self, lHits):
        sHits = '\n'.join(lHits)
        lHits = []
        for sHit in sHits.strip().split('>>')[1:]:
            for sHsp in sHit.split('== domain')[1:]:
                fScore = float(sHsp.split('score:',1)[1].strip().split()[0])
                fEval = float(sHsp.split('E-value:',1)[1].strip().split()[0])
                #print(sHit)
                #lStatusLine = sHit.split(' ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----')[1].strip().split('\n',1)[0].strip().split()
                #fScore = float(lStatusLine[2])
                #fEval = float(lStatusLine[5])
                #iHmmStart = int(lStatusLine[6])
                #iHmmEnd = int(lStatusLine[7])
                #iAliStart = int(lStatusLine[9])
                #iAliEnd = int(lStatusLine[10])
                iAddLineIndex = 0
                if sHit.split('==')[1].split('\n')[1].strip().split()[-1] == 'CS':
                    iAddLineIndex = 1
                lHmmLine = sHit.split('==')[1].split('\n')[1+iAddLineIndex].strip().split()
                sHmmId = lHmmLine[0]
                sHmmSeq = lHmmLine[2]
                iHmmStart = int(lHmmLine[1])
                iHmmEnd = int(lHmmLine[3])
                lSeqLine = sHit.split('==')[1].split('\n')[3+iAddLineIndex].strip().split()
                sSeqId = lSeqLine[0]
                sSeqSeq = lSeqLine[2]
                iAliStart = int(lSeqLine[1])
                iAliEnd = int(lSeqLine[3])
                lHits.append( (sHmmId, sSeqId, iHmmStart, iHmmEnd, iAliStart, iAliEnd, fEval, fScore, sHmmSeq, sSeqSeq) )
                #print((sHmmId, sSeqId, iHmmStart, iHmmEnd, iAliStart, iAliEnd, fEval, fScore, sHmmSeq, sSeqSeq))
        return lHits

