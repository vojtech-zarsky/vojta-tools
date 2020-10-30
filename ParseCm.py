import re
import math

class Hit:
 def __init__(self, sQueryId, sQueryAcc, sHitId, sHitAcc, iQueryStart, iQueryEnd, iHitStart, iHitEnd, sStrand, fGC, fE, fBits):#, sQuerySeq=None, sHitSeq=None):
  #self.db = sDb
  self.queryId = sQueryId
  self.queryAcc = sQueryAcc
  self.hitId = sHitId
  self.hitAcc = sHitAcc
  #self.identity = fIdentity
  #self.alnLen = iAlnLen
  #self.mismatches = iMismatches
  #self.gapOpenings = iGapOpenings
  self.queryStart = int(iQueryStart)
  self.queryEnd = int(iQueryEnd)
  self.hitStart = int(iHitStart)
  self.hitEnd = int(iHitEnd)
  self.strand = sStrand
  self.gc = float(fGC)
  self.eval = float(fE)
  self.score = float(fBits)
  #if sAlnSeq != None:
  # self.alnSeq = sAlnSeq
  #if sQuerySeq != None:
  # self.querySeq = sQuerySeq
  #if sHitSeq != None:
  # self.hitSeq = sHitSeq
 def show(self):
  return self.queryId, self.hitId, self.queryStart, self.queryEnd, self.hitStart, self.hitEnd, self.eval, self.score

class ParseCm:
 def __init__(self, file, cmscan, nooverlap=False):
  self.file = file
  self.sFormat = format
  self.end = 0
  self.iNooverlap = bool(nooverlap)
  self.bCmscan = bool(cmscan)
  self.lFirstLine = []
  for sLine in self.file:
   if sLine[0] == '#' or sLine.strip() == '':
    continue
   self.lFirstLine = re.split('\s+', sLine.strip(), 22)
   break
  if self.lFirstLine == []:
   self.end = 1
 def __iter__(self):
  return self
 def next(self):
  if self.end:
   raise StopIteration
  lHits = [self.lFirstLine]
  try:
   while True:
    lLine = re.split('\s+', self.file.next().strip(), 22)
    if lLine[0] == '#':
     continue
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
   (sHitId, sHitAccession, sQueryId, sQueryAccession, sMdl, sMdlFrom, sMdlTo, sSeqFrom, sSeqTo, sStrand, sTrunc, sPass, sGc, sBias, sScore, sEval, sInc, sHitDesc) = lHit[:18]
   if self.bCmscan:
    hit = Hit(sQueryId, sQueryAccession, sHitId, sHitAccession, sSeqFrom, sSeqTo, sMdlFrom, sMdlTo, sStrand, sGc, sEval, sScore)
   else:
    hit = Hit(sQueryId, sQueryAccession, sHitId, sHitAccession, sMdlFrom, sMdlTo, sSeqFrom, sSeqTo, sStrand, sGc, sEval, sScore)
   lHitsOut.append(hit)
  if self.iNooverlap:
   raise
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



