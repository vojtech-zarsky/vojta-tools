import math
import re
import copy
import numpy as np

class Hit:
 def __init__(self, sDb, sQueryId, sHitId, fIdentity, iAlnLen, iMismatches, iGapOpenings, iQueryStart, iQueryEnd, iHitStart, iHitEnd, fE, fBits, sAlnSeq=None):#, sQuerySeq=None, sHitSeq=None):
  self.db = sDb
  self.queryId = sQueryId
  self.hitId = sHitId
  self.identity = fIdentity
  self.alnLen = iAlnLen
  self.mismatches = iMismatches
  self.gapOpenings = iGapOpenings
  self.queryStart = iQueryStart
  self.queryEnd = iQueryEnd
  self.hitStart = iHitStart
  self.hitEnd = iHitEnd
  self.eval = fE
  self.score = fBits
  if sAlnSeq != None:
   self.alnSeq = sAlnSeq
  #if sQuerySeq != None:
  # self.querySeq = sQuerySeq
  #if sHitSeq != None:
  # self.hitSeq = sHitSeq
 def show(self):
  return self.queryId, self.hitId, self.eval, self.score, self.queryStart, self.queryEnd

class ParseBlast:
 def __init__(self, file, format='table', ncbi=False):
  self.file = file
  self.sFormat = format
  self.end = 0
  self.bNcbi = ncbi
  if self.sFormat == 'table':
   self.initTable()
  elif self.sFormat == 'mytable':
   self.initMytable()
  elif self.sFormat == 'xml':
   self.initXml()
  elif self.sFormat == 'hspsxml':
   self.initXml()
  else:
   raise
 def __iter__(self):
  return self
 def __next__(self):
  return self.next()
 
 def initTable(self):
  self.lFirstLine = []
  for sLine in self.file:
   if sLine.strip() == '':
    continue
   self.lFirstLine = re.split('\s+', sLine.strip())
   break
  if self.lFirstLine == []:
   self.end = 1
 
 def initMytable(self):
  #self.dAminoAcid2Index = {}
  #sAminoAcids = 'ARNDCUQEGHILKMFPSTWYVBZX*-'
  #for iIndex in range(len(sAminoAcids)):
  # self.dAminoAcid2Index[sAminoAcids[iIndex]] = iIndex
  self.lFirstLine = []
  for sLine in self.file:
   if sLine.strip() == '':
    continue
   self.lFirstLine = re.split('\s+', sLine.strip())
   break
  if self.lFirstLine == []:
   self.end = 1
 
 def initXml(self):
  self.sFirstLine = ''
  for sLine in self.file:
   if re.search('<Iteration>', sLine):
    self.sFirstLine = sLine
    break
  if self.sFirstLine == '':
   self.end = 1
 
 def next(self):
  if self.end:
   raise StopIteration()
  if self.sFormat == 'table':
   return self.nextTable()
  elif self.sFormat == 'mytable':
   return self.nextMytable()
  elif self.sFormat == 'xml':
   return self.nextXml()
  elif self.sFormat == 'hspsxml':
   return self.nextHspsXml()
  else:
   raise
 
 def nextHsps(self):
  if self.sFormat == 'xml':
   return self.nextHspsXml()
  else:
   raise
 
 def nextTable(self):
  lHits = [self.lFirstLine]
  try:
   while True:
    lLine = re.split('\s+', self.file.next().strip())
    if lLine[0] != lHits[0][0]:
     self.lFirstLine = lLine
     break
    lHits.append(lLine)
  except:
   self.end = 1
  return self.parseTable(lHits)
 
 def nextMytable(self):
  lHits = [self.lFirstLine]
  try:
   while True:
    lLine = re.split('\s+', self.file.next().strip())
    if lLine == [''] or lLine == ['Search', 'has', 'CONVERGED!']:
     continue
    if lLine[0] != lHits[0][0]:
     self.lFirstLine = lLine
     break
    lHits.append(lLine)
  except:
   self.end = 1
  return self.parseMytable(lHits)
 
 def nextXml(self):
  sHits = self.sFirstLine
  try:
   while True:
    sLine = self.file.next()
    if re.search('<Iteration>', sLine):
     self.sFirstLine = sLine
     break
    sHits = sHits+sLine
  except:
   self.end = 1
  return self.parseXml(sHits)
 
 def nextHspsXml(self):
  sHits = self.sFirstLine
  try:
   while True:
    sLine = self.file.next()
    if re.search('<Iteration>', sLine):
     self.sFirstLine = sLine
     break
    sHits = sHits+sLine
  except:
   self.end = 1
  return self.parseHspsXml(sHits)
 
 def parseTable(self, lHits):
  if len(lHits) == 0 or len(lHits[0]) == 0:
   return []
  lHitsOut = []
  for lHit in lHits:
   #dHit = {}
   #dHit['queryId'] = lHit[0]
   #dHit['hitId'] = lHit[1]
   #dHit['identity'] = float(lHit[2])
   #dHit['alnLen'] = int(lHit[3])
   #dHit['mismatches'] = int(lHit[4])
   #dHit['gapOpenings'] = int(lHit[5])
   #dHit['queryStart'] = int(lHit[6])
   #dHit['queryEnd'] = int(lHit[7])
   #dHit['hitStart'] = int(lHit[8])
   #dHit['hitEnd'] = int(lHit[9])
   #dHit['e'] = float(lHit[10])
   #dHit['logE'] = self.getE(dHit['e'])
   #dHit['bits'] = float(lHit[11])
   lHitsOut.append(Hit( self.file.name, lHit[0], lHit[1], float(lHit[2]), int(lHit[3]), int(lHit[4]), int(lHit[5]), int(lHit[6]), int(lHit[7]), int(lHit[8]), int(lHit[9]), float(lHit[10]), float(lHit[11]) ))
  return lHitsOut
 
 def seq2array(self, sSeq):
  try:
   arraySeq = np.array( map(lambda x:self.dAminoAcid2Index[x],sSeq.upper()) )
  except:
   lTemp = []
   for sAa in sSeq.upper():
    if sAa in self.dAminoAcid2Index:
     lTemp.append(self.dAminoAcid2Index[sAa])
    else:
     lTemp.append(self.dAminoAcid2Index['X'])
   arraySeq = np.array(lTemp)
  return arraySeq
 
 def parseMytable(self, lHits):
  if len(lHits) == 0 or len(lHits[0]) == 0:
   return []
  lHitsOut = []
  for lHit in lHits:
   #print lHit
   #dHit = {}
   #dHit['queryId'] = lHit[0]
   #dHit['hitId'] = lHit[1]
   #dHit['identity'] = float(lHit[2])
   #dHit['alnLen'] = int(lHit[3])
   #dHit['mismatches'] = int(lHit[4])
   #dHit['gapOpenings'] = int(lHit[5])
   #dHit['queryStart'] = int(lHit[6])
   #dHit['queryEnd'] = int(lHit[7])
   #dHit['hitStart'] = int(lHit[8])
   #dHit['hitEnd'] = int(lHit[9])
   #dHit['e'] = float(lHit[10])
   #dHit['logE'] = self.getE(dHit['e'])
   #dHit['bits'] = float(lHit[11])
   #dHit['querySeq'] = lHit[12]
   #dHit['hitSeq'] = lHit[13]
   lHitsOut.append(Hit( self.file.name, lHit[0], lHit[1], float(lHit[2]), int(lHit[3]), int(lHit[4]), int(lHit[5]), int(lHit[6]), int(lHit[7]), int(lHit[8]), int(lHit[9]), float(lHit[10]), float(lHit[11]), lHit[12] ))#, lHit[13] ))
   #lHitsOut.append(dHit)
  return lHitsOut
  
 def traceBack(self, hit, sQuerySeq):
  iQueryIndex = hit.queryStart-1
  sTraceBack = hit.alnSeq
  iTraceBackLen = len(sTraceBack)
  iIndex = 0
  sQueryAlnSeq = ''
  sHitAlnSeq  = ''
  while iIndex < iTraceBackLen:
   if sTraceBack[iIndex].isdigit(): ## number of identities
    iIndex2 = iIndex+1
    for iIndex2 in range(iIndex+2, iTraceBackLen+1):
     if not sTraceBack[iIndex:iIndex2].isdigit():
      iIndex2 = iIndex2-1
      break 
    iIdentityLen = int(sTraceBack[iIndex:iIndex2])
    sHitAlnSeq += sQuerySeq[iQueryIndex:iQueryIndex+iIdentityLen]
    sQueryAlnSeq += sQuerySeq[iQueryIndex:iQueryIndex+iIdentityLen]
    iQueryIndex += iIdentityLen
    iIndex = iIndex2 
   elif sTraceBack[iIndex] == 'X': ##
    sHitAlnSeq += sTraceBack[iIndex+1].lower()
    sQueryAlnSeq += sQuerySeq[iQueryIndex].lower()
    iIndex += 2
    iQueryIndex += 1
   elif sTraceBack[iIndex] == '-': ## query/hit alignment
    sHitAlnSeq += sTraceBack[iIndex+1]
    sQueryAlnSeq += '-'
    iIndex += 2
   else:
    sHitAlnSeq += sTraceBack[iIndex+1]
    sQueryAlnSeq += sTraceBack[iIndex]
    iIndex += 2
    iQueryIndex += 1
  hit.querySeq = sQueryAlnSeq
  hit.hitSeq = sHitAlnSeq
  if len(sQueryAlnSeq) != len(sHitAlnSeq):
   return 0
  return hit
  
 
 def getE(self, fE, fMAX=200.0):
  fE = float(fE)
  fMIN = 0.0
  if fE == 0:
   return fMAX
  fE = -math.log(fE, 10)
  if fE > fMAX:
   fE = fMAX
  elif fE < fMIN:
   fE = fMIN
  return fE
 
 def parseXml(self, sHits):
  lHitsOut = []
  dHitMain = {}
  dHitMain['queryDesc'] = sHits.split('<Iteration_query-def>')[1].split('</')[0].strip().split(' ',1)
  if len(dHitMain['queryDesc']) > 1:
   dHitMain['queryDesc'] = dHitMain['queryDesc'][1]
  else:
   dHitMain['queryDesc'] = ''
  dHitMain['queryId'] = str( sHits.split('<Iteration_query-def>')[1].split('</')[0].strip() ).split()[0].strip()
  dHitMain['queryLen'] = int( sHits.split('<Iteration_query-len>')[1].split('</')[0].strip() )
  for sHit in sHits.split('<Hit>')[1:]:
   dHit = copy.deepcopy(dHitMain)
   dHit['hitDesc'] = str( sHit.split('<Hit_def>')[1].split('</')[0].strip() )
   dHit['hitId'] = str( sHit.split('<Hit_def>')[1].split('</')[0].strip() ).split()[0]
   if self.bNcbi:
    dHit['hitId'] = sHit.split('<Hit_id>')[1].split('</')[0].strip()
   dHit['hitLen'] = int( sHit.split('<Hit_len>')[1].split('</')[0].strip() )
   dHit['hsps'] = []
   for sHsp in sHit.split('<Hsp>')[1:]:
    dHsp = copy.deepcopy(dHit)
    del dHsp['hsps']
    dHsp['bitScore'] = float( sHsp.split('<Hsp_bit-score>')[1].split('</')[0].strip() )
    dHsp['score'] = float( sHsp.split('<Hsp_score>')[1].split('</')[0].strip() )
    dHsp['queryStart'] = int( sHsp.split('<Hsp_query-from>')[1].split('</')[0].strip() )
    dHsp['queryEnd'] = int( sHsp.split('<Hsp_query-to>')[1].split('</')[0].strip() )
    dHsp['hitStart'] = int( sHsp.split('<Hsp_hit-from>')[1].split('</')[0].strip() )
    dHsp['hitEnd'] = int( sHsp.split('<Hsp_hit-to>')[1].split('</')[0].strip() )
    dHsp['identity'] = int( sHsp.split('<Hsp_identity>')[1].split('</')[0].strip() )
    dHsp['positive'] = int( sHsp.split('<Hsp_positive>')[1].split('</')[0].strip() )
    dHsp['gaps'] = 0
    if len(sHsp.split('<Hsp_gaps>')) > 1:
     dHsp['gaps'] = int( sHsp.split('<Hsp_gaps>')[1].split('</')[0].strip() )
    dHsp['alnLen'] = int( sHsp.split('<Hsp_align-len>')[1].split('</')[0].strip() )
    dHsp['querySeq'] = str( sHsp.split('<Hsp_qseq>')[1].split('</')[0].strip() )
    dHsp['hitSeq'] = str( sHsp.split('<Hsp_hseq>')[1].split('</')[0].strip() )
    dHsp['midlineSeq'] = str( sHsp.split('<Hsp_midline>')[1].split('</')[0].strip() )
    dHsp['e'] = float( sHsp.split('<Hsp_evalue>')[1].split('</')[0].strip() )
    dHsp['logE'] = self.getE( dHsp['e'] )
    dHit['hsps'].append(dHsp)
   lHitsOut.append(dHit)
  return lHitsOut
 
 def parseHspsXml(self, sHits):
  lHspsOut = []
  dHitMain = {}
  dHitMain['queryDesc'] = sHits.split('<Iteration_query-def>')[1].split('</')[0].strip().split(' ',1)
  if len(dHitMain['queryDesc']) > 1:
   dHitMain['queryDesc'] = dHitMain['queryDesc'][1]
  else:
   dHitMain['queryDesc'] = ''
  dHitMain['queryId'] = str( sHits.split('<Iteration_query-def>')[1].split('</')[0].strip() ).split()[0].strip()
  dHitMain['queryLen'] = int( sHits.split('<Iteration_query-len>')[1].split('</')[0].strip() )
  for sHit in sHits.split('<Hit>')[1:]:
   dHit = copy.deepcopy(dHitMain)
   dHit['hitDesc'] = str( sHit.split('<Hit_def>')[1].split('</')[0].strip() )
   dHit['hitId'] = str( sHit.split('<Hit_id>')[1].split('</')[0].strip() ).split()[0]
   dHit['hitLen'] = int( sHit.split('<Hit_len>')[1].split('</')[0].strip() )
   #dHit['hsps'] = []
   for sHsp in sHit.split('<Hsp>')[1:]:
    dHsp = copy.deepcopy(dHit)
    #del dHsp['hsps']
    dHsp['bitScore'] = float( sHsp.split('<Hsp_bit-score>')[1].split('</')[0].strip() )
    dHsp['score'] = int( sHsp.split('<Hsp_score>')[1].split('</')[0].strip() )
    dHsp['queryStart'] = int( sHsp.split('<Hsp_query-from>')[1].split('</')[0].strip() )
    dHsp['queryEnd'] = int( sHsp.split('<Hsp_query-to>')[1].split('</')[0].strip() )
    dHsp['hitStart'] = int( sHsp.split('<Hsp_hit-from>')[1].split('</')[0].strip() )
    dHsp['hitEnd'] = int( sHsp.split('<Hsp_hit-to>')[1].split('</')[0].strip() )
    dHsp['identity'] = int( sHsp.split('<Hsp_identity>')[1].split('</')[0].strip() )
    dHsp['positive'] = int( sHsp.split('<Hsp_positive>')[1].split('</')[0].strip() )
    dHsp['gaps'] = 0
    if len(sHsp.split('<Hsp_gaps>')) > 1:
     dHsp['gaps'] = int( sHsp.split('<Hsp_gaps>')[1].split('</')[0].strip() )
    dHsp['alnLen'] = int( sHsp.split('<Hsp_align-len>')[1].split('</')[0].strip() )
    dHsp['querySeq'] = str( sHsp.split('<Hsp_qseq>')[1].split('</')[0].strip() )
    dHsp['hitSeq'] = str( sHsp.split('<Hsp_hseq>')[1].split('</')[0].strip() )
    dHsp['midlineSeq'] = str( sHsp.split('<Hsp_midline>')[1].split('</')[0].strip() )
    dHsp['e'] = float( sHsp.split('<Hsp_evalue>')[1].split('</')[0].strip() )
    dHsp['logE'] = self.getE( dHsp['e'] )
    lHspsOut.append(dHsp)
  return lHspsOut


 
