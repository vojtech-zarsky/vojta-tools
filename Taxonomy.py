import re

class Taxonomy:
 def __init__(self, sTaxdumpDir='./taxdump/', bNewParse=True, bCorrections=True):
  if bNewParse:
   self.newParse(sTaxdumpDir)
  self.dCorrections = {'yersinia':'629', 'morganella':'581', 'bacteria':'2'}
 
 def newParse(self, sTaxdumpDir):
  self.dNamesIds = {}
  self.dIdsNames = {}
  print 'reading names/ids'
  lNonscientificLines = []
  for sLine in open(sTaxdumpDir+'/names.dmp'):
   lLine = sLine.replace('\'-',' -').replace('\'','').lower().strip('	|\n').split('\t|\t') ##
   if lLine[3] == 'scientific name':
    if lLine[1] not in self.dNamesIds:
     self.dNamesIds[lLine[1]] = []
    self.dNamesIds[lLine[1]].append(lLine[0])
   else:
    lNonscientificLines.append(lLine)
   if lLine[0] not in self.dIdsNames:
    self.dIdsNames[lLine[0]] = []
   if lLine[3] == 'scientific name':
    self.dIdsNames[lLine[0]] = [lLine[1]]+self.dIdsNames[lLine[0]]
   else:
    self.dIdsNames[lLine[0]].append(lLine[1])
  
  for lLine in lNonscientificLines:
   if lLine[1] not in self.dNamesIds:
    self.dNamesIds[lLine[1]] = [lLine[0]]
  
  print 'reading nodes'
  self.dNodes = {}
  for sLine in open(sTaxdumpDir+'/nodes.dmp'):
   lLine = sLine.lower().strip('	|\n').split('\t|\t') ##
   self.dNodes[lLine[0]] = {'parent':lLine[1], 'rank':lLine[2]}
  
  print 'reading merges'
  self.dMerges = {}
  for sLine in open(sTaxdumpDir+'/merged.dmp'):
   lLine = map(lambda x:x.strip(), sLine.split('|',2)[:2])
   self.dMerges[lLine[0]] = lLine[1]
 
 def getTaxonomy(self, sTaxName):
  lTaxonomy = []
  sTaxName = sTaxName.replace('\'-',' -').replace('\'','').lower().strip()
  if sTaxName == None or sTaxName == '':
   return lTaxonomy
  if sTaxName not in self.dNamesIds:
   sTaxName = re.sub('^uncultured','',sTaxName)
   sTaxName = re.sub('^candidatus','',sTaxName)
   sTaxName = sTaxName.strip()
   if sTaxName not in self.dNamesIds:
    sTaxName = sTaxName.split()[0]
    if sTaxName not in self.dNamesIds:
     return lTaxonomy
  lTaxIds = self.dNamesIds[sTaxName]
  sTaxId = lTaxIds[0]
  if len(lTaxIds) > 1 and sTaxName in self.dCorrections:
   sTaxId = self.dCorrections[sTaxName]
  
  while sTaxId != '1':
   dNode = self.dNodes[sTaxId]
   lNames = self.dIdsNames[sTaxId]
   lTaxonomy.append((sTaxId, lNames[0], dNode['rank']))
   sTaxId = dNode['parent']
  return lTaxonomy
 
 def getTaxonomyByTaxId(self, sTaxId):
  if sTaxId in self.dMerges:
   sTaxId = self.dMerges[sTaxId]
  lTaxonomy = []
  while sTaxId != '1':
   try:
    dNode = self.dNodes[sTaxId]
   except:
    return lTaxonomy
   lNames = self.dIdsNames[sTaxId]
   lTaxonomy.append((sTaxId, lNames[0], dNode['rank']))
   sTaxId = dNode['parent']
  return lTaxonomy
 
 def isTaxon(self, sTaxName):
  if sTaxName in self.dNamesIds:
   return True
  return False


