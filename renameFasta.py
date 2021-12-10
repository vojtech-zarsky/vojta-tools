## python ../myclasses/renameFasta.py --infile temp.fasta --taxdump ../myclasses/taxdump --taxlevel 4

import sys
import os
import argparse

## for python3

#bgzip

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


class Taxonomy:
    def __init__(self, sTaxdumpDir='./taxdump/', sTaxdumpMod=None, bCorrections=True, bNewParse=True, bAddSupergroups=False, bCleanTaxons=False):
        self.dCorrections = {'yersinia':'629', 'morganella':'581', 'bacteria':'2'}
        self.lSupergroups = [('eukaryota','archaeplastida',['viridiplantae','rhodelphis','rhodophyta','glaucocystophyceae'])]
        self.bAddSupergroups = bAddSupergroups
        self.bCleanTaxons = bCleanTaxons
        if bNewParse:
            self.newParse(sTaxdumpDir, sTaxdumpMod)
    
    def newParse(self, sTaxdumpDir, sTaxdumpMod):
        ## read mod file
        lMods = []
        if sTaxdumpMod != None:
            for sLine in open(sTaxdumpMod):
                if sLine.strip() == '' or sLine[0] == '#':
                    continue
                lLine = sLine.strip().split('\t')
                lMods.append( lLine )
                #if lLine[0] == 'add':
                #    lAdd.append( lLine[1:] )
        self.dNamesIds = {}
        self.dIdsNames = {}
        print('reading names/ids')
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
        print('reading nodes')
        self.dNodes = {}
        for sLine in open(sTaxdumpDir+'/nodes.dmp'):
            lLine = sLine.lower().strip('	|\n').split('\t|\t') ##
            self.dNodes[lLine[0]] = {'parent':lLine[1], 'rank':lLine[2]}
        
        print('reading merges')
        self.dMerges = {}
        for sLine in open(sTaxdumpDir+'/merged.dmp'):
            lLine = list(map(lambda x:x.strip(), sLine.split('|',2)[:2]))
            self.dMerges[lLine[0]] = lLine[1]
        iNewTaxId = 1
        ## add from mod
        for lMod in lMods:
            if lMod[0] == 'add':
                #self.dNamesIds = {}
                #self.dIdsNames = {}
                #self.dNodes
                lTaxonomy = []
                for sTemp in lMod[1:]:
                    sTemp = sTemp.lower().strip()
                    if sTemp.isdigit():
                        sNewTaxId = sTemp
                        sNewTaxName = 'taxid{}'.format(sTemp)
                    elif sTemp.split(' ',1)[0].isdigit():
                        sNewTaxId = sTemp.split(' ',1)[0]
                        sNewTaxName = sTemp.split(' ',1)[1]
                    elif sTemp in self.dNamesIds:
                        sNewTaxId = self.dNamesIds[sTemp][0]
                        sNewTaxName = sTemp
                    else:
                        sNewTaxId = '0{}'.format(iNewTaxId)
                        iNewTaxId += 1
                        sNewTaxName = sTemp
                    lTaxonomy.append( (sNewTaxId, sNewTaxName) )
                for iTaxLevel in range(len(lTaxonomy)):
                    (sNewTaxId, sNewTaxName) = lTaxonomy[iTaxLevel]
                    if sNewTaxId in self.dIdsNames:
                        break
                    self.dIdsNames[sNewTaxId] = [sNewTaxName]
                    self.dNamesIds[sNewTaxName] = [sNewTaxId]
                    sParent = '1'
                    if iTaxLevel+1 <  len(lTaxonomy):
                        sParent = lTaxonomy[iTaxLevel+1][0]
                    self.dNodes[sNewTaxId] = {'parent':sParent, 'rank':''}
        ## add supergroups
        if self.bAddSupergroups:
            iNewTaxId = -1
            for (sParent, sSupergroup, lChildren) in self.lSupergroups:
                ## add names/ids
                self.dIdsNames[str(iNewTaxId)] = [sSupergroup]
                print('new supergroup:', str(iNewTaxId), [sSupergroup])
                self.dNamesIds[sSupergroup] = [str(iNewTaxId)]
                self.dNodes[str(iNewTaxId)] = {'parent':self.dNamesIds[sParent][0], 'rank':''}
                for sChild in lChildren:
                    self.dNodes[self.dNamesIds[sChild][0]]['parent'] = str(iNewTaxId)
            iNewTaxId -= 1

    
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
            sName = lNames[0]
            if self.bCleanTaxons:
                sName = sName.replace(';','').replace(',','')
            lTaxonomy.append((sTaxId, sName, dNode['rank']))
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
            sName = lNames[0]
            if self.bCleanTaxons:
                sName = sName.replace(';','').replace(',','')
            lTaxonomy.append((sTaxId, sName, dNode['rank']))
            sTaxId = dNode['parent']
        return lTaxonomy
    
    def isTaxon(self, sTaxName):
        if sTaxName in self.dNamesIds:
            return True
        return False




parser = argparse.ArgumentParser(prog='Rename FASTA sequences.')
parser.add_argument('--infile', required=True, help='Input FASTA file.')
parser.add_argument('--taxlevel', required=False, default='3', help='Taxlevel [3]')
parser.add_argument('--taxdump', required=True, help='NCBI taxdump folder. (https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz)')
args = parser.parse_args()

args.taxlevel = int(args.taxlevel)

if not os.path.isfile(args.infile):
    print('Cannot find {}!'.format(args.infile))
    raise
lIdsUsed = []
taxonomy = Taxonomy(sTaxdumpDir=args.taxdump)
fileOut = open('{}.names.fasta'.format(args.infile.rsplit('.',1)[0]),'w')

for seq in SeqIO(args.infile):
    seq.id = seq.id.rsplit('.',1)[0] ##
    if seq.id in lIdsUsed:
        continue
    lIdsUsed.append(seq.id)
    lTaxonomy = []
    sOrgn = None
    sTaxon = None
    if re.search('.+\[.+]',seq.desc):
        sOrgn = seq.desc.split('[')[-1].split(']')[0].strip()
        lTaxonomy = list(map( lambda x:x[1], taxonomy.getTaxonomy(sOrgn) ))[::-1]
        if 'candidatus' in sOrgn:
            sOrgn = '_'.join(sOrgn.replace('.','').strip().split()[:3])
        else:
            sOrgn = '_'.join(sOrgn.replace('.','').strip().split()[:2])
        if len(lTaxonomy) > args.taxlevel:
            sTaxon = '_'.join(lTaxonomy[args.taxlevel].replace('.','').replace('(','').replace(')','').strip().split()[:3])
    sId = seq.id
    if sTaxon != None:
        sId = sTaxon+'.'+sId
    if sOrgn != None:
        sId = sId+'.'+sOrgn
    fileOut.write('>{} {} {}\n{}\n'.format(sId, seq.id, seq.desc, seq.seq))
fileOut.close()
