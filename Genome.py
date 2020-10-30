import re

class Seq:
    def __init__(self, sId, sDesc, sSeq):
        self.id = sId
        self.desc = ''
        if len(re.split('\s+', sDesc)) > 1:
            self.desc = re.split('\s+', sDesc, 1)[1]
        self.seq = sSeq
    def getString(self):
        if self.desc != '':
            return '>'+self.id+' '+self.desc+'\n'+self.seq+'\n'
        else:
            return '>'+self.id+'\n'+self.seq+'\n'

class Feature:
    def __init__(self, sParent, sFeatureType, bFeatureStrand, iFeatureStart, iFeatureEnd, sFeatureDesc, fFeatureScore=None):
        self.parent = sParent
        self.type = sFeatureType
        self.strand = int(bFeatureStrand)
        self.start = int(iFeatureStart)
        self.end = int(iFeatureEnd)
        self.desc = sFeatureDesc
        self.score = fFeatureScore
    def show(self):
        return self.parent, self.type, self.desc, self.strand, self.start, self.end, self.score

class Genome:
    def __init__(self, sSeqFile=None, bPrint=False):
        self.readTranslationTable()
        if sSeqFile != None:
            self.sSeqFile = sSeqFile
            self.dFeatureType2FeatureName2FeatureIndex = {}
            self.features = []
            self.dSeqId2FeatureTypes2Features = {}
            self.featureIndex = 0
            self.seqs = []
            self.dSeqId2SeqIndex = {}
            (sSeqId, sSeqDesc, sSeqSeq) = (None, '', '')
            iIndex = 0
            for sLine in open(sSeqFile):
                sLine = sLine.strip()
                if len(sLine) == 0:
                    continue
                if sLine[0] == '>':
                    if sSeqId != None:
                        self.seqs.append( Seq(sSeqId, sSeqDesc, sSeqSeq) )
                        self.dSeqId2SeqIndex[sSeqId] = iIndex
                        self.dSeqId2FeatureTypes2Features[sSeqId] = {}
                        iIndex += 1
                        (sSeqId, sSeqDesc, sSeqSeq) = (None, '', '')
                    lLine = re.split('\s+', sLine, 1)
                    sSeqId = lLine[0][1:].strip()
                    if len(lLine) > 1:
                        sSeqDesc = lLine[1]
                else:
                    sSeqSeq += sLine
            if sSeqId != None:
                self.seqs.append( Seq(sSeqId, sSeqDesc, sSeqSeq) )
                self.dSeqId2SeqIndex[sSeqId] = iIndex
                self.dSeqId2FeatureTypes2Features[sSeqId] = {}
                iIndex += 1
            if bPrint:
                global np
                global plt
                import numpy as np
                import matplotlib.pyplot as plt
    def getStats(self, bGC=False):
        #print('SeqFile: %s'%(self.sSeqFile))
        lLengths = sorted(map(lambda x:len(x.seq), self.seqs))[::-1]
        iSumLen = sum(lLengths)
        iSumCount = 0
        iN50 = None
        for iLen in lLengths:
            iSumCount += iLen
            if iN50 == None and iSumCount*2 >= iSumLen:
                iN50 = iLen
        
        #print('N seqs: %i'%(len(lLengths)))
        #print('Sum Len: %i'%(iSumLen))
        #print('N50: %i'%(iN50))
        if bGC:
            lLengths = list(map(lambda x:len(x.seq), self.seqs))
            lGCs = []
            lATs = []
            lNs = []
            for seq in self.seqs:
                iGCLoc = 0
                iATLoc = 0
                iNLoc = 0
                for sNt in seq.seq.upper():
                    if sNt in ['G','C']:
                        iGCLoc += 1
                    elif sNt in ['A','T']:
                        iATLoc += 1
                    else:
                        iNLoc += 1
                lGCs.append(iGCLoc)
                lATs.append(iATLoc)
                lNs.append(iNLoc)
            fGC = float(sum(lGCs))/(sum(lGCs)+sum(lATs))
            lGCFractions = []
            for i in range(len(lLengths)):
                fGCLoc = float(lGCs[i])/(lGCs[i]+lATs[i])
                lGCFractions.append(fGCLoc)
            #print('Ns: %i / %f'%(sum(lNs), sum(lNs)/float(iSumLen)))
            #print('GC: %f'%(fGC))
            return (len(lLengths), iSumLen, iN50, lLengths, lGCFractions, fGC, sum(lNs))
        return (len(lLengths), iSumLen, iN50, lLengths, [], -1, -1)
    def getSeq(self, sSeqId):
        return self.seqs[ self.dSeqId2SeqIndex[sSeqId.strip()] ]
    def addFeature(self, sSeqId, sFeatureType, bFeatureStrand, iFeatureStart, iFeatureEnd, sFeatureDesc, fFeatureScore=None):
        self.features.append( Feature(sSeqId, sFeatureType, int(bFeatureStrand), int(iFeatureStart), int(iFeatureEnd), sFeatureDesc, fFeatureScore) )
        if sFeatureType not in self.dSeqId2FeatureTypes2Features[sSeqId]: ##
            for (sSeqId, dFeatureTypes2Features) in self.dSeqId2FeatureTypes2Features.items():
                self.dSeqId2FeatureTypes2Features[sSeqId][sFeatureType] = []
        self.dSeqId2FeatureTypes2Features[sSeqId][sFeatureType].append( self.featureIndex )
        if sFeatureType not in self.dFeatureType2FeatureName2FeatureIndex: ##
            self.dFeatureType2FeatureName2FeatureIndex[sFeatureType] = {}
        if sFeatureDesc != None:
            if sFeatureDesc not in self.dFeatureType2FeatureName2FeatureIndex[sFeatureType]:
                self.dFeatureType2FeatureName2FeatureIndex[sFeatureType][sFeatureDesc] = [self.featureIndex]
            else:
                self.dFeatureType2FeatureName2FeatureIndex[sFeatureType][sFeatureDesc].append( self.featureIndex )
        self.featureIndex += 1
        return 0
    def getFeatures(self, sFeatureType, sFeatureDesc):
        return map( lambda x:self.features[x], self.dFeatureType2FeatureName2FeatureIndex[sFeatureType][sFeatureDesc] )
    def readGff(self, sGffFile, sFeatureTypePattern, sFeatureType, sParentFeatureType, parseFunction):
        print('read gff', sGffFile, sFeatureTypePattern, sFeatureType, sParentFeatureType)
        self.dFeatureType2FeatureName2FeatureIndex[sFeatureType] = {}
        self.dFeatureType2FeatureName2FeatureIndex[sParentFeatureType] = {}
        for (sSeqId, dFeatureTypes2Features) in self.dSeqId2FeatureTypes2Features.items():
            dFeatureTypes2Features[sFeatureType] = []
            dFeatureTypes2Features[sParentFeatureType] = []
        iCountFeatures = 0
        dParent2Features = {}
        for sLine in open(sGffFile):
            sLine = sLine.strip()
            if len(sLine) == 0 or sLine[0] == '#':
                continue
            sSeqId, sSource, sType, sStart, sEnd, sScore, sStrand, sPhase, sAttributes = sLine.split('\t')
            if sType != sFeatureTypePattern:
                continue
            iStart, iEnd = int(sStart)-1, int(sEnd)-1
            fScore = None
            if sScore != '.':
                fScore = float(sScore)
            bStrand = 0
            if sStrand == '-':
                bStrand = 1
            sParent = parseFunction(sAttributes)
            if sParent not in dParent2Features:
                dParent2Features[sParent] = []
            dParent2Features[sParent].append( (sSeqId, sFeatureType, bStrand, iStart, iEnd, fScore) )
            iCountFeatures += 1
        print('iCountFeatures', iCountFeatures)
        lFeaturesPerParent = []
        for (sParent, lFeatures) in dParent2Features.items():
            iFeatures = len(lFeatures)
            if iFeatures >= len(lFeaturesPerParent):
                lFeaturesPerParent += [0]*(iFeatures-len(lFeaturesPerParent)+1)
            lFeaturesPerParent[iFeatures] += 1
        print('lFeaturesPerParent', lFeaturesPerParent)
        ## write parent ranges
        for (sParent, lFeatures) in dParent2Features.items():
            (sSeqId, sFeatureType, bStrand, iStart, iEnd, fScore) = lFeatures[0]
            iParentStart = int(min(map(lambda x:x[3], lFeatures)))
            iParentEnd = int(max(map(lambda x:x[4], lFeatures)))
            fAvgScore = None
            fSumScore = 0
            for fScore in map(lambda x:x[5], lFeatures):
                if fScore != None:
                    fSumScore += fScore
            if fSumScore > 0:
                fAvgScore = fSumScore/float(len(lFeatures))
            #print sSeqId, sParentFeatureType, bStrand, iParentStart, iParentEnd, sParent
            self.addFeature(sSeqId, sParentFeatureType, bStrand, iParentStart, iParentEnd, sParent)
            ## write features
            for (sSeqId, sFeatureType, bStrand, iStart, iEnd, fScore) in lFeatures:
                self.addFeature(sSeqId, sFeatureType, bStrand, iStart, iEnd, sParent, fScore)
    def sortFeatures(self):
        ## sort features within sequences
        for (sSeqId, dFeatureTypes2Features) in self.dSeqId2FeatureTypes2Features.items():
            #if sSeqId != 'm51_s00275': #
            # continue
            #print 'sSeqId', sSeqId
            for (sFeatureType, lFeatures) in dFeatureTypes2Features.items():
                #if sFeatureType != 'CDS': #
                # continue
                #print 'feature', sFeatureType, len(lFeatures)
                #dFeatureTypes2Features[sFeatureType] = sorted(lFeatures, key=lambda x:self.features[x].start)
                dFeatureTypes2Features[sFeatureType] = []
                dFeatureTypes2Features[sFeatureType].append( sorted(lFeatures, key=lambda x:self.features[x].start) ) ## sort by start
                dFeatureTypes2Features[sFeatureType].append( sorted(lFeatures, key=lambda x:self.features[x].end) ) ## sort by end
        #for feature in map(lambda x:self.features[x], self.dSeqId2FeatureTypes2Features['m51_s00275']['exon'][1]):
        # print feature.show()
        return 0
    def getFeaturesInRange(self, sSeqId, iStart, iEnd=None, sQueryFeatureType = None):
        iSeqLen = len(self.getSeq(sSeqId).seq)
        if iEnd == None:
            iEnd = iStart
        [iStart, iEnd] = sorted([iStart, iEnd])
        iStart = max(0, iStart)
        iEnd = min(iSeqLen-1, iEnd)
        if iStart > iEnd:
            return []
        lReturn = []
        for (sFeatureType, (lStartSortFeatures, lEndSortFeatures)) in self.dSeqId2FeatureTypes2Features[sSeqId].items():
            if sQueryFeatureType != None and sFeatureType != sQueryFeatureType:
                continue
            dFeaturesSelect = dict( zip( lStartSortFeatures, [None]*len(lStartSortFeatures) ) )
            #print sFeatureType, len(dFeaturesSelect)
            # which features dont overlap on the left? based on ends
            for iFeature in lEndSortFeatures:
                feature = self.features[iFeature]
                if feature.end < iStart:
                    del dFeaturesSelect[iFeature]
                else:
                    break
            #print 'left removal', len(dFeaturesSelect)
            # which features dont overlap on the right? based on starts
            for iFeature in lStartSortFeatures[::-1]:
                feature = self.features[iFeature]
                if feature.start > iEnd:
                    del dFeaturesSelect[iFeature]
                else:
                    break
            lReturn += dFeaturesSelect.keys()
        return map(lambda x:self.features[x], lReturn)
    def iterateFeatures(self, sQueryFeatureType=None):
        for (sFeatureType, dFeatureName2FeatureIndex) in self.dFeatureType2FeatureName2FeatureIndex.items():
            if sQueryFeatureType != None and sQueryFeatureType != sFeatureType:
                continue
            for (sFeatureName, lFeatureIndices) in dFeatureName2FeatureIndex.items():
                for iFeature in lFeatureIndices:
                    yield self.features[iFeature]
    def reverseComplement(self, sSeq):
        sSeqOut = ''
        for sAa in sSeq.upper()[::-1]:
            if sAa == 'C':
                sSeqOut += 'G'
            elif sAa == 'G':
                sSeqOut += 'C'
            elif sAa == 'A':
                sSeqOut += 'T'
            elif sAa == 'T':
                sSeqOut += 'A'
            elif sAa == 'N':
                sSeqOut += 'N'
            else:
                raise
        return sSeqOut
    def getFeatureSeq(self, feature):
        sSeq = self.getSeq(feature.parent).seq[feature.start:feature.end+1]
        if feature.strand == 1:
            sSeq = self.reverseComplement(sSeq)
        return sSeq
    def readTranslationTable(self):
        sAAs = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
        sSts = '---M------**--*----M---------------M----------------------------'
        sBs1 = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
        sBs2 = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
        sBs3 = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'
        self.dCodon2Aa = {}
        for iIndex in range(len(sAAs)):
            self.dCodon2Aa[ sBs1[iIndex]+sBs2[iIndex]+sBs3[iIndex] ] = sAAs[iIndex]
    def translate(self, sNtSeq):
        sNtSeq = sNtSeq.upper()
        sAaSeq = ''
        for iIndex in range(0, len(sNtSeq)-2, 3):
            if sNtSeq[iIndex:iIndex+3] in self.dCodon2Aa:
                sAaSeq += self.dCodon2Aa[ sNtSeq[iIndex:iIndex+3] ]
            else:
                sAaSeq += 'X'
                print(sNtSeq[iIndex:iIndex+3])
        return sAaSeq
    def translateFeatures(self, lFeatures):
        bStrand = lFeatures[0].strand
        sNtSeq = ''
        if bStrand == 0:
            lFeatures = sorted(lFeatures, key=lambda x:x.start)
        elif bStrand == 1:
            lFeatures = sorted(lFeatures, key=lambda x:-x.start)
        for feature in lFeatures:
            sNtSeq += self.getFeatureSeq(feature)
        return self.translate(sNtSeq)
    def addProteinFeature(self, sGeneId, iStart, iEnd, sDesc, fE):
        lExons = self.getFeatures('CDS', sGeneId)
        sSeqId = lExons[0].parent
        bStrand = lExons[0].strand
        if bStrand == 1:
            lExons = lExons[::-1]
        iAaIndex = 0
        iNtIndex = 0
        #iGenomeIndex = 0
        #if bStrand == 1:
        # iNtIndex
        iGenomicPosition = None
        for exon in lExons:
            #if 
            #print 'exon', exon.start, exon.end
            lExonRange = range(exon.start, exon.end+1)
            if bStrand == 1:
                lExonRange = lExonRange[::-1]
            iGenomicStart = None
            iGenomicEnd = None
            for iGenomicIndex in lExonRange:
                if iAaIndex >= iStart and iAaIndex <= iEnd:
                    if iGenomicStart == None:
                        iGenomicStart = iGenomicIndex
                    iGenomicEnd = iGenomicIndex
                iNtIndex += 1
                if iNtIndex%3 == 0:
                    iAaIndex += 1
            #print 'feature', iGenomicStart, iGenomicEnd
            if iGenomicStart != None:
                self.addFeature(sSeqId, 'eggnog', bStrand, iGenomicStart, iGenomicEnd, sDesc, fE)
        return 0
    def myplotRectangle(self, ax, tRange, fWidth=1., tOrigin=(0, 0), sColor='0.5', fOutlineWidth=0):
        tRange = (tRange[0], tRange[1]+1) ## to account for my position notation
        ax.fill( [tOrigin[0]+tRange[0], tOrigin[0]+tRange[1], tOrigin[0]+tRange[1], tOrigin[0]+tRange[0]], [tOrigin[1]+fWidth, tOrigin[1]+fWidth, tOrigin[1]-fWidth, tOrigin[1]-fWidth], color=sColor )
        if fOutlineWidth > 0:
            ax.plot( [tOrigin[0]+tRange[0], tOrigin[0]+tRange[1], tOrigin[0]+tRange[1], tOrigin[0]+tRange[0], tOrigin[0]+tRange[0]], [tOrigin[1]+fWidth, tOrigin[1]+fWidth, tOrigin[1]-fWidth, tOrigin[1]-fWidth, tOrigin[1]+fWidth], color='k' )
        return 0
    def showRegion(self, sSeqId, iStart, iEnd):
        fig,ax = plt.subplots()
        ax.set_ylim(-4, 4)
        iStartRelative = 0
        iEndRelative = abs(iStart-iEnd)+1
        self.myplotRectangle(ax, (iStartRelative, iEndRelative), fWidth=.5, tOrigin=(0, 0), sColor='0.5', fOutlineWidth=0.)
        for cds in self.getFeaturesInRange(sSeqId, iStart, iEnd, 'CDS'):
            self.myplotRectangle(ax, (cds.start-iStart, cds.end-iStart), fWidth=.5, tOrigin=(0, 1), sColor='r', fOutlineWidth=0)
            #if feature.type == 'CDS':
            # self.myplotRectangle(ax, (feature.start-iStart, feature.end-iStart), fWidth=.5, tOrigin=(0, 1), sColor='r', fOutlineWidth=0)
            #elif feature.type == 'exon':
            # self.myplotRectangle(ax, (feature.start-iStart, feature.end-iStart), fWidth=.5, tOrigin=(0, -1), sColor='b', fOutlineWidth=0)
            #elif feature.type == 'eggnog':
            # self.myplotRectangle(ax, (feature.start-iStart, feature.end-iStart), fWidth=.5, tOrigin=(0, -2), sColor='g', fOutlineWidth=0)
        plt.show()
        return 0


"""
genome = Genome('../sequences/mastiga_genome_v5.1.fasta')
print 'len(genome.seqs)', len(genome.seqs)
genome.readGff( '../sequences/Masba_all.gff3', 'CDS', 'CDS', 'gene', lambda x:x.split('Parent=')[1].split(';')[0].rsplit('.',1)[0] )
genome.readGff( '../transcriptomeMapping/mastiga_cDNA_v3-all_variants.gff3_gene', 'exon', 'exon', 'transcript', lambda x:x.split('Name=')[1].split(';')[0] )
genome.sortFeatures()

for feature in genome.iterateFeatures():
    print feature.show()
"""

#genome.getFeaturesInRange('m51_s00275', 10000000, -20)
#genome.getFeaturesInRange('m51_s00275', 12107)







